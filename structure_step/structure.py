# -*- coding: utf-8 -*-

"""Non-graphical part of the Structure step in a SEAMM flowchart"""

import logging
from pathlib import Path
import pkg_resources
import pprint  # noqa: F401
import sys
import time

from tabulate import tabulate
from rdkit.Chem import rdmolops, EnumerateStereoisomers, rdDistGeom, rdMolTransforms

import structure_step
import molsystem
import seamm
from seamm import standard_parameters
from seamm_ase import ASE_mixin
from seamm_geometric import geomeTRIC_mixin
from seamm_util import units_class, getParser
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Structure")

# Add this module's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class Structure(seamm.Node, ASE_mixin, geomeTRIC_mixin):
    """
    The non-graphical part of a Structure step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : StructureParameters
        The control parameters for Structure.

    See Also
    --------
    TkStructure,
    Structure, StructureParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Structure",
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        """A step for optimizing structure in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Structure {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Structure", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title="Structure",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = structure_step.metadata
        self.parameters = structure_step.StructureParameters()
        self._step = 0
        self._file_handler = None
        self._working_configuration = None
        self._working_directory = None
        self._data = {}
        self._results = {}
        self._logfile = None

    @property
    def version(self):
        """The semantic version of this module."""
        return structure_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return structure_step.__git_revision__

    def analyze(self, indent="", **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        table = {
            "": [
                "Converged?",
                "Energy",
                "# steps",
                "Max Force",
                "RMS Force",
                "Max Step",
            ],
            "Value": [
                "Yes" if self._results["converged"] else "**NO**",
                f"{self._results['energy']:.2f}",
                self._results["nsteps"],
                f"{self._results['maximum_gradient']:.6f}",
                f"{self._results['rms_gradient']:.6f}",
                f"{self._results['maximum_step']:.6f}",
            ],
            "Units": [
                "",
                "kJ/mol",
                "",
                "kJ/mol/Å",
                "kJ/mol/Å",
                "Å",
            ],
        }

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
        )
        length = len(tmp.splitlines()[0])
        text = "\n"
        text += "Optimization results".center(length)
        text += "\n"
        text += tmp
        text += "\n"
        printer.important(__(text, indent=11 * " ", wrap=False, dedent=False))

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "structure-step"
        parser = getParser()

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        super().create_parser(name=parser_name)

        if not parser_exists:
            # Any options for diffusivity itself
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        # Now need to walk through the steps in the subflowchart...
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def description_text(self, P=None, short=False, natoms=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P : dict
            An optional dictionary of the current values of the control
            parameters.

        short : bool
            If True, return a short description of the step.

        natoms : int
            The number of atoms in the structure.

        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        result = self.header + "\n"
        text = ""
        target = P["target"]
        if target == "stereoisomers":
            text += "Will generate up to {max stereoisomers} stereoisomers of the "
            text += "structure."
        elif target in ("minimum", "transition state"):
            approach = P["approach"]
            if approach == "Optimization":
                if P["optimizer"].lower().endswith("/geometric"):
                    result += self.describe_geomeTRIC_optimizer(P=P)
                else:
                    text += "The structure will be optimized using the "
                    text += "{optimizer} optimizer, converging to {convergence} "

                max_steps = P["max steps"]
                if (
                    natoms is not None
                    and isinstance(max_steps, str)
                    and "natoms" in max_steps
                ):
                    tmp = max_steps.split()
                    if "natoms" in tmp[0]:
                        max_steps = int(tmp[1]) * natoms
                    else:
                        max_steps = int(tmp[0]) * natoms
                text += f"with a maximum of {max_steps} steps."

                stop = P["continue if not converged"]
                if isinstance(stop, bool) and not stop or stop == "no":
                    text += " The workflow will continue if the structure "
                    text += "does not converge."
            else:
                raise RuntimeError(
                    f"Do not recognize approach '{approach}' for target '{target}'"
                )
        else:
            raise RuntimeError(f"Do not recognize target '{target}'")

        result += str(__(text, **P, indent=4 * " "))

        # Make sure the subflowchart has the data from the parent flowchart
        self.subflowchart.root_directory = self.flowchart.root_directory
        self.subflowchart.executor = self.flowchart.executor
        self.subflowchart.in_jobserver = self.subflowchart.in_jobserver

        if not short and target in ("minimum", "transition state"):
            # Get the first real node
            node = self.subflowchart.get_node("1").next()
            result += "\n\n    The energy and forces will be calculated as follows:\n"

            # Now walk through the steps in the subflowchart...
            while node is not None:
                try:
                    result += str(
                        __(node.description_text(), indent=7 * " ", wrap=False)
                    )
                except Exception as e:
                    print(f"Error describing structure flowchart: {e} in {node}")
                    self.logger.critical(
                        f"Error describing structure flowchart: {e} in {node}"
                    )
                    raise
                except:  # noqa: E722
                    print(
                        "Unexpected error describing structure flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    self.logger.critical(
                        "Unexpected error describing structure flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    raise
                result += "\n"
                node = node.next()

        return result

    def generate_stereoisomers(self, P):
        """Generate the stereoisomers.

        Parameters
        ----------
        P : dict(str, value)
            The control parameters
        """
        _, starting_configuration = self.get_system_configuration()
        molecule = starting_configuration.to_RDKMol()

        options = EnumerateStereoisomers.StereoEnumerationOptions(
            unique=True, tryEmbedding=True, maxIsomers=P["max stereoisomers"]
        )

        n_max = EnumerateStereoisomers.GetStereoisomerCount(molecule, options=options)
        printer.important(
            __(
                f"The upper bound on the number of stereoisomers is {n_max}.",
                indent=4 * " ",
            )
        )

        isomers = tuple(
            EnumerateStereoisomers.EnumerateStereoisomers(molecule, options=options)
        )

        # If only one isomer, may need to embed it
        if len(isomers) == 1 and isomers[0].GetNumConformers() == 0:
            for ring in rdmolops.GetSSSR(molecule):
                if len(ring) <= 4:
                    ps = rdDistGeom.srETKDGv3()
                    break
            else:
                ps = rdDistGeom.ETKDGv3()

            ps.trackFailures = True
            conformer = rdDistGeom.EmbedMolecule(molecule, ps)
            if conformer == -1:
                raise RuntimeError("Could not embed {SMILES}")

            isomers = (molecule,)

        first = True
        for isomer in isomers:
            rdMolTransforms.CanonicalizeMol(isomer)
            conformer = [x for x in isomer.GetConformers()][-1]

            system, configuration = self.get_system_configuration(
                P, first=first, same_as="current"
            )

            configuration.coordinates_from_RDKMol(conformer)
            standard_parameters.set_names(system, configuration, P, _first=first)
            first = False

        printer.important(
            __(
                f"Generated {len(isomers)} stereoisomers.",
                indent=4 * " ",
            )
        )
        printer.important("")

        # Reference
        citations = molsystem.rdkit_citations()
        for i, citation in enumerate(citations, start=1):
            self.references.cite(
                raw=citation,
                alias=f"rdkit_{i}",
                module="structure_step",
                level=1,
                note=f"The principle citation #{i} for RDKit.",
            )

    def optimize(self, P, PP):
        """Optimize the structure to a minimum or transition state

        Parameters
        ----------
        P : dict(str, value)
            The control parameters

        PP : dict(str, value)
            The printable control parameters
        """
        self._data = {
            "step": [],
            "energy": [],
            "max_force": [],
            "rms_force": [],
            "max_step": [],
        }
        self._last_coordinates = None
        self._step = 0

        # Get the final configuration
        _, self._working_configuration = self.get_system_configuration(P)

        tic = time.perf_counter_ns()
        if P["approach"].lower() == "optimization":
            if P["optimizer"].lower().endswith("/ase"):
                self.run_ase_optimizer(P, PP)
            elif P["optimizer"].lower().endswith("/geometric"):
                try:
                    self.run_geomeTRIC_optimizer(P, PP)
                except Exception as e:
                    if "did not converge" in str(e):
                        if not P["continue if not converged"]:
                            raise
                    else:
                        raise
            else:
                raise ValueError(f"Unknown optimizer '{P['optimizer']}' in Structure")
        else:
            raise ValueError(f"Unknown approach '{P['approach']}' in Structure")
        toc = time.perf_counter_ns()
        self._results["t_elapsed"] = round((toc - tic) * 1.0e-9, 3)

        # Print the results
        self.analyze()

        # Store results to db, variables, tables, and json as requested
        self.store_results(
            configuration=self._working_configuration,
            data=self._results,
        )

    def plot(self, E_units="", F_units=""):
        """Generate a plot of the convergence of the geometry optimization."""
        figure = self.create_figure(
            module_path=("seamm",),
            template="line.graph_template",
            title="Geometry optimization convergence",
        )
        plot = figure.add_plot("convergence")

        x_axis = plot.add_axis("x", label="Step", start=0, stop=0.8)
        y_axis = plot.add_axis("y", label=f"Energy ({E_units})")
        y2_axis = plot.add_axis(
            "y",
            anchor=x_axis,
            label=f"Force ({F_units})",
            overlaying="y",
            side="right",
            tickmode="sync",
        )
        y3_axis = plot.add_axis(
            "y",
            anchor=None,
            label="Distance (Å)",
            overlaying="y",
            position=0.9,
            side="right",
            tickmode="sync",
        )
        x_axis.anchor = y_axis

        plot.add_trace(
            color="red",
            name="Energy",
            width=3,
            x=self._data["step"],
            x_axis=x_axis,
            xlabel="step",
            y=self._data["energy"],
            y_axis=y_axis,
            ylabel="Energy",
            yunits=E_units,
        )

        plot.add_trace(
            color="black",
            name="Max Force",
            width=3,
            x=self._data["step"],
            x_axis=x_axis,
            xlabel="step",
            y=self._data["max_force"],
            y_axis=y2_axis,
            ylabel="Max Force",
            yunits=F_units,
        )

        plot.add_trace(
            color="green",
            name="RMS Force",
            width=3,
            x=self._data["step"],
            x_axis=x_axis,
            xlabel="step",
            y=self._data["rms_force"],
            y_axis=y2_axis,
            ylabel="RMS Force",
            yunits=F_units,
        )

        plot.add_trace(
            color="blue",
            name="Max Step",
            width=3,
            x=self._data["step"],
            x_axis=x_axis,
            xlabel="step",
            y=self._data["max_step"],
            y_axis=y3_axis,
            ylabel="Max Step",
            yunits="Å",
        )

        figure.grid_plots("convergence")

        # Write to disk
        path = Path(self.directory) / "Convergence.graph"
        figure.dump(path)

        if "html" in self.options and self.options["html"]:
            path = Path(self.directory) / "Convergence.html"
            figure.template = "line.html_template"
            figure.dump(path)

    def run(self):
        """Run a Structure step.

        Parameters
        ----------
        None

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        _, configuration = self.get_system_configuration()
        n_atoms = configuration.n_atoms

        # Print what we are doing
        printer.important(
            __(
                self.description_text(PP, short=True, natoms=n_atoms),
                indent=self.indent,
            )
        )

        # Just do it!
        target = P["target"]
        if target == "stereoisomers":
            self.generate_stereoisomers(P)
        elif target in ("minimum", "transition state"):
            approach = P["approach"]
            if approach == "Optimization":
                self.optimize(P, PP)
            else:
                raise RuntimeError(
                    f"Do not recognize approach '{approach}' for target '{target}'"
                )
        else:
            raise RuntimeError(f"Do not recognize target '{target}'")

        return next_node

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        self.logger.debug("Setting ids for subflowchart {}".format(self))
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1
