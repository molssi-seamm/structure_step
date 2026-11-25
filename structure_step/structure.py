# -*- coding: utf-8 -*-

"""Non-graphical part of the Structure step in a SEAMM flowchart"""

import contextlib
import importlib
import logging
from pathlib import Path
import pprint  # noqa: F401
import sys
import time

import numpy as np
from tabulate import tabulate
from rdkit.Chem import rdmolops, EnumerateStereoisomers, rdDistGeom, rdMolTransforms

import structure_step
import molsystem
import seamm
from seamm import standard_parameters
from seamm_ase import ASE_mixin
from ase.calculators.calculator import all_changes as ASE_all_changes
from seamm_geometric import geomeTRIC_mixin
from seamm_util import Q_, units_class, getParser
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
from .energy_source import (
    has_substeps,
    is_model_chemistry_step,
    make_driver,
    model_chemistry_session,
    NO_SOURCE,
    SEAMM_NEEDS_MODEL_CHEMISTRY,
    uses_model_chemistry_only,
)

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
path = importlib.resources.files("structure_step") / "data"
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
        # The Model Chemistry's engine session while an optimizer uses it (else
        # the energies and gradients come from the sub-steps)
        self._mc_session = None

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

        if "is_minimum" in self._results:
            table[""].append("Minimum?")
            table["Value"].append("Yes" if self._results["is_minimum"] else "**NO**")
            table["Units"].append("")
            if "lowest_frequency" in self._results:
                table[""].append("Lowest frequency")
                table["Value"].append(f"{self._results['lowest_frequency']:.1f}")
                table["Units"].append("cm^-1")
            for key, label in (
                ("n_escapes", "Saddle escapes"),
                ("n_confirmations", "Modes confirmed"),
            ):
                table[""].append(label)
                table["Value"].append(self._results[key])
                table["Units"].append("")
        if "n_cheap_hessians" in self._results:
            table[""].append("Cheap Hessians")
            table["Value"].append(self._results["n_cheap_hessians"])
            table["Units"].append("")

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

    def apply_convergence_preset(self, P):
        """Set the convergence formula and criteria from the preset convergence.

        The dialog does this when a preset (anything but "custom") is chosen, but a
        flowchart made without it -- by seamm-flowchart or the MCP server -- has only
        the preset's name: the formula was left empty, failing with KeyError '', and
        the criteria were the defaults rather than the preset's.

        Parameters
        ----------
        P : dict
            The values of the parameters, changed in place.

        Returns
        -------
        dict
            P
        """
        convergence = P.get("convergence")
        presets = self.metadata["convergence parameters"]
        if not isinstance(convergence, str) or convergence not in presets:
            return P
        for key, value in presets[convergence].items():
            if key == "convergence formula":
                P[key] = value
            elif "criterion" in key:
                P[key] = Q_(float(value[0]), value[1])
        return P

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
            P = self.apply_convergence_preset(self.parameters.values_to_dict())

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
                elif P["optimizer"].lower().endswith("/seamm"):
                    text += self.describe_seamm_optimizer(P)
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

        if (
            not short
            and target in ("minimum", "transition state")
            and (
                uses_model_chemistry_only(P["optimizer"])
                or not has_substeps(self.subflowchart)
            )
        ):
            result += (
                "\n\n    The energy and forces will come from the flowchart's Model "
                "Chemistry.\n"
            )
        elif not short and target in ("minimum", "transition state"):
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
        # results from an earlier run of this step (e.g. in a loop) must not carry over
        self._results = {}

        # Get the final configuration
        _, self._working_configuration = self.get_system_configuration(P)

        tic = time.perf_counter_ns()
        if P["approach"].lower() == "optimization":
            optimizer = P["optimizer"].lower()
            if optimizer.endswith("/seamm"):
                self.energy_source(P["optimizer"])
                self.run_seamm_optimizer(P, PP)
            elif optimizer.endswith("/ase") or optimizer.endswith("/geometric"):
                with self.energy_engine(P["optimizer"]):
                    if optimizer.endswith("/ase"):
                        self.run_ase_optimizer(P, PP)
                    else:
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

    # ------------------------------------------------------------------
    # Where the energies and gradients come from
    # ------------------------------------------------------------------
    def energy_source(self, optimizer):
        """Where the optimizer's energies and gradients come from.

        The sub-steps when there are any (a note says so if a Model Chemistry is
        upstream too, since it is not used here); else the flowchart's Model
        Chemistry. The SEAMM optimizer uses only the Model Chemistry.

        Parameters
        ----------
        optimizer : str
            The optimizer, e.g. "geomeTRIC/geomeTRIC".

        Returns
        -------
        str
            "sub-steps" or "model chemistry"

        Raises
        ------
        RuntimeError
            When there is neither, naming what is missing.
        """
        have_mc = self.variable_exists("_model_chemistry")
        if uses_model_chemistry_only(optimizer):
            if not have_mc:
                raise RuntimeError(SEAMM_NEEDS_MODEL_CHEMISTRY)
            return "model chemistry"
        if has_substeps(self.subflowchart):
            if have_mc:
                level = self.get_variable("_model_chemistry").get("level", "")
                printer.important(
                    __(
                        "The energy and gradients come from this step's sub-steps; "
                        f"the flowchart's Model Chemistry ({level}) is not used by "
                        "this step.",
                        indent=self.indent + 4 * " ",
                    )
                )
            return "sub-steps"
        if have_mc:
            return "model chemistry"
        raise RuntimeError(NO_SOURCE)

    @contextlib.contextmanager
    def energy_engine(self, optimizer):
        """For the geomeTRIC and ASE optimizers: while open, their energies and
        gradients come from the Model Chemistry's engine (``self._mc_session``)
        when that is the source, else from the sub-steps."""
        if self.energy_source(optimizer) == "sub-steps":
            yield None
            return
        with model_chemistry_session(
            self, self._working_configuration, print_to=printer
        ) as session:
            self._mc_session = session
            try:
                yield session
            finally:
                self._mc_session = None

    def model_chemistry_energy_and_gradients(self, coordinates):
        """The energy (kJ/mol) and gradients (kJ/mol/Å, (n, 3)) at the coordinates
        (Å) from the Model Chemistry's engine; also sets the working
        configuration's coordinates and records the step for the convergence
        table and graph."""
        self._step += 1
        xyz = np.asarray(coordinates, dtype=float).reshape(-1, 3)
        self._working_configuration.atoms.set_coordinates(xyz, fractionals=False)
        engine = self._mc_session.engine
        engine.set_coordinates(xyz, units="Å")
        energy = Q_(engine.energy(units="eV"), "eV").m_as("kJ/mol")
        forces = np.asarray(engine.forces(units="eV/Å"), dtype=float).reshape(-1, 3)
        gradients = -Q_(forces, "eV/Å").m_as("kJ/mol/Å")

        norms = np.linalg.norm(gradients, axis=1)
        max_force = float(np.max(norms))
        rms_force = float(np.sqrt(np.mean(norms**2)))
        if self._step > 1 and self._last_coordinates is not None:
            max_step = float(
                np.max(np.linalg.norm(xyz - self._last_coordinates, axis=1))
            )
        else:
            max_step = 0.0
        self._last_coordinates = np.array(xyz)
        for key, value in (
            ("step", self._step),
            ("energy", energy),
            ("max_force", max_force),
            ("rms_force", rms_force),
            ("max_step", max_step),
        ):
            self._data.setdefault(key, []).append(value)
        self._results.update(
            energy=energy,
            nsteps=self._step,
            maximum_gradient=max_force,
            rms_gradient=rms_force,
            maximum_step=max_step,
        )
        return energy, gradients

    @contextlib.contextmanager
    def _substeps_gave_no_results(self):
        """Say what is missing when the sub-steps ran but left no energy and
        gradients for the optimizer (the connectors' message spoke of a
        properties.json): a sub-step saves them only when its results ask for
        'energy' and 'gradients' in JSON."""
        try:
            yield
        except RuntimeError as e:
            if "There are no energy and gradients in" not in str(e):
                raise
            raise RuntimeError(
                "The sub-steps of the Structure step gave no energy and gradients "
                f"(step {self._step}): in the results of the calculation that gives "
                "them (e.g. the MOPAC Energy sub-step), save 'energy' and "
                "'gradients' to JSON; or remove the sub-steps to use the "
                "flowchart's Model Chemistry."
            ) from e

    def geometric_calculate_gradients(self, coordinates):
        """The energy (E_h) and gradients (E_h/bohr) for geomeTRIC: from the Model
        Chemistry's engine when it is the source, else the sub-steps (the
        geomeTRIC connector's own method)."""
        if self._mc_session is None:
            with self._substeps_gave_no_results():
                return super().geometric_calculate_gradients(coordinates)
        energy, gradients = self.model_chemistry_energy_and_gradients(coordinates)

        self.geometric_read_bibliography()
        if "seamm_geometric" in self._bibliography:
            self.references.cite(
                raw=self._bibliography["seamm_geometric"],
                alias="seamm_geometric",
                module="seamm_geometric",
                level=2,
                note="Main reference for the geomeTRIC connector for SEAMM.",
            )
        return (
            Q_(energy, "kJ/mol").m_as("E_h"),
            Q_(gradients, "kJ/mol/Å").m_as("E_h/a_0"),
        )

    def ase_calculator(
        self, calculator, properties=["energy"], system_changes=ASE_all_changes
    ):
        """The energy and forces for an ASE optimizer: from the Model Chemistry's
        engine when it is the source, else the sub-steps (the ASE connector's own
        method)."""
        if self._mc_session is None:
            with self._substeps_gave_no_results():
                return super().ase_calculator(calculator, properties, system_changes)
        calculator.results = {}
        energy, gradients = self.model_chemistry_energy_and_gradients(
            calculator.atoms.positions
        )
        calculator.results["energy"] = Q_(energy, "kJ/mol").m_as("eV")
        calculator.results["forces"] = -Q_(gradients, "kJ/mol/Å").m_as("eV/Å")

        if self._logfile is not None:
            headers = [
                "Step",
                "E (kJ/mol)",
                "Fmax (kJ/mol/Å)",
                "Frms (kJ/mol/Å)",
                "max step (Å)",
            ]
            with open(self._logfile, "w") as fd:
                fd.write(
                    tabulate(
                        self._data,
                        headers=headers,
                        tablefmt="rounded_outline",
                        floatfmt=".3f",
                    )
                )
                fd.write("\n")
        self.plot(E_units="kJ/mol", F_units="kJ/mol/Å")

        self.ase_read_bibliography()
        if "seamm_ase" in self._bibliography:
            self.references.cite(
                raw=self._bibliography["seamm_ase"],
                alias="seamm_ase",
                module="seamm_ase",
                level=self.citation_level,
                note="The principle citation for the ASE connector for SEAMM.",
            )

    def check_flowchart(self, preceding):
        """Check, before the job runs, that an optimization can get its energy and
        gradients: from sub-steps, or from a Model Chemistry step before this one
        (the SEAMM optimizer only from the Model Chemistry).

        Called by SEAMM's flowchart checks (the editor's Run, ``seamm-flowchart
        validate`` and the MCP server), so such a flowchart is stopped before it
        is queued rather than failing at its first gradient.

        Parameters
        ----------
        preceding : [seamm.Node]
            The steps that run before this one.

        Returns
        -------
        [str]
            The problems found; empty if none.
        """
        P = self.parameters
        target = P["target"].value
        approach = P["approach"].value
        if not self.is_expr(target) and target not in ("minimum", "transition state"):
            return []
        if not self.is_expr(approach) and approach != "Optimization":
            return []
        if any(is_model_chemistry_step(node) for node in preceding):
            return []
        optimizer = P["optimizer"].value
        if self.is_expr(optimizer):
            return []
        if uses_model_chemistry_only(optimizer):
            return [f"'{self.title}': {SEAMM_NEEDS_MODEL_CHEMISTRY}"]
        if not has_substeps(self.subflowchart):
            return [f"'{self.title}': {NO_SOURCE}"]
        return []

    #: the SEAMM optimizer's names for the coordinate systems
    SEAMM_COORDINATES = {
        "TRIC": "TRIC",
        "TRIC-p": "primitive",
        "Cart": "cartesian",
        "DLC": "DLC",
        "HDLC": "HDLC",
    }

    def describe_seamm_optimizer(self, P):
        """The description of the SEAMM optimizer (seamm_optimizer)."""
        text = (
            "The structure will be optimized with the SEAMM optimizer in "
            "{coordinate system}, over the Model Chemistry's engine, until the "
            "largest force on an atom is below {Maximum atomic gradient criterion} "
        )
        cheap = str(P.get("cheap level", "none")).strip()
        if cheap.lower() in ("", "none"):
            text += "(Hessian: the {model hessian} model) "
        else:
            text += f"(Hessian: {cheap}, recomputed when the model fails"
            if P.get("check minimum") in (True, "yes"):
                text += "; the end point is checked to be a minimum"
            text += ") "
        return text

    def run_seamm_optimizer(self, P, PP):
        """Optimize with seamm_optimizer over the flowchart's Model Chemistry (a
        warm MDI engine), with a model or cheap-level Hessian."""
        configuration = self._working_configuration
        n_atoms = configuration.n_atoms
        fmax = P["Maximum atomic gradient criterion"]
        fmax = fmax.m_as("eV/Å") if hasattr(fmax, "m_as") else float(fmax)
        max_steps = P["max steps"]
        if isinstance(max_steps, str) and "natoms" in max_steps:
            tmp = max_steps.split()
            max_steps = int(tmp[1] if "natoms" in tmp[0] else tmp[0]) * n_atoms
        max_steps = int(max_steps)
        coordinates = self.SEAMM_COORDINATES[
            str(P["coordinate system"]).split(":")[0].strip()
        ]
        cheap = str(P.get("cheap level", "none")).strip()
        cheap = None if cheap.lower() in ("", "none") else cheap
        options = {
            "coordinates": coordinates,
            "check_minimum": P.get("check minimum") in (True, "yes")
            and cheap is not None,
        }
        if cheap is None:
            options["hessian"] = (
                "geometric" if P.get("model hessian") == "geomeTRIC guess" else "lindh"
            )

        driver = make_driver(self, print_to=printer)
        with model_chemistry_session(
            self, configuration, driver=driver, print_to=printer
        ) as session:
            converged = driver.optimize_in_blocks(
                session,
                fmax=fmax,
                max_steps=max_steps,
                optimizer="SEAMM",
                optimizer_options=options,
                cheap_level=cheap,
            )
        last = driver.last_optimization or {}
        self._results["converged"] = bool(converged)
        self._results["nsteps"] = last.get("evaluations", 0)
        if "energy (eV)" in last:
            self._results["energy"] = Q_(last["energy (eV)"], "eV").m_as("kJ/mol")
        for key, name in (
            ("max force (eV/Å)", "maximum_gradient"),
            ("rms force (eV/Å)", "rms_gradient"),
        ):
            self._results[name] = Q_(last.get(key, 0.0), "eV/Å").m_as("kJ/mol/Å")
        self._results["maximum_step"] = float(last.get("max step (Å)", 0.0))
        self._results.update(self.minimum_check_results(last))
        if not converged and P["continue if not converged"] not in (True, "yes"):
            raise RuntimeError(f"Optimization did not converge in {max_steps} steps")

    @staticmethod
    def minimum_check_results(last):
        """The results of the SEAMM optimizer's minimum check, from seamm_exec's
        record of the optimization (``last_optimization``): whether the end point is
        a minimum, its lowest frequency, the saddle escapes and the cheap-level modes
        confirmed, when the check ran (a cheap level and "check minimum"); and the
        number of cheap-level Hessians."""
        results = {}
        if last.get("minimum") is not None:
            results["is_minimum"] = bool(last["minimum"])
            if last.get("lowest frequency") is not None:
                results["lowest_frequency"] = float(last["lowest frequency"])
            results["n_escapes"] = int(last.get("escapes", 0))
            results["n_confirmations"] = int(last.get("confirmations", 0))
        if "cheap hessians" in last:
            results["n_cheap_hessians"] = int(last["cheap hessians"])
        return results

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
        self.apply_convergence_preset(P)

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
