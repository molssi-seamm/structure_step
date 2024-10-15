# -*- coding: utf-8 -*-

"""Non-graphical part of the Structure step in a SEAMM flowchart
"""

from datetime import datetime
import json
import logging
from pathlib import Path
import pkg_resources
import pprint  # noqa: F401
import shutil
import sys
import time
import traceback

from ase import Atoms as ASE_Atoms
import ase.optimize as ASE_Optimize
from ase.calculators.calculator import all_changes as ASE_all_changes
import numpy as np
from tabulate import tabulate

import structure_step
import molsystem
import seamm
from seamm_ase import SEAMM_Calculator
from seamm_geometric import geomeTRIC_mixin
from seamm_util import Q_, units_class, getParser
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


class Structure(seamm.Node, geomeTRIC_mixin):
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

    def calculator(
        self,
        calculator,
        properties=["energy"],
        system_changes=ASE_all_changes,
    ):
        """Create a calculator for the structure step.

        Parameters
        ----------
        ase : ase.calculators.calculator.Calculator
            The ASE calculator we are working for
        properties : list of str
            The properties to calculate.
        system_changes : int
            The changes to the system.

        Returns
        -------
        results : dict
            The dictionary of results from the calculation.
        """
        wd = Path(self.directory)
        wd.mkdir(parents=True, exist_ok=True)

        self._step += 1
        self._results["nsteps"] = self._step
        self._data["step"].append(self._step)
        fmt = "05d"

        calculator.results = {}

        n_atoms = len(calculator.atoms)
        self.logger.debug(f"{n_atoms} atoms in the structure")
        positions = calculator.atoms.positions
        self.logger.debug(f"Positions: {positions}")
        cell = calculator.atoms.cell
        self.logger.debug(f"Cell: {cell}")

        # Set the coordinates in the configuration
        self._working_configuration.atoms.set_coordinates(positions, fractionals=False)

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        # Get the first real node
        first_node = self.subflowchart.get_node("1").next()

        # Ensure the nodes have their options
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()

        # And the subflowchart has the executor
        self.subflowchart.executor = self.flowchart.executor

        # Direct most output to iteration.out
        step_id = f"step_{self._step:{fmt}}"
        step_dir = Path(self.directory) / step_id
        step_dir.mkdir(parents=True, exist_ok=True)

        # A handler for the file
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
        path = step_dir / "Step.out"
        path.unlink(missing_ok=True)
        self._file_handler = logging.FileHandler(path)
        self._file_handler.setLevel(printing.NORMAL)
        formatter = logging.Formatter(fmt="{message:s}", style="{")
        self._file_handler.setFormatter(formatter)
        job.addHandler(self._file_handler)

        # Add the step to the ids so the directory structure is reasonable
        self.subflowchart.reset_visited()
        self.set_subids((*self._id, step_id))

        # Run through the steps in the loop body
        node = first_node
        try:
            while node is not None:
                node = node.run()
        except DeprecationWarning as e:
            printer.normal("\nDeprecation warning: " + str(e))
            traceback.print_exc(file=sys.stderr)
            traceback.print_exc(file=sys.stdout)
        except Exception as e:
            printer.job(f"Caught exception in step {self.step}: {str(e)}")
            with open(step_dir / "stderr.out", "a") as fd:
                traceback.print_exc(file=fd)
            raise
        self.logger.debug(f"End of step {self._step}")

        # Remove any redirection of printing.
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        # Get the energy and derivatives
        paths = sorted(step_dir.glob("**/Results.json"))

        if len(paths) == 0:
            raise RuntimeError(
                "There are no energy and gradients in properties.json for step "
                f"{self._step} in {step_dir}."
            )
        else:
            # Find the most recent and assume that is the one wanted
            newest_time = None
            for path in paths:
                with path.open() as fd:
                    data = json.load(fd)
                time = datetime.fromisoformat(data["iso time"])
                if newest_time is None:
                    newest = path
                    newest_time = time
                elif time > newest_time:
                    newest_time = time
                    newest = path
            with newest.open() as fd:
                data = json.load(fd)

        energy = data["energy"]
        if "energy,units" in data:
            units = data["energy,units"]
        else:
            units = "kJ/mol"
        self._data["energy"].append(energy)

        energy *= Q_(1.0, units).to("eV").magnitude
        self._results["energy"] = Q_(energy, "eV").m_as("kJ/mol")

        gradients = data["gradients"]

        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug("\ngradients")
            for i in range(n_atoms):
                self.logger.debug(
                    f"   {gradients[i][0]:8.3f} {gradients[i][1]:8.3f} "
                    f"{gradients[i][2]:8.3f}"
                )

        if "gradients,units" in data:
            funits = data["gradients,units"]
        else:
            funits = "kJ/mol/Å"

        # Get the measures of convergence
        max_force = np.max(np.linalg.norm(gradients, axis=1))
        self._data["max_force"].append(max_force)
        self._results["maximum_gradient"] = Q_(max_force, funits).m_as("kJ/mol/Å")
        rms_force = np.sqrt(np.mean(np.linalg.norm(gradients, axis=1) ** 2))
        self._data["rms_force"].append(rms_force)
        self._results["rms_gradient"] = Q_(rms_force, funits).m_as("kJ/mol/Å")

        if self._step > 1:
            step = positions - self._last_coordinates
            max_step = np.max(np.linalg.norm(step, axis=1))
        else:
            max_step = 0.0
        self._data["max_step"].append(max_step)
        self._results["maximum_step"] = max_step
        self._last_coordinates = np.array(positions)

        # Units!
        gradients = np.array(gradients) * Q_(1.0, funits).to("eV/Å").magnitude

        calculator.results["energy"] = energy
        calculator.results["forces"] = -gradients

        # Log the results
        if self._logfile is not None:
            headers = [
                "Step",
                f"E ({units})",
                f"Fmax ({funits})",
                f"Frms ({funits})",
                "max step (Å)",
            ]
            tmp = tabulate(
                self._data,
                headers=headers,
                tablefmt="rounded_outline",
                disable_numparse=False,
                floatfmt=".3f",
            )
            with open(self._logfile, "w") as fd:
                fd.write(tmp)
                fd.write("\n")

        # and plot the results
        self.plot(E_units=units, F_units=funits)

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
        if P["approach"] == "Optimization":
            text = ""
            if P["optimizer"].lower().endswith("/geometric"):
                result += self.describe_geomeTRIC_optimizer(P=P)
            else:
                text += "The structure will be optimized using the "
                text += "{optimizer} optimizer, converging to {convergence} "

            max_steps = P["max steps"]
            if natoms is not None and "natoms" in max_steps:
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

        result += "\n" + str(__(text, **P, indent=7 * " "))

        # Make sure the subflowchart has the data from the parent flowchart
        self.subflowchart.root_directory = self.flowchart.root_directory
        self.subflowchart.executor = self.flowchart.executor
        self.subflowchart.in_jobserver = self.subflowchart.in_jobserver

        if not short:
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
        self._data = {
            "step": [],
            "energy": [],
            "max_force": [],
            "rms_force": [],
            "max_step": [],
        }
        self._last_coordinates = None
        self._step = 0
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

        # Get the final configuration
        _, self._working_configuration = self.get_system_configuration(P)
        n_atoms = self._working_configuration.n_atoms

        # Print what we are doing
        printer.important(
            __(
                self.description_text(PP, short=True, natoms=n_atoms),
                indent=self.indent,
            )
        )

        if P["approach"].lower() == "optimization":
            if P["optimizer"].lower().endswith("/ase"):
                self.run_ase_optimizer(P, PP)
            elif P["optimizer"].lower().endswith("/geometric"):
                self.run_geomeTRIC_optimizer(P, PP)
            else:
                raise ValueError(f"Unknown optimizer '{P['optimizer']}' in Structure")
        else:
            raise ValueError(f"Unknown approach '{P['approach']}' in Structure")

        # Print the results
        self.analyze()

        # Store results to db, variables, tables, and json as requested
        self.store_results(
            configuration=self._working_configuration,
            data=self._results,
        )

        return next_node

    def run_ase_optimizer(self, P, PP):
        """Run a Structure step.

        Parameters
        ----------
        P : dict
            The current values of the parameters
        PP : dict
            The current values of the parameters, formatted for printing
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

        _, starting_configuration = self.get_system_configuration()
        n_atoms = starting_configuration.n_atoms

        # Print what we are doing
        printer.important(
            __(self.description_text(P, short=True, natoms=n_atoms), indent=self.indent)
        )

        # Create the directory
        wd = Path(self.directory)
        wd.mkdir(parents=True, exist_ok=True)

        # Setup the log file for the optimization
        self._logfile = wd / "optimization.log"

        symbols = starting_configuration.atoms.symbols
        XYZ = starting_configuration.atoms.coordinates

        calculator = SEAMM_Calculator(self)
        atoms = ASE_Atoms("".join(symbols), positions=XYZ, calculator=calculator)

        # The default maximum number of steps may vary depending on the optimizer
        max_steps = P["max steps"]

        # Optimize the structure
        optimizer = P["optimizer"][4:].lower()
        if optimizer == "bfgs":
            optimizer = ASE_Optimize.BFGS(atoms, restart=wd / "bfgs.json", logfile=None)
        elif optimizer == "lbfgs":
            optimizer = ASE_Optimize.LBFGS(
                atoms, restart=wd / "lbfgs.json", logfile=None
            )
        elif optimizer == "fire":
            optimizer = ASE_Optimize.FIRE(atoms, restart=wd / "fire.json", logfile=None)
        elif optimizer == "gpmin":
            optimizer = ASE_Optimize.GPMin(
                atoms, restart=wd / "gpmin.json", logfile=None
            )
        elif optimizer == "mdmin":
            optimizer = ASE_Optimize.MDMin(
                atoms, restart=wd / "mdmin.json", logfile=None
            )
        elif optimizer == "bfgslinesearch":
            optimizer = ASE_Optimize.BFGSLineSearch(
                atoms, restart=wd / "bfgsline.json", logfile=None
            )
        elif optimizer == "lbfgslinesearch":
            optimizer = ASE_Optimize.LBFGSLineSearch(
                atoms, restart=wd / "lbfgsline.json", logfile=None
            )
        else:
            raise ValueError(
                f"Unknown optimizer '{optimizer}' ({P['optimizer']}) in Structure"
            )

        convergence = P["convergence"].m_as("eV/Å")
        if "natoms" in max_steps:
            tmp = max_steps.split()
            if "natoms" in tmp[0]:
                max_steps = int(tmp[1]) * len(atoms)
            else:
                max_steps = int(tmp[0]) * len(atoms)

        # Run the optimization
        exception = None
        tic = time.perf_counter_ns()
        try:
            converged = optimizer.run(fmax=convergence, steps=max_steps)
            print(f"Converged = {converged}")
        except Exception as exception:  # noqa: F841
            print(f"Exception: {exception}")
            converged = False
        finally:
            toc = time.perf_counter_ns()
            self._results["t_elapsed"] = round((toc - tic) * 1.0e-9, 3)
            self._results["converged"] = converged

            # Print the results
            self.analyze()

            # Store results to db, variables, tables, and json as requested
            self.store_results(
                configuration=self._working_configuration,
                data=self._results,
            )

        # Clean up the subdirectories
        if exception is not None or not converged:
            keep = P["on error"]
            if keep == "delete all subdirectories":
                subdirectories = wd.glob("step_*")
                for subdirectory in subdirectories:
                    shutil.rmtree(subdirectory)
            elif keep == "keep last subdirectory":
                subdirectories = wd.glob("step_*")
                subdirectories = sorted(subdirectories)
                for subdirectory in subdirectories[:-1]:
                    shutil.rmtree(subdirectory)
            if not converged:
                raise RuntimeError(
                    f"Optimization did not converge in {max_steps} steps"
                )
            raise exception from None
        else:
            keep = P["on success"]
            if keep == "delete all subdirectories":
                subdirectories = wd.glob("step_*")
                for subdirectory in subdirectories:
                    shutil.rmtree(subdirectory)
            elif keep == "keep last subdirectory":
                subdirectories = wd.glob("step_*")
                subdirectories = sorted(subdirectories)
                for subdirectory in subdirectories[:-1]:
                    shutil.rmtree(subdirectory)

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
