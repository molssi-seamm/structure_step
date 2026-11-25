# -*- coding: utf-8 -*-

"""Where the Structure step's energies and gradients come from.

Either the step's own sub-steps (a subflowchart, e.g. MOPAC with an Energy
sub-step, run at each geometry and read back from its Results.json) or the
flowchart's Model Chemistry (the ``_model_chemistry`` variable published by a
Model Chemistry step upstream). The Model Chemistry is driven as a warm MDI
engine when the program has one here, otherwise as one separate calculation
(task) per geometry through :class:`seamm_exec.Evaluator`. Every optimizer --
geomeTRIC, ASE and SEAMM -- sees the same engine interface
(``set_coordinates``, ``energy``, ``forces``) in an
:class:`seamm_exec.hessian.EngineSession`.
"""

import contextlib
import hashlib

import numpy as np

from seamm_util import Q_

#: The message when a Structure step has no way to get energies and gradients
NO_SOURCE = (
    "The Structure step has no way to calculate the energy and gradients: add a "
    "Model Chemistry step before it, or put the calculation (e.g. MOPAC or ORCA "
    "with an Energy sub-step) in its subflowchart."
)

#: ... and for the SEAMM optimizer, which uses only the Model Chemistry
SEAMM_NEEDS_MODEL_CHEMISTRY = (
    "The SEAMM optimizer gets the energy and gradients from the flowchart's Model "
    "Chemistry, and there is none: add a Model Chemistry step before the Structure "
    "step, or choose another optimizer (e.g. geomeTRIC) to use its sub-steps."
)


def is_model_chemistry_step(node):
    """Whether a node is a Model Chemistry step (without importing that plug-in)."""
    return getattr(node, "extension", None) == "Model Chemistry" or (
        type(node).__name__ == "ModelChemistry"
        and type(node).__module__.startswith("model_chemistry_step")
    )


def has_substeps(flowchart):
    """Whether a subflowchart has any step after its start node."""
    try:
        return flowchart.get_node("1").next() is not None
    except Exception:
        return False


def uses_model_chemistry_only(optimizer):
    """Whether the optimizer works only from the Model Chemistry (SEAMM)."""
    return isinstance(optimizer, str) and optimizer.lower().endswith("/seamm")


class EvaluatorEngine:
    """The interface of an MDI engine over separate calculations: each new geometry
    is one task of a :class:`seamm_exec.Evaluator`, for a Model Chemistry with no
    MDI engine here (or whose program runs only where the job's tasks run).

    The energy and forces of a geometry come from the same task, computed once.
    A task's key carries a hash of its coordinates, so a job that is run again
    reuses only the calculations of the same geometries.

    Parameters
    ----------
    driver : seamm_exec.hessian.HessianDriver
        Its ``evaluate`` runs the tasks.
    evaluator : seamm_exec.Evaluator
    configuration : molsystem._Configuration
        For the elements, charge and multiplicity.
    """

    def __init__(self, driver, evaluator, configuration, prefix="opt"):
        self.driver = driver
        self.evaluator = evaluator
        self.configuration = configuration
        self.prefix = prefix
        self.n_calculations = 0
        self._xyz = None
        self._result = None

    def supports(self, command):
        """No MDI commands, so no analytic Hessian."""
        return False

    def set_coordinates(self, xyz, units="Å"):
        x = np.asarray(xyz, dtype=float).reshape(-1, 3)
        if units not in (None, "Å", "angstrom", "Angstrom"):
            x = Q_(x, units).m_as("Å")
        if self._xyz is None or not np.array_equal(x, self._xyz):
            self._xyz = np.array(x)
            self._result = None

    def _calculate(self):
        if self._xyz is None:
            raise RuntimeError("No coordinates were given to the engine.")
        if self._result is None:
            self.n_calculations += 1
            digest = hashlib.sha1(np.round(self._xyz, 8).tobytes()).hexdigest()[:10]
            key = f"{self.prefix}{self.n_calculations:05d}-{digest}"
            results = self.driver.evaluate(
                self.evaluator,
                self.configuration,
                {key: self._xyz},
                what="the optimization",
            )
            self._result = results[key]  # kJ/mol, E_h/bohr (3N,)
        return self._result

    def energy(self, units="kJ/mol"):
        return Q_(self._calculate()[0], "kJ/mol").m_as(units)

    def forces(self, units="hartree/bohr"):
        gradient = np.asarray(self._calculate()[1], dtype=float).reshape(-1, 3)
        return -Q_(gradient, "hartree/bohr").m_as(units)


def _quiet_driver_class():
    """A HessianDriver whose route choice is silent: its messages are about
    Hessians; the Structure step says what it does itself."""
    from seamm_exec.hessian import HessianDriver

    class _OptimizationDriver(HessianDriver):
        _quiet = False

        def say(self, text):
            if not self._quiet:
                super().say(text)

        def plan(self, configuration):
            self._quiet = True
            try:
                return super().plan(configuration)
            finally:
                self._quiet = False

    return _OptimizationDriver


def make_driver(node, print_to=None):
    """The driver of the Model Chemistry for an optimization: a warm MDI engine
    whenever the program has one here, otherwise separate calculations.

    Returns
    -------
    seamm_exec.hessian.HessianDriver
    """
    mc = node.get_variable("_model_chemistry")
    options = mc.get("options") or {}
    try:
        step = node.flowchart.plugin_manager.get(mc["step"])
    except Exception:
        step = None
    has_mdi = bool(options.get("mdi_capable", False)) and hasattr(
        step, "get_mdi_engine_command"
    )
    # An optimization is one geometry after another, so a warm engine in the job
    # beats a task per step even when the job's tasks go to a queue; the driver
    # still uses tasks when the program is not installed here.
    route = "warm engine" if has_mdi else "separate tasks"
    return _quiet_driver_class()(
        node, mc, route=route, who="The Structure step", print_to=print_to
    )


@contextlib.contextmanager
def model_chemistry_session(node, configuration, driver=None, print_to=None):
    """An :class:`seamm_exec.hessian.EngineSession` over the flowchart's Model
    Chemistry for the configuration: a warm MDI engine, or separate calculations.

    Parameters
    ----------
    node : seamm.Node
        The Structure step.
    configuration : molsystem._Configuration
        The structure to optimize.
    driver : seamm_exec.hessian.HessianDriver, optional
        Default: :func:`make_driver`.
    print_to : printer, optional
        Where the notes for the user go.

    Yields
    ------
    seamm_exec.hessian.EngineSession
    """
    from seamm_exec.hessian import EngineSession

    if driver is None:
        driver = make_driver(node, print_to=print_to)
    mc = driver.model_chemistry()
    level = mc.get("level", "the Model Chemistry")
    if configuration.periodicity != 0:
        raise RuntimeError(
            f"The Structure step cannot optimize a periodic system with the Model "
            f"Chemistry ({level}); put the calculation in its subflowchart instead."
        )
    with contextlib.ExitStack() as stack:
        try:
            session = stack.enter_context(driver.session(configuration))
        except ValueError as e:
            raise RuntimeError(
                f"The Structure step cannot get the energy and gradients from the "
                f"Model Chemistry {level}: it has neither an MDI engine nor separate "
                f"calculations here ({e})"
            ) from e
        if session is None:
            evaluator, _, elements, coords_ang = driver.last_plan
            driver.say(
                f"The energy and gradients come from {level}, as one separate "
                "calculation per optimization step."
            )
            engine = EvaluatorEngine(driver, evaluator, configuration)
            session = EngineSession(
                driver, engine, configuration, evaluator, elements, coords_ang
            )
        else:
            driver.say(f"The energy and gradients come from the {level} MDI engine.")
        yield session
