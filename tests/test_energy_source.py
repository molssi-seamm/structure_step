#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Where the Structure step's energies and gradients come from: the sub-steps or
the flowchart's Model Chemistry, for every optimizer (2026-11-25: geomeTRIC, the
default, failed with "There are no energy and gradients in properties.json" for a
Model Chemistry and no sub-steps)."""

from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

from seamm_util import Q_
import structure_step
from structure_step import energy_source
from structure_step.energy_source import (
    EvaluatorEngine,
    model_chemistry_session,
    NO_SOURCE,
    SEAMM_NEEDS_MODEL_CHEMISTRY,
)

# Water, a little off its minimum, with harmonic bonds and an H-H spring
SYMBOLS = ["O", "H", "H"]
XYZ = np.array([[0.0, 0.0, 0.1], [0.0, 0.85, -0.5], [0.0, -0.75, -0.45]])
SPRINGS = {(0, 1): (0.96, 30.0), (0, 2): (0.96, 30.0), (1, 2): (1.52, 5.0)}


def harmonic(xyz):
    """Energy (eV) and gradients (eV/Å) of the springs."""
    xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
    energy = 0.0
    gradients = np.zeros_like(xyz)
    for (i, j), (r0, k) in SPRINGS.items():
        d = xyz[i] - xyz[j]
        r = np.linalg.norm(d)
        energy += 0.5 * k * (r - r0) ** 2
        g = k * (r - r0) * d / r
        gradients[i] += g
        gradients[j] -= g
    return energy, gradients


class FakeEngine:
    """A warm MDI engine with the harmonic potential."""

    def __init__(self):
        self.xyz = None
        self.n_calls = 0

    def supports(self, command):
        return False

    def set_coordinates(self, xyz, units="Å"):
        self.xyz = Q_(np.asarray(xyz, dtype=float), units).m_as("Å").reshape(-1, 3)

    def energy(self, units="eV"):
        self.n_calls += 1
        return Q_(harmonic(self.xyz)[0], "eV").m_as(units)

    def forces(self, units="eV/Å"):
        return -Q_(harmonic(self.xyz)[1], "eV/Å").m_as(units)


FAKE_MC = {"level": "Fake:SQM@harmonic", "step": "Fake", "method": "x", "options": {}}


@pytest.fixture()
def configuration():
    from molsystem import SystemDB

    db = SystemDB(filename="file:structure_test?mode=memory&cache=shared")
    system = db.create_system(name="water")
    configuration = system.create_configuration(name="start")
    configuration.atoms.append(
        x=XYZ[:, 0].tolist(),
        y=XYZ[:, 1].tolist(),
        z=XYZ[:, 2].tolist(),
        symbol=SYMBOLS,
    )
    yield configuration
    db.close()


class _Structure(structure_step.Structure):
    """A Structure step outside a flowchart: its own configuration and, if
    asked, a Model Chemistry."""

    def __init__(self, tmp_path, configuration, mc=True):
        super().__init__()
        self._test_configuration = configuration
        self._test_mc = mc
        self.flowchart = SimpleNamespace(
            root_directory=str(tmp_path), executor=None, in_jobserver=False
        )
        self._id = ("1",)
        self.references = MagicMock()

    def get_system_configuration(self, P=None, **kwargs):
        return self._test_configuration.system, self._test_configuration

    def variable_exists(self, name):
        return self._test_mc and name == "_model_chemistry"

    def get_variable(self, name):
        assert name == "_model_chemistry"
        return FAKE_MC

    def plot(self, E_units="", F_units=""):
        self.plotted = True


def _parameters(node, optimizer):
    node.parameters["optimizer"].value = optimizer
    P = node.parameters.current_values_to_dict()
    node.apply_convergence_preset(P)
    P["on success"] = "keep all subdirectories"
    return P, dict(P)


@pytest.fixture()
def fake_session(monkeypatch):
    """The Model Chemistry session gives the harmonic engine."""
    import contextlib

    engines = []

    @contextlib.contextmanager
    def session(node, configuration, driver=None, print_to=None):
        engine = FakeEngine()
        engines.append(engine)
        yield SimpleNamespace(engine=engine, configuration=configuration)

    monkeypatch.setattr(structure_step.structure, "model_chemistry_session", session)
    return engines


def _bond_lengths(configuration):
    xyz = np.array(configuration.atoms.get_coordinates(fractionals=False))
    return (
        np.linalg.norm(xyz[0] - xyz[1]),
        np.linalg.norm(xyz[0] - xyz[2]),
        np.linalg.norm(xyz[1] - xyz[2]),
    )


# ----------------------------------------------------------------------
# Choosing the source: 'Energies and gradients from:' (engine)
# ----------------------------------------------------------------------
MC, SUB = "Model Chemistry", "subflowchart"
GEO, SEAMM = "geomeTRIC/geomeTRIC", "SEAMM/SEAMM"


@pytest.mark.parametrize(
    "engine, optimizer, have_mc, have_sub, expected, note",
    [
        (MC, GEO, True, False, MC, None),
        (MC, GEO, True, True, MC, None),  # as chosen, though there are sub-steps
        (MC, GEO, False, True, SUB, energy_source.USING_SUBFLOWCHART),
        (SUB, GEO, True, True, SUB, None),  # as chosen, though there is an MC
        (SUB, GEO, False, True, SUB, None),
        (SUB, GEO, True, False, MC, energy_source.USING_MODEL_CHEMISTRY),
        (SUB, SEAMM, True, True, MC, None),  # SEAMM: only the Model Chemistry
        (MC, "BFGS/ASE", True, False, MC, None),
    ],
)
def test_choose_engine(engine, optimizer, have_mc, have_sub, expected, note):
    assert energy_source.choose_engine(engine, optimizer, have_mc, have_sub) == (
        expected,
        note,
    )


@pytest.mark.parametrize(
    "engine, optimizer, message",
    [
        (MC, GEO, NO_SOURCE),
        (SUB, GEO, NO_SOURCE),
        (SUB, SEAMM, SEAMM_NEEDS_MODEL_CHEMISTRY),
        (MC, SEAMM, SEAMM_NEEDS_MODEL_CHEMISTRY),
    ],
)
def test_choose_engine_neither_names_what_is_missing(engine, optimizer, message):
    with pytest.raises(RuntimeError) as e:
        energy_source.choose_engine(engine, optimizer, False, False)
    assert str(e.value) == message
    assert "properties.json" not in str(e.value)


def test_engine_parameter_mirrors_thermochemistry():
    """The same name, choices and default as the Thermochemistry step's."""
    parameter = structure_step.StructureParameters.parameters["engine"]
    assert parameter["default"] == "Model Chemistry"
    assert parameter["enumeration"] == ("Model Chemistry", "subflowchart")
    assert parameter["description"] == "Energies and gradients from:"


@pytest.mark.parametrize(
    "optimizer, expected", [(GEO, SUB), ("BFGS/ASE", SUB), (SEAMM, MC)]
)
def test_an_older_flowchart_keeps_its_source(optimizer, expected):
    """Saved before 'engine' existed: geomeTRIC and ASE ran the subflowchart and
    the SEAMM optimizer used the Model Chemistry, so they still do."""
    data = {"optimizer": {"value": optimizer, "units": None}}
    P = structure_step.StructureParameters(data=data)
    assert P["engine"].value == expected
    # a flowchart that has 'engine' keeps it; a new step has the default
    data["engine"] = {"value": MC, "units": None}
    assert structure_step.StructureParameters(data=data)["engine"].value == MC
    assert structure_step.StructureParameters()["engine"].value == MC


def test_source_follows_the_choice(tmp_path, configuration, monkeypatch):
    node = _Structure(tmp_path, configuration)
    P = {"engine": MC, "optimizer": GEO}
    assert node.energy_source(P) == MC
    monkeypatch.setattr(structure_step.structure, "has_substeps", lambda f: True)
    assert node.energy_source({"engine": SUB, "optimizer": GEO}) == SUB
    assert node.energy_source(P) == MC


def test_source_neither_names_what_is_missing(tmp_path, configuration):
    node = _Structure(tmp_path, configuration, mc=False)
    with pytest.raises(RuntimeError) as e:
        node.energy_source({"engine": MC, "optimizer": GEO})
    assert str(e.value) == NO_SOURCE


def test_optimize_with_neither_fails_before_any_calculation(tmp_path, configuration):
    node = _Structure(tmp_path, configuration, mc=False)
    P, PP = _parameters(node, "geomeTRIC/geomeTRIC")
    with pytest.raises(RuntimeError, match="Model Chemistry step before it"):
        node.optimize(P, PP)
    assert not (tmp_path / "1" / "geomeTRIC" / "step_00001").exists()


# ----------------------------------------------------------------------
# Every optimizer over the Model Chemistry
# ----------------------------------------------------------------------
@pytest.mark.parametrize(
    "optimizer", ["geomeTRIC/geomeTRIC", "BFGS/ASE", "LBFGS/ASE", "FIRE/ASE"]
)
def test_optimizers_over_the_model_chemistry(
    tmp_path, configuration, fake_session, optimizer
):
    node = _Structure(tmp_path, configuration)
    P, PP = _parameters(node, optimizer)
    if optimizer == "FIRE/ASE":
        P["max steps"] = "60 * natoms"
    node._working_configuration = configuration
    node.optimize(P, PP)

    assert node._results["converged"]
    assert fake_session and fake_session[0].n_calls > 2
    r01, r02, r12 = _bond_lengths(configuration)
    assert r01 == pytest.approx(0.96, abs=0.01)
    assert r02 == pytest.approx(0.96, abs=0.01)
    assert r12 == pytest.approx(1.52, abs=0.02)
    assert node._results["energy"] == pytest.approx(0.0, abs=0.1)  # kJ/mol
    # no sub-step directories: nothing ran the (empty) subflowchart
    assert not list((tmp_path / "1").glob("**/step_00001"))
    assert node._mc_session is None


def test_geometric_gradients_are_in_atomic_units(tmp_path, configuration):
    node = _Structure(tmp_path, configuration)
    node._working_configuration = configuration
    node._data = {}
    node._mc_session = SimpleNamespace(engine=FakeEngine())
    energy, gradients = node.geometric_calculate_gradients(XYZ)
    e_ev, g_ev = harmonic(XYZ)
    assert energy == pytest.approx(Q_(e_ev, "eV").m_as("E_h"))
    assert gradients == pytest.approx(Q_(g_ev, "eV/Å").m_as("E_h/a_0"))
    assert node._results["maximum_gradient"] == pytest.approx(
        Q_(np.max(np.linalg.norm(g_ev, axis=1)), "eV/Å").m_as("kJ/mol/Å")
    )
    assert node._data["step"] == [1]


def test_seamm_optimizer_over_separate_calculations(
    tmp_path, configuration, monkeypatch
):
    """With no MDI engine the SEAMM optimizer now runs over separate calculations
    (it raised 'needs a Model Chemistry that runs as an MDI engine')."""
    import contextlib

    from seamm_exec.hessian import EngineSession

    pytest.importorskip("seamm_optimizer")
    calls = []

    class Driver:
        last_optimization = None

        def say(self, text):
            pass

        def evaluate(self, evaluator, configuration, geometries, what=""):
            calls.extend(geometries)
            out = {}
            for key, x in geometries.items():
                e, g = harmonic(x)
                out[key] = (
                    Q_(e, "eV").m_as("kJ/mol"),
                    Q_(g, "eV/Å").m_as("E_h/a_0").reshape(-1),
                )
            return out

        def optimize_in_blocks(self, session, **kwargs):
            kwargs.pop("cheap_level")
            return session.optimize(**kwargs)

    driver = Driver()

    @contextlib.contextmanager
    def session(node, configuration, driver=None, print_to=None):
        engine = EvaluatorEngine(driver, None, configuration)
        yield EngineSession(driver, engine, configuration, None, [8, 1, 1], XYZ)

    monkeypatch.setattr(structure_step.structure, "model_chemistry_session", session)
    monkeypatch.setattr(structure_step.structure, "make_driver", lambda *a, **k: driver)
    node = _Structure(tmp_path, configuration)
    P, PP = _parameters(node, "SEAMM/SEAMM")
    node._working_configuration = configuration
    node.optimize(P, PP)
    assert node._results["converged"]
    assert len(calls) > 2 and all(k.startswith("opt") for k in calls)
    assert _bond_lengths(configuration)[0] == pytest.approx(0.96, abs=0.01)


# ----------------------------------------------------------------------
# Separate calculations as an engine
# ----------------------------------------------------------------------
def test_evaluator_engine_one_calculation_per_geometry(configuration):
    seen = []

    def evaluate(evaluator, configuration, geometries, what=""):
        seen.append(dict(geometries))
        ((key, x),) = geometries.items()
        return {key: (12.5, np.arange(9.0) * 0.01)}

    driver = SimpleNamespace(evaluate=evaluate)
    engine = EvaluatorEngine(driver, None, configuration)
    assert not engine.supports("<HESSIAN")
    engine.set_coordinates(XYZ, units="Å")
    assert engine.energy(units="kJ/mol") == 12.5
    forces = engine.forces(units="hartree/bohr")
    assert forces.shape == (3, 3)
    assert forces.ravel() == pytest.approx(-np.arange(9.0) * 0.01)
    assert len(seen) == 1  # energy and forces from one task
    engine.set_coordinates(XYZ, units="Å")
    engine.energy()
    assert len(seen) == 1  # same geometry: no new task
    engine.set_coordinates(Q_(XYZ + 0.01, "Å").m_as("bohr"), units="bohr")
    engine.energy()
    assert len(seen) == 2
    (key,) = seen[1]
    assert key.startswith("opt00002-")
    assert seen[1][key] == pytest.approx(XYZ + 0.01)
    # the key changes with the coordinates, so a rerun reuses only the same ones
    assert list(seen[0])[0].split("-")[1] != key.split("-")[1]


def test_session_falls_back_to_separate_calculations(configuration):
    import contextlib

    class Driver:
        last_plan = ("evaluator", None, [8, 1, 1], XYZ)
        said = []

        def model_chemistry(self):
            return FAKE_MC

        def say(self, text):
            self.said.append(text)

        @contextlib.contextmanager
        def session(self, configuration):
            yield None

    driver = Driver()
    with model_chemistry_session(None, configuration, driver=driver) as session:
        assert isinstance(session.engine, EvaluatorEngine)
        assert session.evaluator == "evaluator"
        assert session.n_atoms == 3
    assert "one separate calculation per optimization step" in driver.said[0]


def test_session_errors_name_the_problem(configuration):
    import contextlib

    class Driver:
        def model_chemistry(self):
            return FAKE_MC

        @contextlib.contextmanager
        def session(self, configuration):
            raise ValueError("no evaluator")
            yield

    with pytest.raises(RuntimeError, match="neither an MDI engine nor separate"):
        with model_chemistry_session(None, configuration, driver=Driver()):
            pass


def test_make_driver_prefers_a_warm_engine(monkeypatch):
    class Step:
        def get_mdi_engine_command(self):
            pass

    def node(options):
        return SimpleNamespace(
            get_variable=lambda name: {**FAKE_MC, "options": options},
            flowchart=SimpleNamespace(
                plugin_manager=SimpleNamespace(get=lambda name: Step())
            ),
            title="Structure",
        )

    driver = energy_source.make_driver(node({"mdi_capable": True}))
    assert driver.route == "warm engine"
    driver = energy_source.make_driver(node({}))
    assert driver.route == "separate tasks"


# ----------------------------------------------------------------------
# The flowchart check (the editor's Run, seamm-flowchart validate, MCP)
# ----------------------------------------------------------------------
def test_check_flowchart(monkeypatch):
    from seamm.problems import is_error, fixes

    node = structure_step.Structure()
    model_chemistry = SimpleNamespace(extension="Model Chemistry")
    other = SimpleNamespace(extension="FromSMILES")

    node.parameters["optimizer"].value = GEO
    problems = node.check_flowchart([other])
    assert len(problems) == 1 and NO_SOURCE in problems[0]
    assert is_error(problems[0])
    assert node.check_flowchart([model_chemistry, other]) == []

    node.parameters["optimizer"].value = SEAMM
    problems = node.check_flowchart([other])
    assert len(problems) == 1 and SEAMM_NEEDS_MODEL_CHEMISTRY in problems[0]

    # the subflowchart chosen but empty, with a Model Chemistry: a warning, and
    # the fix chooses the Model Chemistry
    node.parameters["optimizer"].value = "BFGS/ASE"
    node.parameters["engine"].value = SUB
    (problem,) = node.check_flowchart([model_chemistry])
    assert not is_error(problem)
    assert energy_source.USING_MODEL_CHEMISTRY in problem
    (fix,) = fixes(problem)
    fix.apply()
    assert node.parameters["engine"].value == MC
    assert node.check_flowchart([model_chemistry]) == []

    # sub-steps but no Model Chemistry, with the Model Chemistry chosen: a warning
    monkeypatch.setattr(structure_step.structure, "has_substeps", lambda f: True)
    (problem,) = node.check_flowchart([other])
    assert not is_error(problem) and energy_source.USING_SUBFLOWCHART in problem
    node.parameters["engine"].value = SUB
    assert node.check_flowchart([other]) == []
    node.parameters["optimizer"].value = SEAMM
    assert node.check_flowchart([other])  # sub-steps do not help it

    # nothing to optimize: nothing to check
    node.parameters["target"].value = "stereoisomers"
    assert node.check_flowchart([other]) == []


def test_validate_stops_the_flowchart_of_the_report():
    """`Model Chemistry -> from SMILES -> Structure` validates; without the Model
    Chemistry, and with no sub-steps, it does not."""
    import seamm
    from seamm.builder import validate

    def flowchart(with_mc):
        f = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
        previous = f.get_node("1")
        steps = (["Model Chemistry"] if with_mc else []) + ["Structure"]
        for name in steps:
            node = f.create_node(name)
            f.add_node(node)
            f.add_edge(previous, node, edge_type="execution", edge_subtype="next")
            previous = node
        return f

    try:
        good = flowchart(True)
    except Exception as e:  # pragma: no cover
        pytest.skip(f"the Model Chemistry step is not installed: {e}")
    assert not [p for p in validate(good) if "Structure" in p]
    problems = [p for p in validate(flowchart(False)) if "Structure" in p]
    assert problems and NO_SOURCE in problems[0]


def test_description_follows_the_choice():
    node = structure_step.Structure()
    node._id = ("1",)
    node.flowchart = SimpleNamespace(root_directory=None, executor=None)
    text = node.description_text()
    assert "come from the flowchart's Model Chemistry" in text
    node.parameters["engine"].value = SUB
    text = node.description_text()
    assert "will be calculated as follows" in text


@pytest.mark.parametrize("optimizer", ["geomeTRIC", "ASE"])
def test_substeps_without_json_results_say_what_is_missing(
    tmp_path, configuration, monkeypatch, optimizer
):
    """Sub-steps that save no energy and gradients to JSON: the error says to save
    them (it spoke of a properties.json)."""
    import seamm_ase
    import seamm_geometric

    def no_results(self, *args, **kwargs):
        raise RuntimeError(
            "There are no energy and gradients in properties.json for step 1 in x."
        )

    if optimizer == "geomeTRIC":
        monkeypatch.setattr(
            seamm_geometric.geomeTRIC_mixin, "geometric_calculate_gradients", no_results
        )
        name, args = "geometric_calculate_gradients", (XYZ,)
    else:
        monkeypatch.setattr(seamm_ase.ASE_mixin, "ase_calculator", no_results)
        name, args = "ase_calculator", (SimpleNamespace(),)
    node = _Structure(tmp_path, configuration)
    with pytest.raises(RuntimeError) as e:
        getattr(node, name)(*args)
    assert "save 'energy' and 'gradients' to JSON" in str(e.value)
    assert "properties.json" not in str(e.value)
