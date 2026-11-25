#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `structure_step` package."""

from types import SimpleNamespace
import pytest  # noqa: F401
import structure_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = structure_step.Structure()
    assert str(type(result)) == "<class 'structure_step.structure.Structure'>"


def test_convergence_preset_sets_formula_and_criteria():
    """A preset gives the formula and criteria, as the dialog would, for a
    flowchart made without it (the formula was "", the criteria the defaults)."""
    import structure_step

    node = structure_step.Structure()
    P = node.parameters.current_values_to_dict()
    P["convergence"] = "Gaussian (tight)"
    P["convergence formula"] = ""
    node.apply_convergence_preset(P)
    preset = structure_step.metadata["convergence parameters"]["Gaussian (tight)"]
    assert P["convergence formula"] == preset["convergence formula"] != ""
    value, unit = preset["RMS gradient criterion"]
    assert P["RMS gradient criterion"].m_as(unit) == float(value)

    # "custom" keeps what was set
    P["convergence"] = "custom"
    P["convergence formula"] = "E+grad"
    node.apply_convergence_preset(P)
    assert P["convergence formula"] == "E+grad"


def test_preset_gradient_criteria_are_in_atomic_units():
    """The programs' presets define the gradient criteria in E_h/bohr: Gaussian's
    maximum force is 4.5e-4 E_h/bohr = 0.0231 eV/Å (it was labeled E_h/Å, so ASE
    used 0.0123 eV/Å, and a custom criterion reached geomeTRIC 1.89x too loose)."""
    import structure_step

    node = structure_step.Structure()
    P = node.parameters.current_values_to_dict()
    P["convergence"] = "Gaussian"
    node.apply_convergence_preset(P)
    fmax = P["Maximum atomic gradient criterion"]
    assert fmax.m_as("E_h/a_0") == pytest.approx(4.5e-4)
    assert fmax.m_as("eV/Å") == pytest.approx(0.02314, rel=1e-3)
    for name, preset in structure_step.metadata["convergence parameters"].items():
        for key in ("RMS gradient criterion", "Maximum atomic gradient criterion"):
            assert preset[key][1] in ("E_h/a_0", "eV/Å", "kcal/mol/Å"), (name, key)


def test_description_with_an_empty_formula():
    """The description no longer fails with KeyError '' for a preset."""
    import structure_step

    node = structure_step.Structure()
    node.parameters["optimizer"].value = "BFGS/geomeTRIC"
    node.parameters["convergence"].value = "Gaussian (tight)"
    node.parameters["convergence formula"].value = ""
    node._id = ("1",)
    node.flowchart = SimpleNamespace(root_directory=None, executor=None)
    assert "optimized" in node.description_text(short=True)


def test_minimum_check_results():
    """The minimum check's verdict is kept (it was not before 2026.11.25)."""
    from structure_step import Structure

    last = {
        "converged": True,
        "evaluations": 6,
        "cheap hessians": 2,
        "minimum": True,
        "lowest frequency": 41.2,
        "escapes": 0,
        "confirmations": 1,
    }
    assert Structure.minimum_check_results(last) == {
        "is_minimum": True,
        "lowest_frequency": 41.2,
        "n_escapes": 0,
        "n_confirmations": 1,
        "n_cheap_hessians": 2,
    }
    # no check (no cheap level): only what there is
    assert Structure.minimum_check_results({"converged": True, "minimum": None}) == {}
    # an older seamm-exec without the confirmations
    assert Structure.minimum_check_results({"minimum": False, "escapes": 2}) == {
        "is_minimum": False,
        "n_escapes": 2,
        "n_confirmations": 0,
    }
