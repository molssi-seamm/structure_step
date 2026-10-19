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
