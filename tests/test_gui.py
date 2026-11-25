# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialog: create it and re-lay it out for every choice that
drives the layout. Skipped when no display is available."""

import pytest


@pytest.fixture()
def tk_node():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw
    import seamm

    Pmw.initialise(root)
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.tk"
    )
    node = flowchart.create_node("Structure")
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get("Structure")
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    yield tk_node
    root.destroy()


def _shown(tk_node, key):
    return tk_node[key].grid_info() != {}


def test_dialog_layouts_for_every_optimizer(tk_node):
    tk_node.create_dialog()
    tk_node.reset_dialog()
    for optimizer in tk_node.metadata["optimizers"]:
        tk_node["optimizer"].set(optimizer)
        tk_node.reset_dialog()
        has_systems = "coordinate systems" in tk_node.metadata["optimizers"][optimizer]
        assert _shown(tk_node, "coordinate system") == has_systems, optimizer


def test_the_seamm_optimizer_shows_only_its_options(tk_node):
    tk_node.create_dialog()
    tk_node["optimizer"].set("SEAMM/SEAMM")
    tk_node["cheap level"].set("none")
    tk_node.reset_dialog()
    assert _shown(tk_node, "model hessian") and not _shown(tk_node, "check minimum")
    for key in ("calculate hessian", "on success", "on error"):
        assert not _shown(tk_node, key), key
    # only the max-force criterion: the convergence formula is narrowed to ASE's
    assert set(tk_node["convergence formula"].combobox.cget("values")) <= {
        "ASE",
        "custom",
    }
    tk_node["cheap level"].set("xTB:SQM@GFN2-xTB")
    tk_node.reset_dialog()
    assert _shown(tk_node, "check minimum") and not _shown(tk_node, "model hessian")


def _tab_hidden(tk_node):
    return tk_node["notebook"].tab(tk_node["flowchart frame"], "state") == "hidden"


def test_engine_choice_and_what_is_missing(tk_node, monkeypatch):
    """'Energies and gradients from:' as in the Thermochemistry step: the
    Flowchart tab and the subdirectory options only for the subflowchart, and what
    is missing for the choice in red."""
    from structure_step.energy_source import (
        NEEDS_MODEL_CHEMISTRY,
        NEEDS_SUBSTEPS,
        SEAMM_NEEDS_MODEL_CHEMISTRY,
    )

    tk_node.create_dialog()
    tk_node["optimizer"].set("geomeTRIC/geomeTRIC")
    tk_node["engine"].set("Model Chemistry")
    tk_node.reset_dialog()
    label = tk_node["energy source"]
    assert _shown(tk_node, "engine")
    assert label.cget("text") == NEEDS_MODEL_CHEMISTRY
    assert str(label.cget("foreground")) == "red"
    assert _tab_hidden(tk_node)
    assert not _shown(tk_node, "on success") and not _shown(tk_node, "on error")

    tk_node["engine"].set("subflowchart")
    tk_node.reset_dialog()
    assert label.cget("text") == NEEDS_SUBSTEPS
    assert not _tab_hidden(tk_node)  # the place to add the calculation
    assert _shown(tk_node, "on success") and _shown(tk_node, "on error")

    monkeypatch.setattr(type(tk_node), "_has_substeps", lambda s: True)
    tk_node.reset_dialog()
    assert not _shown(tk_node, "energy source")  # works as chosen

    monkeypatch.setattr(type(tk_node), "_upstream_has_model_chemistry", lambda s: True)
    tk_node["engine"].set("Model Chemistry")
    tk_node.reset_dialog()
    assert not _shown(tk_node, "energy source")
    assert _tab_hidden(tk_node)

    # The SEAMM optimizer: only the Model Chemistry, so no choice and no tab
    tk_node["optimizer"].set("SEAMM/SEAMM")
    tk_node.reset_dialog()
    assert not _shown(tk_node, "engine") and _tab_hidden(tk_node)
    assert "Model Chemistry" in label.cget("text")
    assert str(label.cget("foreground")) == "black"
    monkeypatch.setattr(type(tk_node), "_upstream_has_model_chemistry", lambda s: False)
    tk_node.reset_dialog()
    assert label.cget("text") == SEAMM_NEEDS_MODEL_CHEMISTRY

    tk_node["target"].set("stereoisomers")
    tk_node.reset_dialog()
    assert not _shown(tk_node, "energy source") and not _shown(tk_node, "engine")
    assert _tab_hidden(tk_node)


def test_upstream_model_chemistry_is_found(tk_node):
    """A Model Chemistry step connected before the Structure step in the editor."""
    tk_node.create_dialog()
    assert not tk_node._upstream_has_model_chemistry()
    tk_flowchart = tk_node.tk_flowchart
    try:
        plugin = tk_flowchart.plugin_manager.get("Model Chemistry")
    except Exception:
        pytest.skip("the Model Chemistry step is not installed")
    flowchart = tk_node.node.flowchart
    node = flowchart.create_node("Model Chemistry")
    flowchart.add_node(node)
    tk_mc = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=10
    )
    tk_flowchart.graph.add_node(tk_mc)
    tk_flowchart.add_edge(tk_mc, tk_node)
    try:
        assert tk_node._upstream_has_model_chemistry()
        assert not tk_node._has_substeps()
    finally:
        tk_flowchart.remove_node(tk_mc)
