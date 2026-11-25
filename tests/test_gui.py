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


def test_energy_source_shown_and_invalid_case_in_red(tk_node, monkeypatch):
    """No Model Chemistry upstream and no sub-steps: the dialog says so, in red;
    with a Model Chemistry upstream it says where the gradients come from."""
    from structure_step.energy_source import NO_SOURCE, SEAMM_NEEDS_MODEL_CHEMISTRY

    tk_node.create_dialog()
    tk_node["optimizer"].set("geomeTRIC/geomeTRIC")
    tk_node.reset_dialog()
    label = tk_node["energy source"]
    assert _shown(tk_node, "energy source")
    assert label.cget("text") == NO_SOURCE
    assert str(label.cget("foreground")) == "red"
    assert not _tab_hidden(tk_node)  # the place to add the sub-steps

    tk_node["optimizer"].set("SEAMM/SEAMM")
    tk_node.reset_dialog()
    assert label.cget("text") == SEAMM_NEEDS_MODEL_CHEMISTRY
    assert _tab_hidden(tk_node)  # the SEAMM optimizer does not use sub-steps

    monkeypatch.setattr(type(tk_node), "_upstream_has_model_chemistry", lambda s: True)
    tk_node.reset_dialog()
    assert "Model Chemistry" in label.cget("text")
    assert str(label.cget("foreground")) == "black"
    tk_node["optimizer"].set("BFGS/ASE")
    tk_node.reset_dialog()
    assert "from the flowchart's Model Chemistry" in label.cget("text")
    assert not _tab_hidden(tk_node)

    monkeypatch.setattr(type(tk_node), "_has_substeps", lambda s: True)
    tk_node.reset_dialog()
    assert "sub-steps" in label.cget("text")
    assert "not used by this step" in label.cget("text")

    tk_node["target"].set("stereoisomers")
    tk_node.reset_dialog()
    assert not _shown(tk_node, "energy source")
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
