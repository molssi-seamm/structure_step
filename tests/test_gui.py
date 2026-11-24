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
