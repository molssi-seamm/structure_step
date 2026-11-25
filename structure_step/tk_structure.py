# -*- coding: utf-8 -*-

"""The graphical part of a Structure step"""

import pprint  # noqa: F401
import tkinter as tk  # noqa: F401
import tkinter.ttk as ttk

from .energy_source import (
    is_model_chemistry_step,
    NO_SOURCE,
    SEAMM_NEEDS_MODEL_CHEMISTRY,
    uses_model_chemistry_only,
)
from .structure_parameters import StructureParameters
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkStructure(seamm.TkNode):
    """
    The graphical part of a Structure step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Structure_parameters.py

    See Also
    --------
    Structure, TkStructure,
    StructureParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        namespace="org.molssi.seamm.tk",
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.namespace = namespace
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )
        self.create_dialog()

    def create_dialog(self):
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the Structure_parameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkStructure.reset_dialog
        """

        frame = super().create_dialog(title="Structure", widget="notebook")
        # make it large!
        screen_w = self.dialog.winfo_screenwidth()
        screen_h = self.dialog.winfo_screenheight()
        w = int(0.9 * screen_w)
        h = int(0.8 * screen_h)
        x = int(0.05 * screen_w / 2)
        y = int(0.1 * screen_h / 2)

        self.dialog.geometry(f"{w}x{h}+{x}+{y}")

        # Add a frame for the flowchart
        notebook = self["notebook"]
        flowchart_frame = ttk.Frame(notebook)
        self["flowchart frame"] = flowchart_frame
        notebook.add(flowchart_frame, text="Flowchart", sticky=tk.NSEW)

        self.tk_subflowchart = seamm.TkFlowchart(
            master=flowchart_frame,
            flowchart=self.node.subflowchart,
            namespace=self.namespace,
        )
        self.tk_subflowchart.draw()
        # Back from the Flowchart tab, the sub-steps may have changed
        notebook.bind("<<NotebookTabChanged>>", self.reset_dialog, add="+")

        # Where the energy and gradients come from, or that there is nothing to
        # give them (red)
        self["energy source"] = ttk.Label(
            self["frame"], text="", justify=tk.LEFT, wraplength=600
        )

        # Fill in the control parameters
        # Shortcut for parameters
        P = self.node.parameters

        # structure frame to isolate widgets
        frame = self["structure frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Structure Parameters",
            labelanchor="n",
            padding=10,
        )

        for key in StructureParameters.parameters:
            if key not in ("results",):
                self[key] = P[key].widget(frame)

        # A label widget for the convergence equation
        self["convergence text"] = ttk.Label(
            frame, text="Convergence formula here", justify=tk.CENTER
        )

        # and binding to change as needed
        for key in (
            "target",
            "approach",
            "optimizer",
            "convergence formula",
            "convergence",
            "cheap level",
        ):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
        self["cheap level"].combobox.bind("<FocusOut>", self.reset_dialog)

        # structure handling frame to isolate widgets
        frame = self["handling frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Structure Handling",
            labelanchor="n",
            padding=10,
        )

        widgets = []
        row = 0
        for key in (
            "structure handling",
            "subsequent structure handling",
            "system name",
            "configuration name",
        ):
            self[key] = P[key].widget(frame)
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1
        sw.align_labels(widgets, sticky=tk.E)

        # and lay them out
        self.reset_dialog()

        self.setup_results()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        Diffusivity parameters.

        This function simply lays them out row by row with
        aligned labels. You may wish a more complicated layout that
        is controlled by values of some of the control parameters.
        If so, edit or override this method

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkDiffusivity.create_dialog
        """

        # Remove any widgets previously packed
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0

        # Where the energy and gradients come from; the Flowchart tab only when
        # the sub-steps are used
        text, ok, uses_substeps = self.energy_source_text()
        notebook = self["notebook"]
        flowchart_frame = self["flowchart frame"]
        if uses_substeps:
            notebook.add(flowchart_frame)  # shows it again if hidden
        elif notebook.select() != str(flowchart_frame):
            notebook.hide(flowchart_frame)
        if text:
            self["energy source"].configure(
                text=text, foreground="black" if ok else "red"
            )
            self["energy source"].grid(row=row, column=0, sticky=tk.W, pady=5)
            row += 1

        self["structure frame"].grid(row=row, column=0, sticky=tk.EW, pady=10)
        row += 1
        self.reset_structure_frame()

        self["handling frame"].grid(row=row, column=0, sticky=tk.EW, pady=10)
        row += 1
        self.reset_structure_handling_frame()

        return row

    def _has_substeps(self):
        """Whether the subflowchart in the Flowchart tab has a step connected
        after its start node (the graphical one: it is copied to the step only
        when the dialog closes)."""
        try:
            start = self.tk_subflowchart.get_node("1")
            return any(
                edge.edge_type == "execution"
                for edge in self.tk_subflowchart.edges(start, direction="out")
            )
        except Exception:
            return self.node.subflowchart.get_node("1").next() is not None

    def _upstream_has_model_chemistry(self):
        """Whether a Model Chemistry step comes before this one. On any error
        (e.g. the step is not yet linked into the flowchart) False, so the user
        is told what is needed."""
        try:
            return any(is_model_chemistry_step(n) for n in self.previous_nodes())
        except Exception:
            return False

    def energy_source_text(self):
        """Where the optimization's energy and gradients will come from.

        Returns
        -------
        (str, bool, bool)
            The text for the dialog ("" when nothing is optimized), whether it
            works, and whether the sub-steps (the Flowchart tab) are used.
        """
        target = self["target"].get()
        approach = self["approach"].get()
        optimizer = self["optimizer"].get()
        optimizes = (
            self.is_expr(target) or target in ("minimum", "transition state")
        ) and (self.is_expr(approach) or approach == "Optimization")
        if not optimizes:
            return "", True, False
        have_mc = self._upstream_has_model_chemistry()
        if uses_model_chemistry_only(optimizer):
            if have_mc:
                return (
                    "Energy and gradients: from the flowchart's Model Chemistry.",
                    True,
                    False,
                )
            return SEAMM_NEEDS_MODEL_CHEMISTRY, False, False
        if self._has_substeps():
            text = "Energy and gradients: from the sub-steps in the Flowchart tab"
            if have_mc:
                text += " (the flowchart's Model Chemistry is not used by this step)"
            return text + ".", True, True
        if have_mc:
            return (
                "Energy and gradients: from the flowchart's Model Chemistry (or "
                "put the calculation in the Flowchart tab to use it instead).",
                True,
                True,
            )
        return NO_SOURCE, False, True

    def reset_structure_frame(self, widget=None):
        """Layout the widgets in the structure frame
        as needed for the current state"""

        target = self["target"].get()
        approach = self["approach"].get()

        frame = self["structure frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []
        widgets2 = []

        if target == "stereoisomers":
            # Main controls
            for key in ("target", "max stereoisomers"):
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.W)
                widgets.append(self[key])
                row += 1
        elif target in ("minimum", "transition state"):
            # Main controls
            for key in ("target", "approach"):
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.W)
                widgets.append(self[key])
                row += 1

            if approach == "Optimization":
                optimizer = self["optimizer"].get()
                convergence = self["convergence"].get()

                if self.is_expr(target):
                    optimizers = [k for k, v in self.metadata["optimizers"].items()]
                else:
                    optimizers = [
                        k
                        for k, v in self.metadata["optimizers"].items()
                        if target in v["targets"]
                    ]
                self["optimizer"].config(values=optimizers)
                if optimizer not in optimizers:
                    optimizer = optimizers[0]
                    self["optimizer"].set(optimizer)

                if convergence == "custom":
                    convergence_formula = self["convergence formula"].get()
                    self["convergence formula"].state(("!disabled",))
                else:
                    convergence_formula = self.metadata["convergence parameters"][
                        convergence
                    ]["convergence formula"]
                    self["convergence formula"].state(("disabled",))

                convergence_formulas = [
                    *self.metadata["optimizers"][optimizer]["convergence formulas"]
                ] + ["custom"]
                self["convergence formula"].config(values=convergence_formulas)
                if convergence_formula not in convergence_formulas:
                    convergence_formula = convergence_formulas[0]
                self["convergence formula"].set(convergence_formula)

                self["convergence text"].config(
                    text=self.metadata["convergence formulas"][convergence_formula][
                        "text"
                    ]
                )

                # The possible convergence criteria
                convergences = [
                    k
                    for k, v in self.metadata["convergence parameters"].items()
                    if v["convergence formula"] in convergence_formulas
                ]
                self["convergence"].config(values=convergences)
                if convergence not in convergences:
                    convergence = convergences[0]

                for key in (
                    "optimizer",
                    "convergence",
                    "convergence formula",
                ):
                    self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                    widgets.append(self[key])
                    row += 1

                self["convergence text"].grid(row=row, column=1, sticky=tk.E)
                row += 1

                # Grid in the convergence criteria and make them active, or not.
                criteria = self.metadata["convergence formulas"][convergence_formula][
                    "criteria"
                ]
                for key in criteria:
                    self[key].grid(row=row, column=1, sticky=tk.EW)
                    widgets2.append(self[key])
                    row += 1
                    if convergence == "custom":
                        self[key].state(("!disabled",))
                    else:
                        self[key].state(("disabled",))

                # and if not custom, set the values and units
                if convergence != "custom":
                    for key, tmp in self.metadata["convergence parameters"][
                        convergence
                    ].items():
                        if "criterion" in key:
                            self[key].set(tmp[0], unit_string=tmp[1])

                # The coordinate system, for the optimizers that have a choice
                systems = self.metadata["optimizers"][optimizer].get(
                    "coordinate systems"
                )
                if systems:
                    self["coordinate system"].config(values=systems)
                    if self["coordinate system"].get() not in systems:
                        self["coordinate system"].set(systems[0])
                    self["coordinate system"].grid(
                        row=row, column=0, columnspan=2, sticky=tk.EW
                    )
                    widgets.append(self["coordinate system"])
                    row += 1

                if optimizer.lower().endswith("/seamm"):
                    # The SEAMM optimizer: its Hessian options; it keeps no
                    # step subdirectories and has no 'calculate hessian'
                    keys = ["cheap level"]
                    cheap = self["cheap level"].get().strip().lower()
                    if cheap in ("", "none"):
                        keys.append("model hessian")
                    else:
                        keys.append("check minimum")
                    keys += ["max steps", "continue if not converged"]
                else:
                    keys = [
                        "max steps",
                        "calculate hessian",
                        "continue if not converged",
                        "on success",
                        "on error",
                    ]
                for key in keys:
                    self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                    widgets.append(self[key])
                    row += 1

        w1 = sw.align_labels(widgets, sticky=tk.E)
        if len(widgets2) > 0:
            w2 = sw.align_labels(widgets2, sticky=tk.E)
            frame.columnconfigure(0, minsize=w1 - w2 + 50)

    def reset_structure_handling_frame(self, widget=None):
        """Layout the widgets in the structure hadnling frame
        as needed for the current state"""

        target = self["target"].get()

        frame = self["handling frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []

        if target == "stereoisomers":
            # Main controls
            for key in (
                "structure handling",
                "subsequent structure handling",
                "system name",
                "configuration name",
            ):
                self[key].grid(row=row, column=0, sticky=tk.EW)
                widgets.append(self[key])
                row += 1
        elif target in ("minimum", "transition state"):
            # Main controls
            for key in ("structure handling", "system name", "configuration name"):
                self[key].grid(row=row, column=0, sticky=tk.EW)
                widgets.append(self[key])
                row += 1
        sw.align_labels(widgets, sticky=tk.E)

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkStructure.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
