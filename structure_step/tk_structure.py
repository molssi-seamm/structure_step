# -*- coding: utf-8 -*-

"""The graphical part of a Structure step"""

from copy import deepcopy
import pprint  # noqa: F401
import tkinter as tk
import tkinter.ttk as ttk

try:
    import dftbplus_step

    have_dftbplus = True
except Exception:
    have_dftbplus = False
try:
    import mopac_step

    have_mopac = True
except Exception:
    have_mopac = False
import structure_step  # noqa: F401
import seamm
from seamm_util import element_data
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


atno_to_symbol = {d["atomic number"]: symbol for symbol, d in element_data.items()}


def deep_update(a: dict, b: dict) -> dict:
    for bk, bv in b.items():
        av = a.get(bk)
        if isinstance(av, dict) and isinstance(bv, dict):
            a[bk] = deep_update(av, bv)
        else:
            a[bk] = deepcopy(bv)


def expand_range_list(x):
    """Expand a list of integers including ranges into a list.

    Parameters
    ----------
    x : str
        A string giving shorthand for a list, like this '1,2, 10-20, 40,50'

    Returns
    -------
    [int]
        A python list of integers
    """
    result = []
    for part in x.split(","):
        if "-" in part:
            a, b = part.split("-")
            a, b = int(a), int(b)
            result.extend(range(a, b + 1))
        elif ".." in part:
            a, b = part.split("..")
            a, b = int(a), int(b)
            result.extend(range(a, b + 1))
        else:
            a = int(part)
            result.append(a)
    return result


class TkStructure(seamm.TkNode):
    """
    The graphical part of a Structure step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    namespace : str
        The namespace of the current step.
    tk_subflowchart : TkFlowchart
        A graphical Flowchart representing a subflowchart
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

        self._metadata = {}
        self._computational_models = {}
        self._base_models = {}
        self._models = {}
        self._parameterizations = {}
        self._current_base_model = None
        self._current_elements = None
        self._current_model = None
        self._current_parameterization = None

        # Get the information about the computational models
        if have_dftbplus and "computational models" in dftbplus_step.metadata:
            deep_update(self._metadata, dftbplus_step.metadata["computational models"])
        if have_mopac and "computational models" in mopac_step.metadata:
            deep_update(self._metadata, mopac_step.metadata["computational models"])

    def _change_base_model(self, event=None):
        """Handle changing the base model, i.e. Hartree-Fock or DFT."""
        base_model = self["base model"].get()
        if base_model != self._current_base_model:
            # print(f"{base_model=}")
            self._current_base_model = base_model
            model_widget = self["model"]
            model = model_widget.get()
            if base_model == "any":
                models = [*self._base_models.keys()]
            else:
                models = self._base_models[base_model]["models"]
            model_widget.config(values=["any", *sorted(models)])
            if model not in models:
                # If only one model, use it
                if len(models) == 1:
                    model_widget.set(models[0])
                else:
                    model_widget.set("any")
                self._change_model()

    def _change_computational_model(self, event=None, what=None):
        """Handle a change in any part of the computational model.

        Parameters
        ----------
        what : str
            What was changed: "base model", "model", or "parameterization"
        event : Object = None
            Event from the windowing system. Not used.
        """
        new_base_model = self._current_base_model
        new_model = self._current_model
        new_parameterization = self._current_parameterization

        if what == "base model":
            base_model = self["base model"].get()
        elif what == "model":
            model = self["model"].get()
        elif what == "parameterization":
            parameterization = self["parameterization"].get()
        else:
            raise RuntimeError(
                f"Cannot handle changing '{what}' in the computational model"
            )

        while (
            base_model != new_base_model
            or model != new_model
            or parameterization != new_parameterization
        ):
            if base_model != new_base_model:
                # print(f"{base_model=}")
                model_widget = self["model"]
                model = model_widget.get()
                if base_model == "any":
                    models = [*self._base_models.keys()]
                else:
                    models = self._base_models[base_model]["models"]
                if model not in models:
                    # If only one model, use it
                    if len(models) == 1:
                        model = models[0]
                    else:
                        model = "any"

            if model != new_model:
                if model != new_model:
                    # print(f"{model=}")
                    new_model = model
                    widget = self["parameterization"]
                    parameterization = widget.get()
                    if model == "any":
                        base_model = self["base model"].get()
                        if base_model == "any":
                            parameterizations = set()
                            for data in self._models.values():
                                parameterizations |= data["parameterizations"]
                        else:
                            models = self._base_models[base_model]["models"]
                            parameterizations = set()
                            for tmp in models:
                                parameterizations |= self._models[tmp]["parameterizations"]
                    else:
                        parameterizations = self._models[model]["parameterizations"]
                        # Does the base model need to be changed?
                        base_model = self._models[model]["base model"]
                        if base_model != self["base model"].get():
                            self["base model"].set(base_model)
                            self._change_base_model()
                    widget.config(values=["any", *sorted(parameterizations)])
                    if parameterization not in parameterizations:
                        # If only one parameterization use it
                        if len(parameterizations) == 1:
                            widget.set(parameterizations.pop())
                        else:
                            widget.set("any")
                        self._change_parameterization()



    def _change_elements(self, event=None):
        """Handle changing the elements needed."""
        pt = self["elements"]
        elements = set(pt.get())
        print(f"change elements --> {elements}")
        # possibilities = []

    def _change_model(self, event=None):
        """Handle changing the model, i.e. PM7, MP2, CCSD, or DFT functional."""
        model = self["model"].get()
        if model != self._current_model:
            # print(f"{model=}")
            self._current_model = model
            widget = self["parameterization"]
            parameterization = widget.get()
            if model == "any":
                base_model = self["base model"].get()
                if base_model == "any":
                    parameterizations = set()
                    for data in self._models.values():
                        parameterizations |= data["parameterizations"]
                else:
                    models = self._base_models[base_model]["models"]
                    parameterizations = set()
                    for tmp in models:
                        parameterizations |= self._models[tmp]["parameterizations"]
            else:
                parameterizations = self._models[model]["parameterizations"]
                # Does the base model need to be changed?
                base_model = self._models[model]["base model"]
                if base_model != self["base model"].get():
                    self["base model"].set(base_model)
                    self._change_base_model()
            widget.config(values=["any", *sorted(parameterizations)])
            if parameterization not in parameterizations:
                # If only one parameterization use it
                if len(parameterizations) == 1:
                    widget.set(parameterizations.pop())
                else:
                    widget.set("any")
                self._change_parameterization()

    def _change_parameterization(self, event=None):
        """Handle changing the basis set or parameterization."""
        parameterization = self["parameterization"].get()
        if parameterization != self._current_parameterization:
            # print(f"{parameterization=}")
            self._current_parameterization = parameterization
            if parameterization != "any":
                # Does the model need to be changed?
                model = self._parameterizations[parameterization]["model"]
                if model != self["model"].get():
                    self["model"].set(model)
                    self._change_model()

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

        frame = super().create_dialog(title="Structure(s)")
        # Shortcut for parameters
        P = self.node.parameters

        # Create frames for the various content
        frame0 = self["element frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Elements that are needed",
            labelanchor="n",
            padding=10,
        )
        # and add the widgets
        self["elements"] = P["elements"].widget(frame0)
        self["elements"].grid(row=0, column=0, sticky=tk.EW)

        frame1 = self["main frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="What to do",
            labelanchor="n",
            padding=10,
        )
        # and add the widgets
        for key in structure_step.StructureParameters.parameters:
            if key[0] != "_" and key not in ("results", "elements"):
                self[key] = P[key].widget(frame1)

        # Show which elements are available
        available = set()
        for base_model, base_model_data in self._metadata.items():
            self._base_models[base_model] = {"models": set()}
            if "models" in base_model_data:
                for model, model_data in base_model_data["models"].items():
                    if "parameterizations" in model_data:
                        self._models[model] = {
                            "base_model": base_model,
                            "parameterizations": set(),
                        }
                        for parameterization, data in model_data[
                            "parameterizations"
                        ].items():
                            if "optimization" in data and not data["optimization"]:
                                continue
                            if "elements" not in data:
                                continue
                            elements = expand_range_list(data["elements"])
                            computational_model = f"{model}/{model}/{parameterization}"
                            tmp = self._computational_models[
                                computational_model
                            ] = deepcopy(data)
                            tmp["elements"] = elements
                            tmp["symbols"] = {atno_to_symbol[z] for z in elements}

                            self._base_models[base_model]["models"].add(model)
                            self._models[model]["parameterizations"].add(
                                parameterization
                            )
                            self._models[model]["base model"] = base_model
                            self._parameterizations[parameterization] = {
                                "base model": base_model,
                                "model": model,
                                "elements": elements,
                                "symbols": tmp["symbols"],
                                "computational_model": computational_model,
                                "data": tmp,
                            }
                            available |= tmp["symbols"]

        pt = self["elements"]
        elements = set(pt.elements)
        not_available = elements - available
        pt.disable(not_available)

        # Update the dialog as the user selects elements
        pt.command = self._change_elements

        # And as the user changes the computational model
        w = self["base model"]
        self._current_base_model = w.get()
        w.bind(
            "<<ComboboxSelected>>",
            lambda event, c="base model": self._change_computational_model(event, c),
        )
        w.bind("<Return>", self._change_base_model)
        w.bind("<FocusOut>", self._change_base_model)

        tmp = [*self._base_models.keys()]
        w.combobox.config(values=["any", *sorted(tmp)])

        w = self["model"]
        self._current_model = w.get()
        w.bind(
            "<<ComboboxSelected>>",
            lambda event, c="model": self._change_computational_model(event, c),
        )
        w.bind("<Return>", self._change_model)
        w.bind("<FocusOut>", self._change_model)

        tmp = [*self._models.keys()]
        w.combobox.config(values=["any", *sorted(tmp)])

        w = self["parameterization"]
        self._current_parameterization = w.get()
        w.bind(
            "<<ComboboxSelected>>",
            lambda event, c="parameterization": self._change_computational_model(
                event, c
            ),
        )
        w.bind("<Return>", self._change_parameterization)
        w.bind("<FocusOut>", self._change_parameterization)

        tmp = set()
        for data in self._models.values():
            tmp |= data["parameterizations"]
        w.combobox.config(values=["any", *sorted(tmp)])

        # Optimization control
        frame2 = self["optking frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Control of the optimization",
            labelanchor="n",
            padding=10,
        )
        # and add the widgets
        row = 0
        widgets = []
        for key in structure_step.StructureParameters.optking_parameters:
            self[key] = P[key].widget(frame2)
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1
        # Align the labels
        sw.align_labels(widgets, sticky=tk.E)

        # Structure handling
        frame3 = self["handling frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="How to handle the structure(s)",
            labelanchor="n",
            padding=10,
        )
        # and add the widgets
        row = 0
        widgets = []
        for key in seamm.standard_parameters.structure_handling_parameters:
            self[key] = P[key].widget(frame3)
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1
        # Align the labels
        sw.align_labels(widgets, sticky=tk.E)

        # and lay them out
        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        Structure_parameter.

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
        TkStructure.create_dialog
        """

        # Remove any widgets previously packed
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        # Shortcut for parameters
        P = self.node.parameters

        # What are we doing?
        approach = self["approach"].get()

        # Setup the periodic table
        eframe = self["element frame"]
        eframe.grid(row=0, column=0, sticky=tk.EW)

        # Setup the main frame
        main = self["main frame"]
        main.grid(row=1, column=0, sticky=tk.EW)
        for slave in main.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []
        for key in ("base model", "model", "parameterization", "approach"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1
        # Align the labels
        sw.align_labels(widgets, sticky=tk.E)

        if approach == "optimization":
            self["optimizer"].grid(row=row, column=1, sticky=tk.EW)
            row += 1
            main.columnconfigure(0, minsize=40)

            # optimizer = self["optimizer"].get()
            # if optimizer == "default" or optimizer == "OptKing":
            #     self["optking frame"].grid(row=0, column=1, sticky=tk.EW)

        self["handling frame"].grid(row=2, column=0, sticky=tk.EW)

        # Setup the results if there are any
        # NB. TO-DO may need to set self.calculation and self.methods and re-setup the
        # results.
        have_results = (
            "results" in self.node.metadata and len(self.node.metadata["results"]) > 0
        )
        if have_results and "results" in P:
            self.setup_results()

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
