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
        self._available_elements = set()
        self._computational_models = {}
        self._base_models = {}
        self._models = {}
        self._parameterizations = {}
        self._current_base_model = None
        self._current_elements = []
        self._current_model = None
        self._current_parameterization = None

        # Get the information about the computational models
        if have_dftbplus and "computational models" in dftbplus_step.metadata:
            deep_update(self._metadata, dftbplus_step.metadata["computational models"])
        if have_mopac and "computational models" in mopac_step.metadata:
            deep_update(self._metadata, mopac_step.metadata["computational models"])

    def _change_computational_model(self, event=None, what=None):
        """Handle a change in any part of the computational model.

        Parameters
        ----------
        what : str
            What was changed: "base model", "model", or "parameterization"
        event : Object = None
            Event from the windowing system. Not used.
        """
        # print(60 * "-")
        # print(f"{event=} {what=}")
        # print(60 * "-")

        # Any elements that the user requires
        pt = self["elements"]
        required_elements = set(pt.get())

        base_model = self["base model"].get()
        model = self["model"].get()
        parameterization = self["parameterization"].get()

        if what == "base model":
            pass
        elif what == "model":
            pass
        elif what == "parameterization":
            pass
        else:
            raise RuntimeError(
                f"Cannot handle changing '{what}' in the computational model"
            )

        print(f"\t     base model = {base_model}")
        print(f"\t          model = {model}")
        print(f"\tparametrization = {parameterization}")
        print("")

        if what == "base model" and base_model != self._current_base_model:
            # The base model was changed so reset the model and parameterization
            model = "any"
            parameterization = "any"

        if what == "model" and model != self._current_model:
            # The model was changed so reset the parameterization
            parameterization = "any"

        if base_model == "any":
            base_models = set(self._base_models.keys())
            tmp_models = set()
            for values in self._base_models.values():
                tmp_models |= values["models"]
        else:
            tmp_models = set(self._base_models[base_model]["models"])

        if model != "any":
            if model not in tmp_models:
                raise RuntimeError(
                    f"Model '{model}' not in the models for base model {base_model}"
                )
            tmp_models = set([model])

        # Check for required elements

        base_models = set()
        models = set()
        parameterizations = set()

        coverage = set()
        if parameterization == "any":
            for tmp_model in tmp_models:
                for tmp_param in self._models[tmp_model]["parameterizations"]:
                    data = self._parameterizations[tmp_param]
                    if required_elements <= data["symbols"]:
                        base_models.add(data["base model"])
                        models.add(data["model"])
                        parameterizations.add(tmp_param)
                        coverage |= self._parameterizations[tmp_param]["symbols"]
        else:
            data = self._parameterizations[parameterization]
            if required_elements <= data["symbols"]:
                base_models.add(data["base model"])
                models.add(data["model"])
                parameterizations.add(parameterization)
                coverage |= self._parameterizations[parameterization]["symbols"]

        if coverage < required_elements:
            # The computational models do not cover the required elements
            computational_model = f"{base_model}/{model}/{parameterization}"
            missing = ", ".join(sorted(required_elements - coverage))
            msg = (
                f"The chosen computational model '{computational_model}' does not "
                f"cover all the elements requested. Missing are {missing}"
            )
            tk.messagebox.showwarning(title="Not all elements covered", message=msg)

            self["base model"].set(self._current_base_model)
            self["model"].set(self._current_model)
            self["parameterization"].set(self._current_parameterization)

            return

        if base_model != "any" and base_model not in base_models:
            raise RuntimeError(
                f"Base model '{model}' not in the base models. How is that possible?"
            )

        # If only one base model, use it
        if base_model == "any" and len(base_models) == 1:
            base_model = [*base_models][0]

        # If only one model, use it
        print(f"{models=}")
        if what != "model" and len(models) == 1:
            model = [*models][0]

        if what != "parameterization" and len(parameterizations) == 1:
            parameterization = [*parameterizations][0]

        # print(f"\t          models = {models}")
        # print(f"\tparametrizations = {parameterizations}")

        self["base model"].set(base_model)

        self["model"].config(values=["any", *sorted(models)])
        self["model"].set(model)

        self["parameterization"].config(values=["any", *sorted(parameterizations)])
        self["parameterization"].set(parameterization)

        pt.set_text_color("all", "black")
        pt.set_text_color(coverage, "green")

        self._current_base_model = base_model
        self._current_model = model
        self._current_parameterization = parameterization

    def _change_computational_model_sv(self, event=None, what=None):
        """Handle a change in any part of the computational model.

        Parameters
        ----------
        what : str
            What was changed: "base model", "model", or "parameterization"
        event : Object = None
            Event from the windowing system. Not used.
        """
        # print(60 * "-")
        # print(f"{event=} {what=}")
        # print(60 * "-")

        base_model = new_base_model = self._current_base_model
        model = new_model = self._current_model
        parameterization = new_parameterization = self._current_parameterization

        # Any elements that the user requires
        pt = self["elements"]
        required_elements = set(pt.get())

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
            print(f"\t     base model = {base_model}\t{new_base_model}")
            print(f"\t          model = {model}\t{new_model}")
            print(f"\tparametrization = {parameterization}\t{new_parameterization}")
            print("")

            if base_model != new_base_model:
                if base_model == "any":
                    tmp_models = [*self._base_models.keys()]
                    model = "any"
                    parameterization = "any"
                else:
                    tmp_models = self._base_models[base_model]["models"]

                if model not in tmp_models:
                    model = "any"
                if model != "any":
                    tmp_models = [model]

                # Check for required elements
                base_models = set()
                models = set()
                parameterizations = set()
                for tmp_model in tmp_models:
                    for tmp_param in self._models[tmp_model]["parameterizations"]:
                        data = self._parameterizations[tmp_param]
                        if required_elements <= data["symbols"]:
                            base_models.add(data["base model"])
                            models.add(data["model"])
                            parameterizations.add(tmp_param)

                if base_model not in base_models:
                    # If only one base model, use it
                    if len(base_models) == 1:
                        model = [*base_models][0]
                    else:
                        model = "any"
                new_base_model = base_model

                if model not in models:
                    # If only one model, use it
                    if len(models) == 1:
                        model = [*models][0]
                    else:
                        model = "any"

            if model != new_model:
                new_model = model
                if model == "any":
                    if base_model == "any":
                        tmp_models = [*self._base_models.keys()]
                    else:
                        tmp_models = self._base_models[base_model]["models"]
                    parameterization = "any"
                else:
                    tmp_models = [model]

                # Check for required elements
                base_models = set()
                models = set()
                parameterizations = set()
                for tmp_model in tmp_models:
                    for tmp_param in self._models[tmp_model]["parameterizations"]:
                        data = self._parameterizations[tmp_param]
                        if required_elements <= data["symbols"]:
                            base_models.add(data["base model"])
                            models.add(data["model"])
                            parameterizations.add(tmp_param)

                if parameterization not in parameterizations:
                    # If only one parameterization use it
                    if len(parameterizations) == 1:
                        parameterization = [*parameterizations][0]
                    else:
                        parameterization = "any"

            if parameterization != new_parameterization:
                new_parameterization = parameterization
                if parameterization != "any":
                    # Does the model need to be changed?
                    model = self._parameterizations[parameterization]["model"]

        print(f"\t     base model = {new_base_model}")
        print(f"\t          model = {new_model}")
        print(f"\tparametrization = {new_parameterization}")

        # Get the models and parameterizations given the model
        if base_model == "any":
            models = set()
            for data in self._base_models.values():
                for tmp_model in data["models"]:
                    found = False
                    for tmp in self._models[tmp_model]["parameterizations"]:
                        if required_elements <= self._parameterizations[tmp]["symbols"]:
                            found = True
                            break
                    if found:
                        models.add(tmp_model)
        else:
            models = set()
            for tmp_model in self._base_models[base_model]["models"]:
                found = False
                for tmp in self._models[tmp_model]["parameterizations"]:
                    if required_elements <= self._parameterizations[tmp]["symbols"]:
                        found = True
                        break
                if found:
                    models.add(tmp_model)

        parameterizations = set()
        if model == "any":
            for tmp_model in models:
                for tmp in self._models[tmp_model]["parameterizations"]:
                    if required_elements <= self._parameterizations[tmp]["symbols"]:
                        parameterizations.add(tmp)
        else:
            if parameterization == "any":
                for tmp in self._models[model]["parameterizations"]:
                    if required_elements <= self._parameterizations[tmp]["symbols"]:
                        parameterizations.add(tmp)
            else:
                if (
                    required_elements
                    <= self._parameterizations[parameterization]["symbols"]
                ):
                    parameterizations.add(parameterization)

        # print(f"\t          models = {models}")
        # print(f"\tparametrizations = {parameterizations}")

        # Get the elements the selection covers ... and check it meets the requirements.
        coverage = set()
        base_models = set()
        models = set()
        for tmp in parameterizations:
            base_models.add(self._parameterizations[tmp]["base model"])
            models.add(self._parameterizations[tmp]["model"])
            coverage |= self._parameterizations[tmp]["symbols"]

        if required_elements <= coverage:
            # Works, so proceed to set the widgets and current state
            self["base model"].config(values=["any", *sorted(base_models)])
            self["base model"].set(base_model)

            self["model"].config(values=["any", *sorted(models)])
            self["model"].set(model)

            self["parameterization"].config(values=["any", *sorted(parameterizations)])
            self["parameterization"].set(parameterization)

            pt.set_text_color("all", "black")
            pt.set_text_color(coverage, "green")

            self._current_base_model = base_model
            self._current_model = model
            self._current_parameterization = parameterization
        else:
            # The computational models do not cover the required elements
            computational_model = f"{base_model}/{model}/{parameterization}"
            missing = ", ".join(sorted(required_elements - coverage))
            msg = (
                f"The chosen computational model '{computational_model}' does not "
                f"cover all the elements requested. Missing are {missing}"
            )
            tk.messagebox.showwarning(title="Not all elements covered", message=msg)

            self["base model"].set(self._current_base_model)
            self["model"].set(self._current_model)
            self["parameterization"].set(self._current_parameterization)

    def _change_elements(self, event=None):
        """Handle changing the elements needed."""
        pt = self["elements"]
        elements = set(pt.get())

        base_model = self._current_base_model
        model = self._current_model
        parameterization = self._current_parameterization

        base_models = set()
        models = set()
        parameterizations = set()
        for tmp, data in self._parameterizations.items():
            if elements <= data["symbols"]:
                parameterizations.add(tmp)
                models.add(data["model"])
                base_models.add(data["base model"])

        if len(parameterizations) == 0:
            msg = (
                "No model covers all the requested elements, so backing out "
                "the last element selected."
            )
            tk.messagebox.showwarning(
                title="No model covers all the requested elements", message=msg
            )
            pt.set(self._current_elements)
            return

        if parameterization != "any" and parameterization not in parameterizations:
            msg = (
                "The currently selected model, "
                f"'{base_model}/{model}/{parameterization}' does not cover the "
                "requested elements, so am resetting the model to ones that do."
            )
            tk.messagebox.showwarning(
                title="The current model does not cover the requested elements",
                message=msg,
            )

            if base_model != "any":
                if base_model not in base_models:
                    base_model = "any"

            if model != "any":
                if model not in models:
                    model = "any"

            if parameterization != "any":
                if parameterization not in parameterizations:
                    parameterization = "any"

        # Get the models and parameterizations given the model
        if base_model == "any":
            if len(base_models) == 1:
                base_model = [*base_models][0]

        if base_model == "any":
            models = set()
            for data in self._base_models.values():
                models |= data["models"]
        else:
            models = self._base_models[base_model]["models"]

        if model == "any":
            if len(models) == 1:
                model = [*models][0]

        if model == "any":
            parameterizations = set()
            for tmp in models:
                parameterizations |= self._models[tmp]["parameterizations"]
        else:
            if parameterization == "any":
                parameterizations = self._models[model]["parameterizations"]
            else:
                parameterizations = [parameterization]

        if parameterization == "any":
            if len(parameterizations) == 1:
                parameterization = [*parameterizations][0]

        # Get the elements the selection covers ... and check it meets the requirements.
        coverage = set()
        for tmp in parameterizations:
            coverage |= self._parameterizations[tmp]["symbols"]

        # Set the widgets
        self["base model"].config(values=["any", *sorted(base_models)])
        self["base model"].set(base_model)

        self["model"].config(values=["any", *sorted(models)])
        self["model"].set(model)

        self["parameterization"].config(values=["any", *sorted(parameterizations)])
        self["parameterization"].set(parameterization)

        pt.set_text_color("all", "black")
        pt.set_text_color(coverage, "green")

        self._current_base_model = base_model
        self._current_model = model
        self._current_parameterization = parameterization

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
                            computational_model = (
                                f"{base_model}/{model}/{parameterization}"
                            )
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

        self._available_elements = available

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
        w.bind(
            "<Return>",
            lambda event, c="base model": self._change_computational_model(event, c),
        )
        w.bind(
            "<FocusOut>",
            lambda event, c="base model": self._change_computational_model(event, c),
        )

        tmp = [*self._base_models.keys()]
        w.combobox.config(values=["any", *sorted(tmp)])

        w = self["model"]
        self._current_model = w.get()
        w.bind(
            "<<ComboboxSelected>>",
            lambda event, c="model": self._change_computational_model(event, c),
        )
        w.bind(
            "<Return>",
            lambda event, c="model": self._change_computational_model(event, c),
        )
        w.bind(
            "<FocusOut>",
            lambda event, c="model": self._change_computational_model(event, c),
        )

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
        w.bind(
            "<Return>",
            lambda event, c="parameterization": self._change_computational_model(
                event, c
            ),
        )
        w.bind(
            "<FocusOut>",
            lambda event, c="parameterization": self._change_computational_model(
                event, c
            ),
        )

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
