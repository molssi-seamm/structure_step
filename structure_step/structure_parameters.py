# -*- coding: utf-8 -*-
"""
Control parameters for the Structure step in a SEAMM flowchart
"""

import logging
import seamm

logger = logging.getLogger(__name__)


class StructureParameters(seamm.Parameters):
    """
    The control parameters for Structure.

    You need to replace the "time" entry in dictionary below these comments with the
    definitions of parameters to control this step. The keys are parameters for the
    current plugin,the values are dictionaries as outlined below.

    Examples
    --------
    ::

        parameters = {
            "time": {
                "default": 100.0,
                "kind": "float",
                "default_units": "ps",
                "enumeration": tuple(),
                "format_string": ".1f",
                "description": "Simulation time:",
                "help_text": ("The time to simulate in the dynamics run.")
            },
        }

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
        Each key of the dictionary is a dictionary that contains the
        the following keys:

    parameters["default"] :
        The default value of the parameter, used to reset it.

    parameters["kind"] : enum()
        Specifies the kind of a variable. One of  "integer", "float", "string",
        "boolean", or "enum"

        While the "kind" of a variable might be a numeric value, it may still have
        enumerated custom values meaningful to the user. For instance, if the parameter
        is a convergence criterion for an optimizer, custom values like "normal",
        "precise", etc, might be adequate. In addition, any parameter can be set to a
        variable of expression, indicated by having "$" as the first character in the
        field. For example, $OPTIMIZER_CONV.

    parameters["default_units"] : str
        The default units, used for resetting the value.

    parameters["enumeration"] : tuple
        A tuple of enumerated values.

    parameters["format_string"] : str
        A format string for "pretty" output.

    parameters["description"] : str
        A short string used as a prompt in the GUI.

    parameters["help_text"] : str
        A longer string to display as help for the user.

    See Also
    --------
    Structure, TkStructure, Structure StructureParameters, StructureStep
    """

    parameters = {
        "approach": {
            "default": "Optimization",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("Optimization",),
            "format_string": "",
            "description": "Approach:",
            "help_text": "The approach or method for determining the structure.",
        },
        "target": {
            "default": "minimum",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "minimum",
                "transition state",
            ),
            "format_string": "",
            "description": "Target:",
            "help_text": "The type of structure that is the target.",
        },
        "optimizer": {
            "default": "geomeTRIC/geomeTRIC",
            "kind": "enum",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Method:",
            "help_text": "The optimizer to use.",
        },
        "coordinate system": {
            "default": "TRIC: translation-rotation internal coordinates",
            "kind": "integer",
            "default_units": "",
            "enumeration": (
                "TRIC: translation-rotation internal coordinates",
                "TRIC-p: primitive (redundant) TRIC",
                "Cart: Cartesian coordinates",
                "DLC: delocalized internal coordinates",
                "HDLC: hybrid delocalized internal coordinates",
            ),
            "format_string": "",
            "description": "The coordinate system to use:",
            "help_text": "The coordinates system to use in the calculation.",
        },
        "convergence formula": {
            "default": "",
            "kind": "enum",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Convergence formula:",
            "help_text": "The formula of the criteria for convergence.",
        },
        "convergence": {
            "default": "Gaussian",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "custom",
                "ASE",
                "Gaussian (loose)",
                "Gaussian",
                "Gaussian (tight)",
                "Gaussian (very tight)",
                "Inter-fragment (tight)",
                "MolPro",
                "NWChem (loose)",
                "NWChem",
                "NWChem (tight)",
                "QChem",
                "Turbomole",
            ),
            "format_string": "s",
            "description": "Convergence criterion:",
            "help_text": "The criterion for convergence of the optimizer.",
        },
        "Energy change criterion": {
            "default": 1.0e-6,
            "kind": "float",
            "default_units": "E_h",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "Energy:",
            "help_text": "The criterion for the change of the energy between steps.",
        },
        "RMS gradient criterion": {
            "default": 3.0e-4,
            "kind": "float",
            "default_units": "E_h/Å",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "RMS gradient:",
            "help_text": "The criterion for the RMS gradient.",
        },
        "Maximum atomic gradient criterion": {
            "default": 4.5e-4,
            "kind": "float",
            "default_units": "E_h/Å",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "Maximum atomic gradient:",
            "help_text": "The criterion for the maximum gradient of an atom.",
        },
        "RMS step criterion": {
            "default": 3.0e-4,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "RMS step:",
            "help_text": "The criterion for the RMS displacement between steps.",
        },
        "Maximum atomic step criterion": {
            "default": 4.5e-4,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "Maximum atomic step:",
            "help_text": (
                "The criterion for the maximum displacement of an atom between steps."
            ),
        },
        "calculate hessian": {
            "default": "never",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("never", "first", "each", "first+last"),
            "format_string": "",
            "description": "Recalculate Hessian:",
            "help_text": "How often to recalculate the Hessian.",
        },
        "max steps": {
            "default": "12 * natoms",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("6 * natoms", "12 * natoms", "18 * natoms"),
            "format_string": "",
            "description": "Maximum # of steps:",
            "help_text": "The maximum number of steps to take.",
        },
        "continue if not converged": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Continue if not converged:",
            "help_text": "Whether to stop if the optimizer does not converge.",
        },
        "on success": {
            "default": "keep last subdirectory",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "keep last subdirectory",
                "keep all subdirectories",
                "delete all subdirectories",
            ),
            "format_string": "",
            "description": "On success:",
            "help_text": "Which subdirectories to keep.",
        },
        "on error": {
            "default": "keep all subdirectories",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "keep last subdirectory",
                "keep all subdirectories",
                "delete all subdirectories",
            ),
            "format_string": "",
            "description": "On error:",
            "help_text": "Which subdirectories to keep if there is an error.",
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("StructureParameters.__init__")

        super().__init__(
            defaults={
                **StructureParameters.parameters,
                **seamm.standard_parameters.structure_handling_parameters,
                **defaults,
            },
            data=data,
        )
