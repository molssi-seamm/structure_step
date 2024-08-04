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
        "optimizer": {
            "default": "BFGS",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "BFGS",
                "BFGSLineSearch",
                "LBFGS",
                "LBFGSLineSearch",
                "GPMin",
                "MDMin",
                "FIRE",
            ),
            "format_string": "",
            "description": "Method:",
            "help_text": "The optimizer to use.",
        },
        "convergence": {
            "default": 0.001,
            "kind": "float",
            "default_units": "kJ/mol/Å",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "Convergence criterion:",
            "help_text": "The criterion for convergence of the optimizer.",
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
