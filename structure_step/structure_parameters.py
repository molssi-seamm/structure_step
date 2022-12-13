# -*- coding: utf-8 -*-
"""
Control parameters for the Structure step in a SEAMM flowchart
"""

import logging
import seamm
import pprint  # noqa: F401

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

    parameters["enumeration"]: tuple
        A tuple of enumerated values.

    parameters["format_string"]: str
        A format string for "pretty" output.

    parameters["description"]: str
        A short string used as a prompt in the GUI.

    parameters["help_text"]: str
        A longer string to display as help for the user.

    See Also
    --------
    Structure, TkStructure, Structure StructureParameters, StructureStep
    """

    parameters = {
        "approach": {
            "default": "optimization",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "optimization",
                "place holder",
            ),
            "format_string": "s",
            "description": "Approach:",
            "help_text": "The approach to generating structures.",
        },
        "elements": {
            "default": "",
            "kind": "periodic table",
            "default_units": None,
            "enumeration": None,
            "format_string": "",
            "description": "Elements:",
            "help_text": "The elements to include.",
        },
        "base model": {
            "default": "default",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Base model:",
            "help_text": "The base model for the simulation.",
        },
        "model": {
            "default": "default",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "model:",
            "help_text": "The model for the simulation.",
        },
        "basis set": {
            "default": "default",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Basis set/parameterization:",
            "help_text": "The basis set or parameterization for the simulation.",
        },
        "optimizer": {
            "default": "default",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "default",
                "native",
                "OptKing",
            ),
            "format_string": "s",
            "description": "Optimizer:",
            "help_text": "The optimizer to use.",
        },
        # Results handling
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    optking_parameters = {
        "OptKing optimization method": {
            "default": "RFO",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "RFO",
                "P_RFO",
                "NR",
                "SD",
                "LINESEARCH",
            ),
            "format_string": "s",
            "description": "Method:",
            "help_text": "The optimization method to use.",
        },
        "OptKing max geometry steps": {
            "default": "default",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "",
            "description": "Maximum steps:",
            "help_text": (
                "The maximum number of steps to take in the optimization. "
                "'default' is based on the system size, giving a reasonable "
                "limit in most cases."
            ),
        },
        "OptKing geometry convergence": {
            "default": "QCHEM",
            "kind": "float",
            "default_units": "",
            "enumeration": (
                "QCHEM",
                "MOLPRO",
                "GAU",
                "GAU_LOOSE",
                "GAU_TIGHT",
                "GAU_VERYTIGHT",
                "TURBOMOLE",
                "CFOUR",
                "NWCHEM_LOOSE",
            ),
            "format_string": "",
            "description": "Convergence criteria:",
            "help_text": "The criteria to use for convergence.",
        },
        "OptKing recalculate hessian": {
            "default": "every step",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("every step", "at beginning", "never"),
            "format_string": "",
            "description": "Recalculate Hessian:",
            "help_text": (
                "How often to recalculate the Hessian (in steps). Smaller "
                "values help convergence but are expensive."
            ),
        },
        "OptKing hessian update": {
            "default": "bfgs",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "bfgs",
                "ms",
                "powell",
                "none",
            ),
            "format_string": "s",
            "description": "Hessian update:",
            "help_text": "The algorithm for updating the Hessian.",
        },
        "OptKing coordinates": {
            "default": "Internal",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "Internal",
                "Delocalized",
                "Natural",
                "Cartesian",
                "Both",
            ),
            "format_string": "s",
            "description": "Type of coordinates:",
            "help_text": "The type of coordinates to use in the minimization.",
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
                **StructureParameters.optking_parameters,
                **seamm.standard_parameters.structure_handling_parameters,
                **defaults,
            },
            data=data,
        )
