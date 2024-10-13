# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from Structure
"""

metadata = {}

# Optimizers
metadata["optimizers"] = {
    "geomeTRIC/geomeTRIC": {
        "convergence formulas": (
            "E+grad+step",
            "ASE",
            "MolPro",
            "MOPAC",
            "QChem",
        ),
        "coordinate systems": (
            "TRIC: translation-rotation internal coordinates",
            "TRIC-p: primitive (redundant) TRIC",
            "Cart: Cartesian coordinates",
            "DLC: delocalized internal coordinates",
            "HDLC: hybrid delocalized internal coordinates",
        ),
        "targets": (
            "minimum",
            "transition state",
        ),
        "default": True,
    },
    "BFGS/ASE": {
        "convergence formulas": ("ASE",),
        "targets": ("minimum",),
    },
    "BFGSLineSearch/ASE": {
        "convergence formulas": ("ASE",),
        "targets": ("minimum",),
    },
    "LBFGS/ASE": {
        "convergence formulas": ("ASE",),
        "targets": ("minimum",),
    },
    "LBFGSLineSearch/ASE": {
        "convergence formulas": ("ASE",),
        "targets": ("minimum",),
    },
    "GPMin/ASE": {
        "convergence formulas": ("ASE",),
        "targets": ("minimum",),
    },
    "MDMin/ASE": {
        "convergence formulas": ("ASE",),
        "targets": ("minimum",),
    },
    "FIRE/ASE": {
        "convergence formulas": ("ASE",),
        "targets": ("minimum",),
    },
}

# The various formulas for convergence criteria
metadata["convergence formulas"] = {
    "E+grad+step": {
        "text": (
            "  _energy change_ < Energy change criterion          \n"
            "   _rms gradient_ < RMS gradient criterion           \n"
            "_atomic gradient_ < Maximum atomic gradient criterion\n"
            "       _rms step_ < RMS step criterion               \n"
            "    _atomic step_ < Maximum atomic step criterion    "
        ),
        "criteria": (
            "Energy change criterion",
            "RMS gradient criterion",
            "Maximum atomic gradient criterion",
            "RMS step criterion",
            "Maximum atomic step criterion",
        ),
    },
    "ASE": {
        "text": "_atomic gradient_ < Maximum atomic gradient criterion",
        "criteria": ("Maximum atomic gradient criterion",),
    },
    "MolPro": {
        "text": (
            "_atomic gradient_ < Maximum atomic gradient criterion\n"
            "and\n"
            "  _energy change_ < Energy change criterion          \n"
            "or\n"
            "    _atomic step_ < Maximum atomic step criterion    "
        ),
        "criteria": (
            "Energy change criterion",
            "Maximum atomic gradient criterion",
            "Maximum atomic step criterion",
        ),
    },
    "MOPAC": {
        "text": "_rms gradient_ < RMS gradient criterion",
        "criteria": ("RMS gradient criterion",),
    },
    "QChem": {
        "text": (
            "_atomic gradient_ < Maximum atomic gradient criterion\n"
            "and\n"
            "  _energy change_ < Energy change criterion          \n"
            "or\n"
            "    _atomic step_ < Maximum atomic step criterion    "
        ),
        "criteria": (
            "Energy change criterion",
            "Maximum atomic gradient criterion",
            "Maximum atomic step criterion",
        ),
    },
}

# The convergence criteria
metadata["convergence parameters"] = {
    "ASE": {
        "Energy change criterion": ("99.9", "E_h"),
        "RMS gradient criterion": ("99.9", "E_h/Å"),
        "Maximum atomic gradient criterion": ("0.05", "eV/Å"),
        "RMS step criterion": ("99.9", "Å"),
        "Maximum atomic step criterion": ("99.9", "Å"),
        "convergence formula": "ASE",
    },
    "Gaussian (loose)": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("1.7e-4", "E_h/Å"),
        "Maximum atomic gradient criterion": ("2.5e-3", "E_h/Å"),
        "RMS step criterion": ("6.7e-3", "Å"),
        "Maximum atomic step criterion": ("1.0e-2", "Å"),
        "convergence formula": "E+grad+step",
    },
    "Gaussian": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("3.0e-4", "E_h/Å"),
        "Maximum atomic gradient criterion": ("4.5e-4", "E_h/Å"),
        "RMS step criterion": ("1.2e-3", "Å"),
        "Maximum atomic step criterion": ("1.8e-3", "Å"),
        "convergence formula": "E+grad+step",
    },
    "Gaussian (tight)": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("1.0e-6", "E_h/Å"),
        "Maximum atomic gradient criterion": ("2.0e-6", "E_h/Å"),
        "RMS step criterion": ("4.0e-5", "Å"),
        "Maximum atomic step criterion": ("6.0e-5", "Å"),
        "convergence formula": "E+grad+step",
    },
    "Gaussian (very tight)": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("1.0e-6", "E_h/Å"),
        "Maximum atomic gradient criterion": ("2.0e-6", "E_h/Å"),
        "RMS step criterion": ("4.0e-6", "Å"),
        "Maximum atomic step criterion": ("6.0e-6", "Å"),
        "convergence formula": "E+grad+step",
    },
    "Inter-fragment (tight)": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("1.0e-5", "E_h/Å"),
        "Maximum atomic gradient criterion": ("1.5e-5", "E_h/Å"),
        "RMS step criterion": ("4.0e-4", "Å"),
        "Maximum atomic step criterion": ("6.0e-4", "Å"),
        "convergence formula": "E+grad+step",
    },
    "MolPro": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("99.9", "E_h/a_0"),
        "Maximum atomic gradient criterion": ("3.0e-4", "E_h/a_0"),
        "RMS step criterion": ("99.9", "a_0"),
        "Maximum atomic step criterion": ("3.0e-4", "a_0"),
        "convergence formula": "MolPro",
    },
    "MOPAC": {
        "Energy change criterion": ("99.9", "kcal/mol"),
        "RMS gradient criterion": ("1.0", "kcal/mol/Å"),
        "Maximum atomic gradient criterion": ("99.9", "kcal/mol/Å"),
        "RMS step criterion": ("99.9", "Å"),
        "Maximum atomic step criterion": ("99.9", "Å"),
        "convergence formula": "MOPAC",
    },
    "MOPAC (precise)": {
        "Energy change criterion": ("99.9", "kcal/mol"),
        "RMS gradient criterion": ("0.01", "kcal/mol/Å"),
        "Maximum atomic gradient criterion": ("99.9", "kcal/mol/Å"),
        "RMS step criterion": ("99.9", "Å"),
        "Maximum atomic step criterion": ("99.9", "Å"),
        "convergence formula": "MOPAC",
    },
    "NWChem (loose)": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("3.0e-3", "E_h/Å"),
        "Maximum atomic gradient criterion": ("4.5e-3", "E_h/Å"),
        "RMS step criterion": ("3.6e-3", "Å"),
        "Maximum atomic step criterion": ("5.4e-3", "Å"),
        "convergence formula": "E+grad+step",
    },
    "NWChem": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("3.0e-4", "E_h/Å"),
        "Maximum atomic gradient criterion": ("4.5e-4", "E_h/Å"),
        "RMS step criterion": ("1.2e-3", "Å"),
        "Maximum atomic step criterion": ("1.8e-3", "Å"),
        "convergence formula": "E+grad+step",
    },
    "NWChem (tight)": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("1.0e-6", "E_h/Å"),
        "Maximum atomic gradient criterion": ("2.0e-6", "E_h/Å"),
        "RMS step criterion": ("4.0e-5", "Å"),
        "Maximum atomic step criterion": ("6.0e-5", "Å"),
        "convergence formula": "E+grad+step",
    },
    "QChem": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("99.9", "E_h/Å"),
        "Maximum atomic gradient criterion": ("3.0e-4", "E_h/Å"),
        "RMS step criterion": ("99.9", "Å"),
        "Maximum atomic step criterion": ("1.2e-3", "Å"),
        "convergence formula": "QChem",
    },
    "Turbomole": {
        "Energy change criterion": ("1.0e-6", "E_h"),
        "RMS gradient criterion": ("5.0e-4", "E_h/Å"),
        "Maximum atomic gradient criterion": ("1.0e-3", "E_h/Å"),
        "RMS step criterion": ("5.0e-4", "Å"),
        "Maximum atomic step criterion": ("1.0e-3", "Å"),
        "convergence formula": "E+grad+step",
    },
}

"""Description of the computational models for Structure.

Hamiltonians, approximations, and basis set or parameterizations,
only if appropriate for this code. For example::

    metadata["computational models"] = {
        "Hartree-Fock": {
            "models": {
                "PM7": {
                    "parameterizations": {
                        "PM7": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": True,
                            "code": "mopac",
                        },
                        "PM7-TS": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": False,
                            "code": "mopac",
                        },
                    },
                },
            },
        },
    }
"""
# metadata["computational models"] = {
# }

"""Description of the Structure keywords.

(Only needed if this code uses keywords)

Fields
------
description : str
    A human readable description of the keyword.
takes values : int (optional)
    Number of values the keyword takes. If missing the keyword takes no values.
default : str (optional)
    The default value(s) if the keyword takes values.
format : str (optional)
    How the keyword is formatted in the MOPAC input.

For example::
    metadata["keywords"] = {
        "0SCF": {
            "description": "Read in data, then stop",
        },
        "ALT_A": {
            "description": "In PDB files with alternative atoms, select atoms A",
            "takes values": 1,
            "default": "A",
            "format": "{}={}",
        },
    }
"""
# metadata["keywords"] = {
# }

"""Properties that Structure produces.
`metadata["results"]` describes the results that this step can produce. It is a
dictionary where the keys are the internal names of the results within this step, and
the values are a dictionary describing the result. For example::

    metadata["results"] = {
        "total_energy": {
            "calculation": [
                "energy",
                "optimization",
            ],
            "description": "The total energy",
            "dimensionality": "scalar",
            "methods": [
                "ccsd",
                "ccsd(t)",
                "dft",
                "hf",
            ],
            "property": "total energy#Psi4#{model}",
            "type": "float",
            "units": "E_h",
        },
    }

Fields
______

calculation : [str]
    Optional metadata describing what subtype of the step produces this result.
    The subtypes are completely arbitrary, but often they are types of calculations
    which is why this is name `calculation`. To use this, the step or a substep
    define `self._calculation` as a value. That value is used to select only the
    results with that value in this field.

description : str
    A human-readable description of the result.

dimensionality : str
    The dimensions of the data. The value can be "scalar" or an array definition
    of the form "[dim1, dim2,...]". Symmetric tringular matrices are denoted
    "triangular[n,n]". The dimensions can be integers, other scalar
    results, or standard parameters such as `n_atoms`. For example, '[3]',
    [3, n_atoms], or "triangular[n_aos, n_aos]".

methods : str
    Optional metadata like the `calculation` data. `methods` provides a second
    level of filtering, often used for the Hamiltionian for *ab initio* calculations
    where some properties may or may not be calculated depending on the type of
    theory.

property : str
    An optional definition of the property for storing this result. Must be one of
    the standard properties defined either in SEAMM or in this steps property
    metadata in `data/properties.csv`.

type : str
    The type of the data: string, integer, or float.

units : str
    Optional units for the result. If present, the value should be in these units.
"""
metadata["results"] = {
    "energy": {
        "description": "The total energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": "kJ/mol",
    },
    "nsteps": {
        "description": "The number of steps in the optimization",
        "dimensionality": "scalar",
        "type": "int",
        "units": "",
    },
    "converged": {
        "description": "Whether the optimization converged",
        "dimensionality": "scalar",
        "type": "bool",
        "units": "",
    },
    "rms_gradient": {
        "description": "The root mean square of the gradients",
        "dimensionality": "scalar",
        "type": "float",
        "units": "kJ/mol/Å",
    },
    "maximum_step": {
        "description": "The maximum step size of any atom",
        "dimensionality": "scalar",
        "type": "float",
        "units": "Å",
    },
    "maximum_gradient": {
        "description": "The maximum gradient on any atom",
        "dimensionality": "scalar",
        "type": "float",
        "units": "kJ/mol/Å",
    },
    "t_elapsed": {
        "description": "The elapsed time for the optimization",
        "dimensionality": "scalar",
        "type": "float",
        "units": "s",
    },
}
