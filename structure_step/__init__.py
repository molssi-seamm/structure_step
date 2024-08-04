# -*- coding: utf-8 -*-

"""
structure_step
A SEAMM plug-in for energy optimized structures
"""

# Bring up the classes so that they appear to be directly in
# the structure_step package.

from .structure import Structure  # noqa: F401, E501
from .structure_parameters import StructureParameters  # noqa: F401, E501
from .structure_step import StructureStep  # noqa: F401, E501
from .tk_structure import TkStructure  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
