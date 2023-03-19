# -*- coding: utf-8 -*-
"""
pyddscat is a set of tools for setting up DDSCAT jobs and analaysing the results.


It includes a number of submodules:
  -core: the central definiton for a ddscat run including the run settings
  -targets: target definition classes
  -results: classes for manipulating the ddscat output files
  -ranges: classes for defining DDScat wavelength, size, and other ranges
  -utils: a handful of common utilities

"""
from . import core
from . import utils

# Create profile folder in user home directory if one does not already exist.
try:
    utils.make_profile()
except IOError:
    pass

from . import fileio
from . import ranges
from . import results
import importlib

importlib.reload(results)
from . import targets

from .core import DDscat, Settings, set_config
from .core import pol_cL, pol_cR, pol_lV, pol_lH

__all__ = [
    "DDscat",
    "Settings",
    "set_config",
    "ranges",
    "results",
    "targets",
    "fileio",
    "utils",
    "pol_cL",
    "pol_cR",
    "pol_lV",
    "pol_lH",
]

from .__about__ import __author__, __description__, __version__
