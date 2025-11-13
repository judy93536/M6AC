"""
M6AC - MODTRAN6-based Atmospheric Correction for AVIRIS-NG

A physics-based atmospheric correction method combining proven retrieval
algorithms with MODTRAN6's native spectral response function handling.
"""

__version__ = '0.1.0'
__author__ = 'Judy Northrop'

from .core import M6AC
from .retrieval import water_vapor_cibr, aerosol_ddv
from .lut import LUTGenerator, LUTInterpolator
from .utils import extract_scene_metadata

__all__ = [
    'M6AC',
    'water_vapor_cibr',
    'aerosol_ddv',
    'LUTGenerator',
    'LUTInterpolator',
    'extract_scene_metadata',
]
