"""Digital halftoning in Python.

Ports of the MATLAB code in this repository (with the fixes described in
REVIEW.md):

>>> from halftoning import bayers, ulichney, dotdiff, floyd_serp, errdiff, edbs
>>> from halftoning import point_process as pps

All halftoning functions take a grayscale or RGB image (uint8 0..255 or
float 0..1) and return a 0/1 array of the same size (1 = white).
"""

from ._utils import HAVE_NUMBA, im2double, rgb2gray, to_gray_double
from .ordered import Bayers, Ulichney, bayers, ulichney, ordered_dither
from .error_diffusion import FLOYD_STEINBERG, JARVIS, STUCKI, errdiff, floyd_serp
from .dot_diffusion import dotdiff
from .dbs import EDBS, edbs
from . import point_process

__all__ = [
    "bayers", "ulichney", "ordered_dither", "Bayers", "Ulichney",
    "errdiff", "floyd_serp", "FLOYD_STEINBERG", "JARVIS", "STUCKI",
    "dotdiff", "edbs", "EDBS", "point_process",
    "im2double", "rgb2gray", "to_gray_double", "HAVE_NUMBA",
]
__version__ = "1.0.0"
