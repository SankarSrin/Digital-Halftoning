"""Ordered dithering: Bayer and Ulichney threshold arrays.

Python ports of ``Matlab/Bayers.m`` and ``Matlab/Ulichney.m``.
"""

import numpy as np

from ._utils import tile_to, to_gray_double

# Bayer dispersed dot (Matlab/Bayers.m, met=1)
BAYER_DISPERSED = np.array([
    [0, 48, 12, 60, 3, 51, 15, 63],
    [32, 16, 44, 28, 35, 19, 47, 31],
    [8, 56, 4, 52, 11, 59, 7, 55],
    [40, 24, 36, 20, 43, 27, 39, 23],
    [2, 50, 14, 62, 1, 49, 13, 61],
    [34, 18, 46, 30, 33, 17, 45, 29],
    [10, 58, 6, 54, 9, 57, 5, 53],
    [42, 26, 38, 22, 41, 25, 37, 21],
]) / 64

# Bayer clustered dot (Matlab/Bayers.m, met=2)
BAYER_CLUSTERED = np.array([
    [24, 10, 12, 26, 35, 47, 49, 37],
    [8, 0, 2, 14, 45, 59, 61, 51],
    [22, 6, 4, 16, 43, 57, 63, 53],
    [30, 20, 18, 28, 33, 41, 55, 39],
    [34, 46, 48, 36, 25, 11, 13, 27],
    [44, 58, 60, 50, 9, 1, 3, 15],
    [42, 56, 62, 52, 23, 7, 5, 17],
    [32, 40, 54, 38, 31, 21, 19, 29],
]) / 64

# Ulichney 4x4 clustered dot (met=1). Thresholds 0..15/16 (was /15: white never printed fully white).
ULICHNEY_CLUSTERED_4 = np.array([
    [12, 5, 6, 13],
    [4, 0, 1, 7],
    [11, 3, 2, 8],
    [15, 10, 9, 14],
]) / 16

# Ulichney 8x8 45-degree clustered dot (met=2)
ULICHNEY_CLUSTERED_8 = np.array([
    [0.567, 0.635, 0.608, 0.514, 0.424, 0.365, 0.392, 0.486],
    [0.847, 0.878, 0.910, 0.698, 0.153, 0.122, 0.090, 0.302],
    [0.820, 0.969, 0.941, 0.667, 0.180, 0.031, 0.059, 0.333],
    [0.725, 0.788, 0.757, 0.545, 0.275, 0.212, 0.243, 0.455],
    [0.424, 0.365, 0.392, 0.486, 0.567, 0.635, 0.608, 0.514],
    [0.153, 0.122, 0.090, 0.302, 0.847, 0.878, 0.910, 0.698],
    [0.180, 0.031, 0.059, 0.333, 0.820, 0.969, 0.941, 0.667],
    [0.275, 0.212, 0.243, 0.455, 0.725, 0.788, 0.757, 0.545],
])

# Ulichney 14x14 void-and-cluster dispersed dot (met=3)
ULICHNEY_VOID_CLUSTER_14 = np.array([
    [131, 187, 8, 78, 50, 18, 134, 89, 155, 102, 29, 95, 184, 73],
    [22, 86, 113, 171, 142, 105, 34, 166, 9, 60, 151, 128, 40, 110],
    [168, 137, 45, 28, 64, 188, 82, 54, 124, 189, 80, 13, 156, 56],
    [7, 61, 186, 121, 154, 6, 108, 177, 24, 100, 38, 176, 93, 123],
    [83, 148, 96, 17, 88, 133, 44, 145, 69, 161, 139, 72, 30, 181],
    [115, 27, 163, 47, 178, 65, 164, 14, 120, 48, 5, 127, 153, 52],
    [190, 58, 126, 81, 116, 21, 106, 77, 173, 92, 191, 63, 99, 12],
    [76, 144, 4, 185, 37, 149, 192, 39, 135, 23, 117, 31, 170, 132],
    [35, 172, 103, 66, 129, 79, 3, 97, 57, 159, 70, 141, 53, 94],
    [114, 20, 49, 158, 19, 146, 169, 122, 183, 11, 104, 180, 2, 165],
    [152, 87, 182, 118, 91, 42, 67, 25, 84, 147, 43, 85, 125, 68],
    [16, 136, 71, 10, 193, 112, 160, 138, 51, 111, 162, 26, 194, 46],
    [174, 107, 41, 143, 33, 74, 1, 101, 195, 15, 75, 140, 109, 90],
    [32, 62, 157, 98, 167, 119, 179, 59, 36, 130, 175, 55, 0, 150],
]) / 196


def ordered_dither(im, threshold):
    """Threshold a grayscale image against a tiled threshold array.

    Returns a boolean halftone (True = white / 1), same size as ``im``.
    """
    g = to_gray_double(im)
    return g > tile_to(threshold, g.shape)


def bayers(im, met=1):
    """Bayer ordered dither (``Bayers.m``).

    met=1: Bayer dispersed dot, met=2 (or anything else): Bayer clustered dot.
    ``im`` may be uint8/uint16/float (0..1), grayscale or RGB, any size.
    """
    return ordered_dither(im, BAYER_DISPERSED if met == 1 else BAYER_CLUSTERED)


def ulichney(im, met=1):
    """Ulichney ordered dither (``Ulichney.m``).

    met=1: 4x4 clustered dot, met=2: 8x8 clustered dot,
    met=3: 14x14 void-and-cluster dispersed dot.
    """
    masks = {1: ULICHNEY_CLUSTERED_4, 2: ULICHNEY_CLUSTERED_8, 3: ULICHNEY_VOID_CLUSTER_14}
    if met not in masks:
        raise ValueError("ulichney: met must be 1, 2 or 3")
    return ordered_dither(im, masks[met])


# MATLAB-style aliases
Bayers = bayers
Ulichney = ulichney
