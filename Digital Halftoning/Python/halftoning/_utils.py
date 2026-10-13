"""Shared helpers: MATLAB-compatible image conversion and an optional JIT."""

import math
import warnings

import numpy as np

try:  # numba makes the pixel loops (error diffusion, DBS, point statistics) fast
    from numba import njit as _numba_njit

    HAVE_NUMBA = True

    def jit(func):
        return _numba_njit(cache=True)(func)

except ImportError:  # pragma: no cover - exercised only without numba
    HAVE_NUMBA = False
    warnings.warn(
        "numba is not installed: halftoning will run in pure Python and the "
        "iterative methods (EDBS, point-process statistics) will be very slow. "
        "Install it with `pip install numba`.",
        RuntimeWarning,
        stacklevel=2,
    )

    def jit(func):
        return func


# Coefficients MATLAB's rgb2gray uses for floating-point images.
_RGB2GRAY = np.array([0.298936021293775, 0.587043074451121, 0.114020904255103])


def im2double(im):
    """MATLAB ``im2double``: uint8 -> /255, uint16 -> /65535, bool -> 0/1,
    floating point unchanged (copied as float64)."""
    a = np.asarray(im)
    if a.dtype == np.uint8:
        return a.astype(np.float64) / 255.0
    if a.dtype == np.uint16:
        return a.astype(np.float64) / 65535.0
    if a.dtype == np.int16:
        return (a.astype(np.float64) + 32768.0) / 65535.0
    if a.dtype == np.bool_:
        return a.astype(np.float64)
    if np.issubdtype(a.dtype, np.floating):
        return a.astype(np.float64)
    raise TypeError(f"im2double: unsupported dtype {a.dtype}")


def rgb2gray(im):
    """MATLAB ``rgb2gray`` for a float RGB image of shape (rows, cols, 3)."""
    a = np.asarray(im, dtype=np.float64)
    return a[..., 0] * _RGB2GRAY[0] + a[..., 1] * _RGB2GRAY[1] + a[..., 2] * _RGB2GRAY[2]


def to_gray_double(im):
    """``im2double`` followed by ``rgb2gray`` for 3-channel input.

    A 4th (alpha) channel is ignored. Returns a 2-D float64 array.
    """
    a = im2double(im)
    if a.ndim == 3:
        if a.shape[2] == 1:
            a = a[..., 0]
        elif a.shape[2] in (3, 4):
            a = rgb2gray(a[..., :3])
        else:
            raise ValueError(f"expected a grayscale or RGB image, got shape {a.shape}")
    if a.ndim != 2:
        raise ValueError(f"expected a 2-D image, got shape {a.shape}")
    return a


def tile_to(mask, shape):
    """``repmat`` a threshold/class matrix over ``shape`` and crop to it."""
    mask = np.asarray(mask)
    r, c = shape
    br, bc = mask.shape
    return np.tile(mask, (math.ceil(r / br), math.ceil(c / bc)))[:r, :c]


def as_binary(X, name="X"):
    """Return a uint8 0/1 copy of a logical/0-1 point-process image."""
    a = np.asarray(X)
    if a.dtype == np.bool_:
        return a.astype(np.uint8)
    if not np.all((a == 0) | (a == 1)):
        raise ValueError(f"{name} must be logical or contain only zeros and ones")
    return a.astype(np.uint8)
