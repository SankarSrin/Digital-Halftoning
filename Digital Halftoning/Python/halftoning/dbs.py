"""Efficient Direct Binary Search (DBS). Python port of ``Matlab/EDBS.m``."""

import math

import numpy as np
from scipy.signal import convolve2d

from ._utils import jit, to_gray_double

HALF = 6  # half size of the 13x13 autocorrelation CPP


def gaussian_cpp(fs=7):
    """The Gaussian HVS filter of EDBS.m and its autocorrelation (13x13)."""
    d = (fs - 1) / 6
    g = (fs - 1) // 2
    GF = np.zeros((fs, fs))
    for k in range(-g, g + 1):
        for l in range(-g, g + 1):
            c = (k * k + l * l) / (2 * d * d)
            GF[k + g, l + g] = math.exp(-c) / (2 * 3.14 * d * d)  # 3.14 as in EDBS.m
    return GF, convolve2d(GF, GF)


@jit
def _dbs_kernel(dst, CEP, CPP, maxit):
    rows, cols = dst.shape
    H = HALF
    c00 = CPP[H, H]
    it = 0
    while it < maxit:
        it += 1
        count = 0
        for i in range(rows):
            for j in range(cols):
                a0c = 0.0
                a1c = 0.0
                cpx = 0
                cpy = 0
                eps_min = 0.0
                for y in range(-1, 2):
                    if i + y < 0 or i + y >= rows:
                        continue
                    for x in range(-1, 2):
                        if j + x < 0 or j + x >= cols:
                            continue
                        if y == 0 and x == 0:  # toggle
                            a0 = -1.0 if dst[i, j] == 1 else 1.0
                            a1 = 0.0
                        elif dst[i + y, j + x] != dst[i, j]:  # swap
                            a0 = -1.0 if dst[i, j] == 1 else 1.0
                            a1 = -a0
                        else:
                            a0 = 0.0
                            a1 = 0.0
                        eps = ((a0 * a0 + a1 * a1) * c00
                               + 2 * a0 * a1 * CPP[H + y, H + x]
                               + 2 * a0 * CEP[i + H, j + H]
                               + 2 * a1 * CEP[i + y + H, j + x + H])
                        if eps_min > eps:
                            eps_min = eps
                            a0c = a0
                            a1c = a1
                            cpx = x
                            cpy = y
                if eps_min < 0:
                    for y in range(-H, H + 1):
                        for x in range(-H, H + 1):
                            CEP[i + y + H, j + x + H] = CEP[i + y + H, j + x + H] + a0c * CPP[y + H, x + H]
                    for y in range(-H, H + 1):
                        for x in range(-H, H + 1):
                            CEP[i + y + cpy + H, j + x + cpx + H] = (
                                CEP[i + y + cpy + H, j + x + cpx + H] + a1c * CPP[y + H, x + H])
                    dst[i, j] = dst[i, j] + a0c
                    dst[i + cpy, j + cpx] = dst[i + cpy, j + cpx] + a1c
                    count += 1
        if count == 0:
            break
    return it


def edbs(im, init=None, max_iter=None, seed=None, return_iterations=False):
    """Efficient Direct Binary Search halftoning (``EDBS.m``).

    Parameters
    ----------
    im       : grayscale/RGB image (uint8 0..255, or float 0..1).
    init     : optional starting halftone (0/1). Default: random 50% pattern.
    max_iter : optional cap on the number of full passes (default: until no
               toggle/swap lowers the error).
    seed     : seed for the random starting pattern.

    Returns the float 0/1 halftone (and the number of passes if
    ``return_iterations``).
    """
    imr = to_gray_double(im)
    if init is None:
        rng = np.random.default_rng(seed)
        dst = (rng.random(imr.shape) > 0.5).astype(np.float64)
    else:
        dst = (np.asarray(init) != 0).astype(np.float64)
        if dst.shape != imr.shape:
            raise ValueError("init must have the same shape as the image")
    _, CPP = gaussian_cpp(7)
    CEP = convolve2d(dst - imr, CPP[::-1, ::-1])  # == MATLAB xcorr2(Err, CPP)
    maxit = np.iinfo(np.int64).max if max_iter is None else int(max_iter)
    it = _dbs_kernel(dst, CEP, CPP, maxit)
    return (dst, it) if return_iterations else dst


EDBS = edbs
