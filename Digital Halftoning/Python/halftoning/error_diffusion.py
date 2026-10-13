"""Error diffusion. Python ports of ``Matlab/errdiff.m`` and ``floyd_serp.m``."""

import numpy as np

from ._utils import jit, to_gray_double

ORIGIN = -99  # marker used inside error filters for the current pixel

FLOYD_STEINBERG = np.array([[0, ORIGIN * 16, 7], [3, 5, 1]]) / 16
JARVIS = np.array([[0, 0, ORIGIN * 48, 7, 5],
                   [3, 5, 7, 5, 3],
                   [1, 3, 5, 3, 1]]) / 48
STUCKI = np.array([[0, 0, ORIGIN * 42, 8, 4],
                   [2, 4, 8, 4, 2],
                   [1, 2, 4, 2, 1]]) / 42


@jit
def _errdiff_kernel(inp, fc, r0, c0, lam, serpentine):
    ri, ci = inp.shape
    rm = fc.shape[0] - 1
    cm = fc.shape[1] - 1
    # Pad by the full filter width on both sides so the filter can be
    # mirrored for serpentine scanning wherever its origin column is.
    inpad = np.zeros((ri + rm, ci + 2 * cm))
    inpad[r0:r0 + ri, cm:cm + ci] = inp
    out = np.zeros((ri, ci))
    qn = np.zeros((ri, ci))
    f = fc.copy()
    c = c0
    left_to_right = True
    for y in range(ri):
        for t in range(ci):
            x = t if left_to_right else ci - 1 - t
            inpix = inpad[y + r0, x + cm]
            outpix = 1.0 if (inpix + lam * inp[y, x]) >= 0.5 else 0.0
            out[y, x] = outpix
            qerr = outpix - inpix
            qn[y, x] = qerr
            cs = x + cm - c
            for a in range(rm + 1):
                for b in range(cm + 1):
                    inpad[y + a, cs + b] = inpad[y + a, cs + b] - qerr * f[a, b]
        if serpentine:
            left_to_right = not left_to_right
            f = f[:, ::-1].copy()
            c = cm - c
    return out, qn


def errdiff(in_, fc=None, l=0.0, dir=1, v=False):
    """Error diffusion halftoning (``errdiff.m``).

    Parameters
    ----------
    in_ : 2-D array, values 0..1 (not rescaled, as in the MATLAB code).
    fc  : error filter; the current pixel is marked with ``-99`` (times the
          normaliser), e.g. ``FLOYD_STEINBERG`` (the default).
    l   : threshold modulation parameter (Eschbach); 0 = standard.
    dir : 1 = raster scan, -1 = serpentine.
    v   : ignored (kept for signature compatibility).

    Returns
    -------
    out : binary halftone (float 0/1)
    qn  : quantisation error per pixel
    k   : the ``k`` statistic of errdiff.m
    """
    inp = np.asarray(in_, dtype=np.float64)
    if inp.ndim != 2:
        raise ValueError("errdiff expects a 2-D array; convert RGB to gray first")
    fc = FLOYD_STEINBERG if fc is None else fc
    fc = np.array(fc, dtype=np.float64)
    # MATLAB: find(fc==-99) after the filter was scaled, i.e. the origin is
    # the entry whose value equals -99 exactly.
    pos = np.argwhere(fc == ORIGIN)
    if len(pos) != 1:
        raise ValueError("error filter must contain exactly one -99 origin marker")
    r0, c0 = (int(p) for p in pos[0])
    fc[r0, c0] = 0.0
    out, qn = _errdiff_kernel(inp, fc, r0, c0, float(l), dir == -1)
    xp = (out - 0.5 - qn).ravel(order="F")
    k = np.sum(np.abs(xp)) / (2 * np.sum(xp ** 2))
    return out, qn, k


def floyd_serp(im):
    """Floyd-Steinberg error diffusion with serpentine scan (``floyd_serp.m``).

    Accepts uint8/float/RGB input. Returns ``(halftone, qerr, k)``.
    """
    return errdiff(to_gray_double(im), FLOYD_STEINBERG, 0, -1, False)
