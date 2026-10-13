"""Point-process statistics for halftones (Lau & Arce, *Modern Digital Halftoning*).

Python ports of the ``Point Process Statistics`` folder:

=============  ===========================================================
``pp``         pp.m      - minority-pixel point process of a halftone
``intnsty``    intnsty.m - intensity (fraction of points) per channel
``bayer``      bayer.c   - Lau's 8x8 Bayer dither
``rapsd``      rapsd.c   - radially averaged power spectral density
``anis``       anis.c    - anisotropy (dB)
``pc``         pc.c      - pair correlation function
``vpc``        vpc.c     - (cross-)pair correlation between channels
``rsmm``       rsmm.c    - reduced second moment measure
``ddf``        ddf.c     - directional distribution function
=============  ===========================================================

Each function reproduces the corrected C code exactly (same sort order, same
bin edges, same normalisation); see REVIEW.md for what was corrected.
Images are indexed ``X[row, col]`` like MATLAB.
"""

import math

import numpy as np
from scipy import ndimage

from ._utils import as_binary, jit

PI_C = 3.14159265358979  # the constant #defined in the C sources


# --------------------------------------------------------------------------
# pp.m / intnsty.m
# --------------------------------------------------------------------------
def pp(X):
    """Point process of a binary halftone: the minority pixels of each channel.

    Channels with more than half ones are inverted. Returns a bool array.
    """
    a = np.asarray(X)
    if not np.all((a == 0) | (a == 1)):
        raise ValueError("X must contain only zeros and ones!")
    a = a != 0
    if a.ndim == 2:
        return a.copy() if a.mean() <= 0.5 else ~a
    out = np.empty_like(a)
    for c in range(a.shape[2]):
        out[..., c] = a[..., c] if a[..., c].mean() <= 0.5 else ~a[..., c]
    return out


def intnsty(X):
    """Intensity (fraction of points) of each channel; returns shape (C,)."""
    a = np.asarray(X, dtype=np.float64)
    if a.ndim == 2:
        a = a[..., None]
    M, N = a.shape[:2]
    return np.array([a[..., c].sum() / (M * N) for c in range(a.shape[2])])


# --------------------------------------------------------------------------
# bayer.c
# --------------------------------------------------------------------------
_LAU_BAYER = np.array([
    0.0234, 0.5234, 0.1484, 0.6484, 0.0547, 0.5547, 0.1797, 0.6797,
    0.9297, 0.2734, 0.7734, 0.3984, 0.8984, 0.3047, 0.8047, 0.4297,
    0.2422, 0.7422, 0.0859, 0.5859, 0.2109, 0.7109, 0.1172, 0.6172,
    0.8672, 0.4922, 0.9922, 0.3359, 0.8359, 0.4609, 0.9609, 0.3672,
    0.0391, 0.5391, 0.1641, 0.6641, 0.0078, 0.5078, 0.1328, 0.6328,
    0.8828, 0.2891, 0.7891, 0.4141, 0.9141, 0.2578, 0.7578, 0.3828,
    0.1953, 0.6953, 0.1016, 0.6016, 0.2266, 0.7266, 0.0703, 0.5703,
    0.8203, 0.4453, 0.9453, 0.3516, 0.8516, 0.4766, 0.9766, 0.3203,
]).reshape(8, 8).T  # the C array is column-major: mask[m%8 + 8*(n%8)]


def bayer(X):
    """Lau's ``bayer.c``: 8x8 Bayer dither of a double image in 0..1 (bool out)."""
    x = np.asarray(X, dtype=np.float64)
    if x.ndim != 2:
        raise ValueError("bayer: X must be a 2-D double image")
    r, c = x.shape
    mask = np.tile(_LAU_BAYER, (-(-r // 8), -(-c // 8)))[:r, :c]
    return x > mask


# --------------------------------------------------------------------------
# rapsd.c / anis.c
# --------------------------------------------------------------------------
def _spectral_args(X, M, FM, shift, rapsd_rule):
    x = np.asarray(X, dtype=np.float64)
    if x.ndim != 2:
        raise ValueError("X must be a 2-D image")
    spat = int(M)
    if spat == 0:
        raise ValueError("Spatial Window Dimensions cannot be 0")
    spec = spat
    sr = sc = spat
    if FM is not None:
        spec = int(FM)
        if spec == 0 or (rapsd_rule and spec < 32):  # rapsd ignores FM < 32
            spec = spat
    if shift is not None:
        s = np.atleast_1d(np.asarray(shift, dtype=np.float64)).ravel()
        if s.size == 1:
            sr = sc = int(s[0])
        elif s.size == 2:
            sr, sc = int(s[0]), int(s[1])
        else:
            raise ValueError("[DM DN] must have 1 or 2 elements")
        if sr * sc == 0:
            sr = sc = spat
    if spec < spat:
        raise ValueError("FM must be at least M (spectral window >= spatial window)")
    if spat > x.shape[0] or spat > x.shape[1]:
        raise ValueError("Spatial window M is larger than the image")
    fr_max = int(math.floor(spec / math.sqrt(2)))
    if spec % 2 == 1:
        fr_max = int(math.floor((spec - 1) / math.sqrt(2)))
    return x, spat, spec, sr, sc, fr_max


def _power_estimate(x, spat, spec, sr, sc):
    rows, cols = x.shape
    P = np.zeros((spec, spec))
    win = np.zeros((spec, spec))
    count = 0
    m = 0
    while m + spat <= rows:
        n = 0
        while n + spat <= cols:
            count += 1
            win[:spat, :spat] = x[m:m + spat, n:n + spat]
            F = np.fft.fft2(win)
            P += (F.real * F.real + F.imag * F.imag) / (spec * spec)
            n += sc
        m += sr
    return P / count, count


def _radial_index(spec, start=0):
    i = np.arange(start, spec)[:, None].astype(np.float64)
    j = np.arange(start, spec)[None, :].astype(np.float64)
    f = np.minimum.reduce([
        np.sqrt(i * i + j * j),
        np.sqrt(i * i + (j - spec) ** 2),
        np.sqrt((i - spec) ** 2 + (j - spec) ** 2),
        np.sqrt((i - spec) ** 2 + j * j),
    ])
    return (np.floor(f) - 1).astype(np.int64)


def rapsd(X, M, FM=None, shift=None):
    """Radially averaged power spectral density (``rapsd.c``).

    ``P, fr, count = rapsd(X, M, FM=None, shift=None)``

    X     : image (double; a 0/1 halftone)
    M     : spatial window size; FM : spectral (zero-padded) window size,
            ignored if < 32 as in the C code; shift : scalar or [DM, DN]
            window displacement (default M, i.e. non-overlapping).
    Returns P normalised by g(1-g), the radial frequencies fr (cycles/pixel)
    and the number of periodograms averaged.
    """
    x, spat, spec, sr, sc, fr_max = _spectral_args(X, M, FM, shift, True)
    P, count = _power_estimate(x, spat, spec, sr, sc)
    g = x.ravel(order="F").sum() / x.size
    variance = g * (1 - g)
    ri = _radial_index(spec).ravel()
    ok = ri >= 0
    annulus = np.zeros(fr_max)
    samples = np.zeros(fr_max)
    np.add.at(annulus, ri[ok], P.ravel()[ok])
    np.add.at(samples, ri[ok], 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        annulus = annulus / (samples * variance)
    fr = (0.5 + np.arange(fr_max)) / spec
    return annulus, fr, count


def anis(X, M, FM=None, shift=None):
    """Anisotropy in dB (``anis.c``). Same arguments/returns as :func:`rapsd`."""
    x, spat, spec, sr, sc, fr_max = _spectral_args(X, M, FM, shift, False)
    P, count = _power_estimate(x, spat, spec, sr, sc)
    ri = _radial_index(spec).ravel()
    ok = ri >= 0
    annulus = np.zeros(fr_max)
    samples = np.zeros(fr_max)
    np.add.at(annulus, ri[ok], P.ravel()[ok])
    np.add.at(samples, ri[ok], 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        annulus = annulus / samples
        # variance loop starts at i=1, j=1 as in anis.c
        ri1 = _radial_index(spec, start=1).ravel()
        P1 = P[1:, 1:].ravel()
        ok1 = ri1 >= 0
        A = np.zeros(fr_max)
        d = P1[ok1] - annulus[ri1[ok1]]
        np.add.at(A, ri1[ok1], d * d)
        A = A / (annulus * annulus * (samples - 1))
        A = 10 * np.log10(A)
    fr = (0.5 + np.arange(fr_max)) / spec
    return A, fr, count


# --------------------------------------------------------------------------
# Sorted offset ("cost") array shared by pc / vpc / rsmm / ddf
# --------------------------------------------------------------------------
@jit
def _hpsort_perm(cost):
    """Numerical Recipes heapsort, exactly as in the C files, on an index
    permutation (ties keep the C code's order)."""
    n = cost.shape[0]
    ra = np.arange(n)
    if n < 2:
        return ra
    l = (n >> 1) + 1
    ir = n
    while True:
        if l > 1:
            l -= 1
            rra = ra[l - 1]
        else:
            rra = ra[ir - 1]
            ra[ir - 1] = ra[0]
            ir -= 1
            if ir == 1:
                ra[0] = rra
                break
        i = l
        j = l + l
        while j <= ir:
            if j < ir and cost[ra[j - 1]] < cost[ra[j]]:
                j += 1
            if cost[rra] < cost[ra[j - 1]]:
                ra[i - 1] = ra[j - 1]
                i = j
                j <<= 1
            else:
                j = ir + 1
        ra[i - 1] = rra
    return ra


_COST_CACHE = {}


def _cost_array(rows, cols):
    """Offsets (row, col) of an rows x cols window sorted by distance."""
    key = (rows, cols)
    if key not in _COST_CACHE:
        m0 = -(rows // 2) + 1
        n0 = -(cols // 2) + 1
        jj, ii = np.meshgrid(np.arange(cols), np.arange(rows))
        r = (m0 + ii).ravel(order="F").astype(np.int64)  # column-major like C
        c = (n0 + jj).ravel(order="F").astype(np.int64)
        cost = np.sqrt((r * r + c * c).astype(np.float64))
        p = _hpsort_perm(cost)
        _COST_CACHE[key] = (r[p], c[p], cost[p])
        if len(_COST_CACHE) > 8:
            _COST_CACHE.pop(next(iter(_COST_CACHE)))
    return _COST_CACHE[key]


def _bin_index(cost, dr):
    q = np.ceil(cost / dr).astype(np.int64) - 1
    q[q == -1] = 0
    return q


# --------------------------------------------------------------------------
# rsmm.c
# --------------------------------------------------------------------------
@jit
def _second_moment_counts(Xc, Xn, orow, ocol, ocost, maxd, start):
    rows, cols = Xc.shape
    npts = np.zeros(orow.shape[0], dtype=np.int64)
    nmin = np.zeros(orow.shape[0], dtype=np.int64)
    for n in range(cols):
        for m in range(rows):
            if Xc[m, n] == 1:
                for r in range(start, orow.shape[0]):
                    if ocost[r] > maxd:
                        break
                    if (orow[r] >= m or orow[r] >= rows - m
                            or ocol[r] >= n or ocol[r] >= cols - n
                            or m + orow[r] < 0 or n + ocol[r] < 0):
                        break
                    if ocost[r] >= maxd:
                        break
                    npts[r] += 1
                    if Xn[m + orow[r], n + ocol[r]] == 1:
                        nmin[r] += 1
    return npts, nmin


def rsmm(X, max_distance=None):
    """Reduced second moment measure (``rsmm.c``).

    Returns an array the size of X; element ``[R//2-1 + dy, C//2-1 + dx]`` is
    the measure at offset (dy, dx). Offsets that were never sampled are NaN.
    """
    x = as_binary(X)
    if x.ndim != 2:
        raise ValueError("rsmm: X must be 2-D")
    rows, cols = x.shape
    maxd = float(rows if max_distance is None else max_distance)
    orow, ocol, ocost = _cost_array(rows, cols)
    npts, nmin = _second_moment_counts(x, x, orow, ocol, ocost, maxd, 1)
    gray = float(x.sum()) / (float(rows) * float(cols))
    out = np.zeros((rows, cols))
    xr = orow + rows // 2 - 1
    yc = ocol + cols // 2 - 1
    ok = (xr < rows) & (yc < cols)
    with np.errstate(divide="ignore", invalid="ignore"):
        val = nmin.astype(np.float64) / (npts.astype(np.float64) * gray)
    out[xr[ok], yc[ok]] = val[ok]
    return out


# --------------------------------------------------------------------------
# vpc.c
# --------------------------------------------------------------------------
def vpc(X, C=(1, 1), dr=0.5, max_distance=None):
    """Pair correlation between channels C=[I, J] of a (multi-channel)
    point process (``vpc.c``). Points of channel J are the centres; points of
    channel I are counted around them.

    Returns ``(g, r)``: correlation per radial bin and the bin centres.
    """
    a = as_binary(X)
    if a.ndim == 2:
        a = a[..., None]
    if a.ndim != 3:
        raise ValueError("Input image X must be 2 or 3 dimensional!")
    rows, cols, chn = a.shape
    ci, cj = int(C[0]) - 1, int(C[1]) - 1
    if not (0 <= ci < chn and 0 <= cj < chn):
        raise ValueError("Non-existing channels specified by C!")
    maxd = float(rows // 2) if max_distance is None else float(max_distance)
    if maxd > cols:
        maxd = float(cols // 2)
    if maxd > rows:
        maxd = float(rows // 2)
    nb = int(math.ceil(maxd / dr)) - 1
    orow, ocol, ocost = _cost_array(rows, cols)
    q = _bin_index(ocost, dr)
    Xi = np.ascontiguousarray(a[..., ci])
    Xj = np.ascontiguousarray(a[..., cj])
    npts, nmin = _second_moment_counts(Xj, Xi, orow, ocol, ocost, maxd, 0)
    gray = float(Xi.sum()) / (float(rows) * float(cols))
    stop = np.argmax(q >= nb) if np.any(q >= nb) else q.size
    g = np.zeros(max(nb, 0))
    cnt = np.zeros(max(nb, 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        term = nmin[:stop].astype(np.float64) / (gray * npts[:stop].astype(np.float64))
        np.add.at(g, q[:stop], term)  # sequential, in sorted order, like C
        np.add.at(cnt, q[:stop], 1.0)
        g = g / cnt
    r = (np.arange(max(nb, 0)) + 0.5) * dr
    return g, r


# --------------------------------------------------------------------------
# pc.c
# --------------------------------------------------------------------------
@jit
def _pc_counts(prow, pcol, crow, ccol, nb, dR, rows, cols):
    hist = np.zeros(nb)
    npts = np.zeros(nb)
    for k in range(crow.shape[0]):
        Mc = crow[k]
        Nc = ccol[k]
        maxd = nb * dR
        if maxd > Mc:
            maxd = Mc
        if maxd > Nc:
            maxd = Nc
        if maxd > cols - Nc:
            maxd = cols - Nc
        if maxd > rows - Mc:
            maxd = rows - Mc
        q = int(math.ceil(maxd / dR)) - 1
        if q == -1:
            q = 0
        for n in range(q):
            npts[n] += 1
        maxd = dR * q
        for m in range(prow.shape[0]):
            dy = prow[m] - Mc
            dx = pcol[m] - Nc
            cost = math.sqrt(dy * dy + dx * dx)
            if cost <= maxd:
                b = int(math.ceil(cost / dR)) - 1
                if b == -1:
                    b = 0
                hist[b] += 1
    return hist, npts


def pc(X, dr, max_distance=None, mode=None):
    """Pair correlation function (``pc.c``).

    ``g, r = pc(X, dr, max_distance=None, mode=None)``

    X    : logical point process (e.g. ``pp(halftone)``)
    dr   : radial bin width; max_distance : largest radius
    mode : ``'cluster'`` (anything starting with 'clus') uses the centroids of
           4-connected clusters as centres; otherwise every point.
    """
    x = as_binary(X)
    if x.ndim != 2:
        raise ValueError("pc: X must be 2-D")
    rows, cols = x.shape
    dR = float(dr)
    if max_distance is not None:
        maxd = float(max_distance)
        if rows < cols and (rows // 2) < maxd:
            nb = int(math.ceil(rows / (2.0 * dR))) - 1
        elif cols <= rows and (cols // 2) < maxd:
            nb = int(math.ceil(cols / (2.0 * dR))) - 1
        else:
            nb = int(math.ceil(maxd / dR)) - 1
    else:
        if rows < cols:
            nb = int(math.ceil(rows / (2.0 * dR))) - 1
        else:
            nb = int(math.ceil(cols / (2.0 * dR))) - 1
    nb = max(nb, 0)
    r = (np.arange(nb) + 0.5) * dR

    # area (number of lattice offsets) in each annulus
    _, _, ocost = _cost_array(rows, cols)
    q = _bin_index(ocost, dR)
    volume = np.bincount(q[q < nb], minlength=nb).astype(np.float64)[:nb]

    col_idx, row_idx = np.nonzero(x.T)  # point list in column-major order
    prow = row_idx.astype(np.float64)
    pcol = col_idx.astype(np.float64)
    intensity = float(prow.size) / float(rows) / float(cols)

    if isinstance(mode, str) and mode[:4] == "clus":
        lab, ncl = ndimage.label(x, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
        ii, jj = np.nonzero(lab)
        lv = lab[ii, jj]
        size = np.bincount(lv, minlength=ncl + 1)[1:].astype(np.float64)
        crow = np.bincount(lv, weights=ii, minlength=ncl + 1)[1:] / size
        ccol = np.bincount(lv, weights=jj, minlength=ncl + 1)[1:] / size
    else:
        crow, ccol = prow, pcol

    hist, npts = _pc_counts(prow, pcol, crow, ccol, nb, dR, rows, cols)
    g = np.zeros(nb)
    ok = (volume != 0) & (npts != 0)
    g[ok] = ((hist[ok] / volume[ok]) / npts[ok]) / intensity
    return g, r


# --------------------------------------------------------------------------
# ddf.c
# --------------------------------------------------------------------------
@jit
def _ddf_kernel(x, orow, ocol, ocost, inner, outer, nang, dA, pi_c):
    rows, cols = x.shape
    nr = inner.shape[0]
    T = np.zeros((nr, nang))
    area_bins = np.zeros((nr, nang))
    area_rings = np.zeros(nr)
    total = orow.shape[0]
    for m in range(1, total):
        cost = ocost[m]
        angle = math.atan2(-1.0 * orow[m], 1.0 * ocol[m])
        if angle < -dA / 2.0:
            angle = angle + 2 * pi_c
        angle += dA / 2.0
        for n in range(nr):
            if inner[n] <= cost and cost < outer[n]:
                area_rings[n] += 1
                q = int(math.floor(angle / dA))
                if q < 0 or q >= nang:
                    raise ValueError("Error in angular bin")
                area_bins[n, q] += 1
    for r in range(nr):
        for n in range(cols):
            for m in range(rows):
                if x[m, n] != 0:
                    if (m <= outer[r] or n <= outer[r]
                            or (rows - m) <= outer[r] or (cols - n) <= outer[r]):
                        continue
                    j = 1
                    while j < total and ocost[j] < outer[r]:
                        if ocost[j] >= inner[r]:
                            md = m + orow[j]
                            nd = n + ocol[j]
                            if x[md, nd] != 0:
                                angle = math.atan2(-1.0 * orow[j], 1.0 * ocol[j])
                                if angle < -dA / 2.0:
                                    angle = angle + 2 * pi_c
                                angle += dA / 2.0
                                q = int(math.floor(angle / dA))
                                if q < 0 or q >= nang:
                                    raise ValueError("Error in angular bin")
                                T[r, q] += 1
                        j += 1
        ring = 0.0
        for n in range(nang):
            ring += T[r, n]
        if ring > 0.0:
            for n in range(nang):
                if area_bins[r, n] > 0:
                    T[r, n] = (T[r, n] / area_bins[r, n]) / (ring / area_rings[r])
                else:
                    T[r, n] = 0.0
    return T


def ddf(X, ri, ro, N=72):
    """Directional distribution function (``ddf.c``).

    ``T, a = ddf(X, ri, ro, N=72)``: for each ring ri[k] <= |x-y| < ro[k],
    the normalised distribution of neighbour directions in N angular bins
    (T has shape (len(ri), N)); ``a`` holds the bin centre angles (radians).
    """
    x = as_binary(X)
    if x.ndim != 2:
        raise ValueError("ddf: X must be 2-D")
    inner = np.atleast_1d(np.asarray(ri, dtype=np.float64)).ravel()
    outer = np.atleast_1d(np.asarray(ro, dtype=np.float64)).ravel()
    if inner.size != outer.size:
        raise ValueError("Vectors Ri and Ro must be same length!")
    nang = int(N)
    dA = 2.0 * PI_C / nang
    rows, cols = x.shape
    orow, ocol, ocost = _cost_array(rows, cols)
    T = _ddf_kernel(x, orow, ocol, ocost, inner, outer, nang, dA, PI_C)
    a = np.arange(nang) * dA
    return T, a
