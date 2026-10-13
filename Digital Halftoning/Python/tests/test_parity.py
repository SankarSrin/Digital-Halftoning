"""Python vs MATLAB/C parity.

``tests/data/octave_reference.npz`` holds test inputs and the outputs of the
(fixed) MATLAB functions and C MEX functions of this repository, computed with
GNU Octave 8.4 (see ``tests/octave/README.md`` to regenerate). Every Python
port must reproduce them: binary halftones exactly, statistics to within
floating-point rounding (FFT / summation order).
"""

import os

import numpy as np
import pytest

import halftoning as h
from halftoning import point_process as P

DATA = np.load(os.path.join(os.path.dirname(__file__), "data", "octave_reference.npz"))
IN = {k[3:]: DATA[k] for k in DATA.files if k.startswith("in_")}
REF = {k[4:]: DATA[k] for k in DATA.files if k.startswith("ref_")}
X, Xc, Xr, X3 = (IN[k].astype(bool) for k in ("X", "Xc", "Xr", "X3"))
ASYM = np.array([[0, 0, -99 * 16, 7], [1, 3, 5, 0]]) / 16


def check(name, value, exact=False):
    ref = np.asarray(REF[name], dtype=float)
    val = np.asarray(value, dtype=float)
    if val.shape != ref.shape:
        ref = ref.reshape(val.shape)
    if exact:
        np.testing.assert_array_equal(val, ref, err_msg=name)
    else:
        np.testing.assert_allclose(val, ref, rtol=1e-10, atol=1e-12, equal_nan=True, err_msg=name)


IMAGES = ["pep", "crop", "ramp", "rgb"]


@pytest.mark.parametrize("img", IMAGES)
@pytest.mark.parametrize("met", [1, 2])
def test_bayers(img, met):
    check(f"bayers{met}_{img}", h.bayers(IN[img], met), exact=True)


@pytest.mark.parametrize("img", IMAGES)
@pytest.mark.parametrize("met", [1, 2, 3])
def test_ulichney(img, met):
    check(f"ulichney{met}_{img}", h.ulichney(IN[img], met), exact=True)


@pytest.mark.parametrize("img", IMAGES)
@pytest.mark.parametrize("met", [1, 2, 3])
def test_dotdiff(img, met):
    check(f"dotdiff{met}_{img}", h.dotdiff(IN[img], met), exact=True)


@pytest.mark.parametrize("img", IMAGES)
def test_floyd_serp(img):
    check(f"floyd_{img}", h.floyd_serp(IN[img])[0], exact=True)


def test_errdiff_variants():
    c = h.im2double(IN["crop"])
    o, q, k = h.errdiff(c, h.JARVIS, 0, 1)
    check("ed_jarvis_raster", o, exact=True)
    check("ed_jarvis_raster_q", q, exact=True)
    check("ed_jarvis_raster_k", k)
    o, q, k = h.errdiff(c, h.JARVIS, 0, -1)
    check("ed_jarvis_serp", o, exact=True)
    check("ed_jarvis_serp_k", k)
    check("ed_asym_serp", h.errdiff(c, ASYM, 0, -1)[0], exact=True)
    check("ed_asym_raster", h.errdiff(c, ASYM, 0, 1)[0], exact=True)
    o, q, k = h.errdiff(c, h.FLOYD_STEINBERG, 0.5, -1)
    check("ed_fs_l05", o, exact=True)
    check("ed_fs_l05_k", k)


@pytest.mark.parametrize("which", ["a", "b"])
def test_edbs(which):
    out = h.edbs(IN[f"dbs_{which}"], init=IN[f"init_{which}"])
    check(f"edbs_{which}", out, exact=True)


def test_pp_intnsty_bayer():
    check("pp_X3", P.pp(X3), exact=True)
    check("pp_inv", P.pp(~X), exact=True)
    check("int_X3", P.intnsty(X3))
    check("lau_bayer", P.bayer(IN["dimg"]), exact=True)


SPECTRAL = [(32,), (64, 128, [16, 32]), (33,), (16, 20), (64, 64, 8), (48, 64)]


@pytest.mark.parametrize("k,args", list(enumerate(SPECTRAL, 1)))
def test_rapsd_anis(k, args):
    D = X.astype(float)
    p, f, cnt = P.rapsd(D, *args)
    check(f"rapsd{k}", p)
    check(f"rapsd{k}_f", f)
    check(f"rapsd{k}_c", cnt)
    a, _, cnt = P.anis(D, *args)
    check(f"anis{k}", a)
    check(f"anis{k}_c", cnt)


def test_rapsd_rect():
    check("rapsd_rect", P.rapsd(Xr.astype(float), 30)[0])


def test_pc():
    g, r = P.pc(X, 1)
    check("pc1", g)
    check("pc1_r", r)
    check("pc2", P.pc(X, 0.5, 30)[0])
    check("pc3", P.pc(Xc, 1, 40, "clus")[0])
    check("pc4", P.pc(Xr, 0.7)[0])
    check("pc5", P.pc(Xr.T, 1, 100)[0])
    check("pc6", P.pc(X, 1, 32, "xyz")[0])
    check("pc7", P.pc(Xc, 2, 60, "cluster")[0])


def test_vpc():
    g, r = P.vpc(X, (1, 1), 0.5)
    check("vpc1", g)
    check("vpc1_r", r)
    check("vpc2", P.vpc(X3, (1, 2), 1)[0])
    check("vpc3", P.vpc(X3, (3, 1), 1, 20)[0])
    check("vpc4", P.vpc(Xr, (1, 1), 0.75)[0])


def test_rsmm():
    check("rsmm1", P.rsmm(X))
    check("rsmm2", P.rsmm(Xc, 20))
    check("rsmm3", P.rsmm(Xr))


def test_ddf():
    T, a = P.ddf(X, [1, 3, 5], [3, 5, 8])
    check("ddf1", T)
    check("ddf1_a", a)
    check("ddf2", P.ddf(Xc, [2, 6], [6, 12], 36)[0])
    check("ddf3", P.ddf(Xr, [0.5], [4], 12)[0])
