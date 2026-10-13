"""Behavioural tests that do not depend on MATLAB/Octave reference data."""

import numpy as np
import pytest

import halftoning as h
from halftoning import point_process as P
from halftoning.dot_diffusion import CM_KNUTH_8, CM_MESE_8, CM_MESE_16
from halftoning.ordered import (BAYER_CLUSTERED, BAYER_DISPERSED, ULICHNEY_CLUSTERED_4,
                                ULICHNEY_VOID_CLUSTER_14)

HALFTONERS = {
    "bayer1": lambda im: h.bayers(im, 1),
    "bayer2": lambda im: h.bayers(im, 2),
    "ulichney1": lambda im: h.ulichney(im, 1),
    "ulichney2": lambda im: h.ulichney(im, 2),
    "ulichney3": lambda im: h.ulichney(im, 3),
    "dotdiff1": lambda im: h.dotdiff(im, 1),
    "dotdiff2": lambda im: h.dotdiff(im, 2),
    "dotdiff3": lambda im: h.dotdiff(im, 3),
    "floyd": lambda im: h.floyd_serp(im)[0],
    "edbs": lambda im: h.edbs(im, seed=1),
}


@pytest.mark.parametrize("name", HALFTONERS)
def test_shape_binary_and_extremes(name):
    f = HALFTONERS[name]
    for shape in [(37, 53), (16, 16)]:
        out = np.asarray(f(np.full(shape, 0.37)))
        assert out.shape == shape
        assert set(np.unique(out)) <= {0, 1}
    assert np.all(np.asarray(f(np.zeros((32, 40)))) == 0)
    assert np.all(np.asarray(f(np.ones((32, 40)))) == 1)


@pytest.mark.parametrize("name", HALFTONERS)
@pytest.mark.parametrize("gray", [0.25, 0.5, 0.75])
def test_mean_tone(name, gray):
    size = 48 if name == "edbs" else 112  # multiple of 4, 8, 14 and 16
    out = np.asarray(HALFTONERS[name](np.full((size, size), gray)), dtype=float)
    assert abs(out.mean() - gray) < 0.02, out.mean()


def test_uint8_rgb_inputs():
    rng = np.random.default_rng(0)
    rgb = rng.integers(0, 256, (24, 30, 3), dtype=np.uint8)
    gray = h.to_gray_double(rgb)
    assert gray.shape == (24, 30)
    np.testing.assert_array_equal(h.bayers(rgb), h.bayers(gray))
    np.testing.assert_array_equal(h.dotdiff(rgb, 2), h.dotdiff(gray, 2))
    np.testing.assert_allclose(h.im2double(np.array([0, 255], np.uint8)), [0.0, 1.0])


@pytest.mark.parametrize("m", [CM_MESE_8, CM_KNUTH_8, CM_MESE_16, ULICHNEY_VOID_CLUSTER_14 * 196,
                               BAYER_DISPERSED * 64, BAYER_CLUSTERED * 64, ULICHNEY_CLUSTERED_4 * 16])
def test_class_and_threshold_matrices_are_permutations(m):
    m = np.rint(m).astype(int)
    assert sorted(m.ravel()) == list(range(m.size))


def _naive_errdiff(img, weights, serpentine):
    """Independent error diffusion: weights = {(dy, dx): w} relative to the
    current pixel; the x offsets are mirrored on right-to-left rows."""
    a = img.astype(float).copy()
    rows, cols = a.shape
    out = np.zeros_like(a)
    for y in range(rows):
        rtl = serpentine and y % 2 == 1
        for x in (range(cols - 1, -1, -1) if rtl else range(cols)):
            o = 1.0 if a[y, x] >= 0.5 else 0.0
            e = o - a[y, x]
            out[y, x] = o
            for (dy, dx), w in weights.items():
                yy, xx = y + dy, x + (-dx if rtl else dx)
                if 0 <= yy < rows and 0 <= xx < cols:
                    a[yy, xx] -= e * w
    return out


@pytest.mark.parametrize("serp", [False, True])
def test_errdiff_against_independent_implementation(serp):
    rng = np.random.default_rng(3)
    img = rng.random((23, 31))
    # asymmetric filter whose origin is not the middle column
    fc = np.array([[0, 0, -99 * 16, 7], [1, 3, 5, 0]]) / 16
    weights = {(0, 1): 7 / 16, (1, -2): 1 / 16, (1, -1): 3 / 16, (1, 0): 5 / 16}
    out = h.errdiff(img, fc, 0, -1 if serp else 1)[0]
    np.testing.assert_array_equal(out, _naive_errdiff(img, weights, serp))


def test_dbs_converges_and_lowers_error():
    rng = np.random.default_rng(0)
    img = np.clip(np.linspace(0, 1, 40)[None, :] + 0.0 * rng.random((32, 40)), 0, 1)
    out, it = h.edbs(img, seed=0, return_iterations=True)
    # a further run from the result must not change anything
    again = h.edbs(img, init=out)
    np.testing.assert_array_equal(out, again)
    assert it >= 2


def test_pp_and_intensity():
    X = np.zeros((4, 4))
    X[:3] = 1
    Y = P.pp(X)
    assert Y.dtype == bool and Y.sum() == 4
    np.testing.assert_allclose(P.intnsty(Y), [0.25])
    with pytest.raises(ValueError):
        P.pp(np.full((2, 2), 0.5))


def test_pair_correlation_of_random_pattern_is_flat():
    rng = np.random.default_rng(5)
    X = rng.random((128, 128)) < 0.2
    g, r = P.pc(X, 2)
    assert np.all(np.abs(g[1:-1] - 1) < 0.1)
    assert g[-1] == 0  # pc.c quirk: the last radial bin is never filled (see REVIEW.md)
    g, r = P.vpc(X, (1, 1), 2)
    assert np.all(np.abs(g[1:] - 1) < 0.1)


def test_rapsd_of_white_noise_is_flat():
    rng = np.random.default_rng(6)
    X = (rng.random((256, 256)) < 0.3).astype(float)
    p, fr, n = P.rapsd(X, 64)
    assert n == 16 and fr.size == 45
    assert abs(p[2:].mean() - 1) < 0.05


def test_spectral_argument_checks():
    X = np.zeros((64, 64))
    X[::3, ::5] = 1
    with pytest.raises(ValueError):
        P.anis(X, 32, 16)  # spectral window smaller than spatial window
    with pytest.raises(ValueError):
        P.rapsd(X, 128)  # window larger than image
    a = P.rapsd(X, 32, None, [8, 16])[2]
    b = P.rapsd(X, 32, None, 8)[2]
    assert (a, b) == (15, 25)  # [DM DN] honoured: 5x3 vs 5x5 windows
