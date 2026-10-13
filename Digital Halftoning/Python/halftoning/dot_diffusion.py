"""Dot diffusion (Knuth 1987). Python port of ``Matlab/dotdiff.m`` (fixed)."""

import numpy as np

from ._utils import jit, tile_to, to_gray_double

# Mese-Vaidyanathan 8x8 class matrix (met=1)
CM_MESE_8 = np.array([
    [47, 31, 51, 24, 27, 45, 5, 21],
    [37, 63, 53, 11, 22, 4, 1, 33],
    [61, 0, 57, 16, 26, 29, 46, 8],
    [20, 14, 9, 62, 18, 41, 38, 6],
    [17, 13, 25, 15, 55, 48, 52, 58],
    [3, 7, 2, 32, 30, 34, 56, 60],
    [28, 40, 36, 39, 49, 43, 35, 10],
    [54, 23, 50, 12, 42, 59, 44, 19],
])

# Knuth 8x8 class matrix (met=2); row 7, col 8 is 44 (the MATLAB had 42 twice)
CM_KNUTH_8 = np.array([
    [34, 48, 40, 32, 29, 15, 23, 31],
    [42, 58, 56, 53, 21, 5, 7, 10],
    [50, 62, 61, 45, 13, 1, 2, 18],
    [38, 46, 54, 37, 25, 17, 9, 26],
    [28, 14, 22, 30, 35, 49, 41, 33],
    [20, 4, 6, 11, 43, 59, 57, 52],
    [12, 0, 3, 19, 51, 63, 60, 44],
    [24, 16, 8, 27, 39, 47, 55, 36],
])

# Mese-Vaidyanathan 16x16 class matrix (met=3)
CM_MESE_16 = np.array([
    [207, 0, 13, 17, 28, 55, 18, 102, 81, 97, 74, 144, 149, 169, 170, 172],
    [3, 6, 23, 36, 56, 50, 65, 87, 145, 130, 137, 158, 182, 184, 195, 221],
    [7, 14, 24, 37, 67, 69, 86, 5, 106, 152, 150, 165, 183, 192, 224, 1],
    [15, 26, 43, 53, 51, 101, 115, 131, 139, 136, 166, 119, 208, 223, 226, 4],
    [22, 39, 52, 71, 84, 103, 164, 135, 157, 173, 113, 190, 222, 225, 227, 16],
    [40, 85, 72, 83, 104, 117, 167, 133, 168, 180, 200, 219, 231, 228, 12, 21],
    [47, 120, 54, 105, 123, 132, 146, 176, 179, 202, 220, 230, 245, 2, 20, 41],
    [76, 73, 127, 109, 138, 134, 178, 181, 206, 196, 229, 244, 246, 19, 42, 49],
    [80, 99, 112, 147, 142, 171, 177, 203, 218, 232, 243, 248, 247, 33, 48, 68],
    [108, 107, 140, 143, 185, 163, 204, 217, 233, 242, 249, 255, 44, 45, 70, 79],
    [110, 141, 88, 75, 175, 205, 214, 234, 241, 250, 254, 38, 46, 77, 116, 100],
    [111, 148, 160, 174, 201, 215, 235, 240, 251, 252, 253, 61, 62, 93, 94, 125],
    [151, 159, 189, 199, 197, 216, 236, 239, 25, 31, 60, 82, 92, 95, 124, 114],
    [156, 188, 191, 209, 213, 237, 238, 29, 32, 59, 64, 91, 118, 78, 128, 155],
    [187, 194, 198, 212, 9, 10, 30, 35, 58, 63, 90, 96, 122, 129, 154, 161],
    [193, 210, 211, 8, 11, 27, 34, 57, 66, 89, 98, 121, 126, 153, 162, 186],
])

# Diffusion weights (normalised per pixel by the weights that are still open)
DM = np.array([[1, 2, 1], [2, 0, 2], [1, 2, 1]], dtype=np.float64)


@jit
def _dotdiff_kernel(im, cm, rows, cols, dm):
    H = np.zeros(im.shape)
    for t in range(rows.shape[0]):
        r = rows[t]
        c = cols[t]
        ii = cm[r, c]
        if im[r, c] > 0.5:
            H[r, c] = 1.0
        qerr = im[r, c] - H[r, c]
        sw = 0.0
        for a in range(3):
            for b in range(3):
                if cm[r - 1 + a, c - 1 + b] > ii:
                    sw += dm[a, b]
        if sw > 0:
            f = qerr / sw
            for a in range(3):
                for b in range(3):
                    w = dm[a, b] if cm[r - 1 + a, c - 1 + b] > ii else 0.0
                    im[r - 1 + a, c - 1 + b] = im[r - 1 + a, c - 1 + b] + f * w
    return H


def dotdiff(im, met=1, class_matrix=None):
    """Dot diffusion halftoning (``dotdiff.m``).

    met=1: Mese-Vaidyanathan 8x8, met=2: Knuth 8x8, otherwise Mese-Vaidyanathan
    16x16. A custom ``class_matrix`` (permutation of 0..n-1) overrides ``met``.
    Returns a float 0/1 halftone the same size as ``im``.
    """
    g = to_gray_double(im)
    if class_matrix is None:
        class_matrix = CM_MESE_8 if met == 1 else CM_KNUTH_8 if met == 2 else CM_MESE_16
    s1, s2 = g.shape
    pad = np.zeros((s1 + 2, s2 + 2))
    pad[1:-1, 1:-1] = g
    cm = np.full((s1 + 2, s2 + 2), -1, dtype=np.int64)  # border: never receives error
    cm[1:-1, 1:-1] = tile_to(class_matrix, (s1, s2))
    # Visit order = MATLAB's: class by class, column-major within a class.
    rr, cc = np.nonzero(cm >= 0)
    order = np.lexsort((rr, cc, cm[rr, cc]))
    H = _dotdiff_kernel(pad, cm, rr[order].astype(np.int64), cc[order].astype(np.int64), DM)
    return H[1:-1, 1:-1]
