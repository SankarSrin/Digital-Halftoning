#Digital Halftoning
#Mr.Sankarasrinivasan S
#Ref:  Digital Halftone Database (DHD): A Comprehensive Analysis on Halftone Types
"""Demo: halftone an image with every method and compare their statistics.

    python main.py                      # uses peppers.tif, writes ./output
    python main.py photo.png -o out --no-dbs

Writes one PNG per halftone and, if matplotlib is installed, a figure with
the radially averaged power spectrum (RAPSD) and pair correlation of each.
"""

import argparse
import os
import time

import numpy as np
from PIL import Image

import halftoning as h
from halftoning import point_process as pps

METHODS = {
    "bayer_dispersed": lambda im: h.bayers(im, 1),
    "bayer_clustered": lambda im: h.bayers(im, 2),
    "ulichney_clustered_4x4": lambda im: h.ulichney(im, 1),
    "ulichney_clustered_8x8": lambda im: h.ulichney(im, 2),
    "ulichney_void_cluster": lambda im: h.ulichney(im, 3),
    "dotdiff_mese_8x8": lambda im: h.dotdiff(im, 1),
    "dotdiff_knuth_8x8": lambda im: h.dotdiff(im, 2),
    "dotdiff_mese_16x16": lambda im: h.dotdiff(im, 3),
    "floyd_steinberg_serpentine": lambda im: h.floyd_serp(im)[0],
    "edbs": lambda im: h.edbs(im, seed=0),
}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("image", nargs="?", default=os.path.join(here, "peppers.tif"))
    ap.add_argument("-o", "--out", default="output")
    ap.add_argument("--no-dbs", action="store_true", help="skip direct binary search (slowest)")
    args = ap.parse_args()

    img = np.array(Image.open(args.image))
    gray = h.to_gray_double(img)
    os.makedirs(args.out, exist_ok=True)
    print(f"input {args.image}: {gray.shape[0]}x{gray.shape[1]}, mean gray {gray.mean():.4f}")

    results = {}
    for name, fn in METHODS.items():
        if name == "edbs" and args.no_dbs:
            continue
        t = time.time()
        out = np.asarray(fn(img), dtype=bool)
        results[name] = out
        Image.fromarray((out * 255).astype(np.uint8)).save(os.path.join(args.out, f"{name}.png"))
        print(f"{name:28s} {time.time() - t:6.2f} s   mean {out.mean():.4f}")

    # Point-process statistics on a 256x256 centre crop
    r0 = max(0, gray.shape[0] // 2 - 128)
    c0 = max(0, gray.shape[1] // 2 - 128)
    stats = {}
    for name, out in results.items():
        X = pps.pp(out[r0:r0 + 256, c0:c0 + 256])
        P, fr, _ = pps.rapsd(X.astype(float), 64)
        g, r = pps.pc(X, 1)
        stats[name] = (fr, P, r[:-1], g[:-1])  # pc's last bin is always empty
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed: skipping statistics figure")
        return
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for name, (fr, P, r, g) in stats.items():
        ax[0].plot(fr, P, label=name, lw=1)
        ax[1].plot(r, g, label=name, lw=1)
    ax[0].set(xlabel="radial frequency (cycles/pixel)", ylabel="RAPSD / g(1-g)", title="Radially averaged power spectrum")
    ax[1].set(xlabel="radius (pixels)", ylabel="pair correlation", title="Pair correlation of minority pixels", xlim=(0, 30))
    ax[0].set_yscale("log")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "statistics.png"), dpi=120)
    print(f"wrote {args.out}/statistics.png")


if __name__ == "__main__":
    main()
