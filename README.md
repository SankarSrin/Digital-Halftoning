# Digital-Halftoning
The repository contains MATLAB and Python code for the digital halftoning techniques

Code provided for

**MATLAB** (`Digital Halftoning/Matlab`):
- Bayers (dispersed and clustered dot)
- Ulichney (4x4 and 8x8 clustered dot, 14x14 void-and-cluster)
- Dot-Diffusion (dotdiff.m : can perform Knuth and Mese-Vaidhyanathan dot-diffusion)
- Error-Diffusion (floyd_serp.m  : can perform floyd-steinberg error diffusion, errdiff.m : will be called during exectution)
- Direct binary search (DBS) (EDBS.m : efficient DBS)

**Python** (`Digital Halftoning/Python`): every MATLAB function above plus the
Point Process Statistics, as the `halftoning` package:

```python
import numpy as np
from PIL import Image
import halftoning as h
from halftoning import point_process as pps

img = np.array(Image.open("peppers.tif"))      # uint8 or float 0..1, gray or RGB, any size
b  = h.bayers(img, 1)                          # Bayer dispersed (2 = clustered)
u  = h.ulichney(img, 3)                        # Ulichney 1/2 = clustered, 3 = void-and-cluster
d  = h.dotdiff(img, 2)                         # 1 = Mese 8x8, 2 = Knuth, 3 = Mese 16x16
fs, qerr, k = h.floyd_serp(img)                # Floyd-Steinberg, serpentine
o, qerr, k = h.errdiff(h.im2double(img), h.JARVIS, 0, -1)
dbs = h.edbs(img, seed=0)                      # efficient direct binary search

X = pps.pp(fs)                                 # minority-pixel point process
P, fr, n = pps.rapsd(X.astype(float), 64)      # radially averaged power spectrum
g, r = pps.pc(X, 1)                            # pair correlation (also vpc, rsmm, ddf, anis)
```

```bash
cd "Digital Halftoning/Python"
pip install -r requirements.txt
python main.py            # halftones peppers.tif with every method -> ./output
python -m pytest tests    # 106 tests, incl. output-for-output parity with the MATLAB/C code
```

`numba` is optional but makes error diffusion, dot diffusion, DBS and the
point-process statistics fast.

**Point Process Statistics** (`Digital Halftoning/Point Process Statistics`):
the code for Point Process Statistics is also included.
The source code is taken from Modern Digital Halftoning: refer to the software guide for details.
Recompile the C files with `mex file.c` (the bundled `.mexw64` binaries predate the fixes in REVIEW.md).

**Code review:** see [REVIEW.md](REVIEW.md) for the bugs found and fixed in the
MATLAB and C code (2026) and how the Python port was verified.


For details refer to blog:

https://imageprocessing-sankarsrin.blogspot.com/2016/12/digital-half-toning-ordered-dithering.html



Ref Paper: 

**__Guo, Jing Ming, and S. Sankarasrinivasan. "Digital Halftone Database (DHD): A Comprehensive Analysis on Halftone Types." 2018 Asia-Pacific Signal and Information Processing Association Annual Summit and Conference (APSIPA ASC). IEEE, 2018.__**
