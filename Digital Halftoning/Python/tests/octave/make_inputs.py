"""Write inputs.mat (test images and point processes) for ref.m."""
import os, sys, numpy as np, scipy.io as sio
HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, PYDIR)
from PIL import Image
import halftoning as h, halftoning.point_process as P
rng = np.random.default_rng(7)
pep = np.array(Image.open(os.path.join(PYDIR, 'peppers.tif')))
crop = pep[30:91, 100:177].copy()                       # 61x77, odd size
ramp = np.tile(np.linspace(0, 1, 50), (40, 1))           # double 40x50
rgb = np.stack([pep[:48,:56], pep[100:148,40:96], pep[200:248,150:206]], 2).astype(np.uint8)
dbs_a = pep[100:148, 60:108].copy()                     # 48x48
dbs_b = crop[:32, :40].copy()                           # 32x40
init_a = (rng.random(dbs_a.shape) > 0.5).astype(float)
init_b = (rng.random(dbs_b.shape) > 0.5).astype(float)
X = P.pp(h.floyd_serp(pep)[0])                           # FS point process 256x256
Xc = P.pp(h.ulichney(pep, 1))                            # clustered
Xr = rng.random((64, 90)) < 0.15                         # random rect
X3 = np.stack([X[:128,:128], P.pp(h.bayers(pep))[:128,:128], Xr[:64,:64].repeat(2,0).repeat(2,1)], 2)
dimg = rng.random((37, 45))
sio.savemat(os.path.join(HERE, 'inputs.mat'), dict(pep=pep, crop=crop, ramp=ramp, rgb=rgb, dbs_a=dbs_a, dbs_b=dbs_b,
            init_a=init_a, init_b=init_b, X=X, Xc=Xc, Xr=Xr, X3=X3, dimg=dimg), do_compression=True)
print('ok', X.mean(), Xc.mean(), Xr.mean())
