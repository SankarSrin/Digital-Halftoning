"""Pack inputs.mat + ref_fixed.mat into ../data/octave_reference.npz."""
import os
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
r = sio.loadmat(os.path.join(HERE, "ref_fixed.mat"))["R"][0, 0]
I = sio.loadmat(os.path.join(HERE, "inputs.mat"))
d = {"in_" + k: I[k] for k in I if not k.startswith("__")}
for k in r.dtype.names:
    v = r[k]
    if v.dtype.kind == "U":
        raise SystemExit(f"{k}: MATLAB/C raised an error: {v[0]}")
    d["ref_" + k] = np.asarray(v, dtype=float)
np.savez_compressed(os.path.join(HERE, "..", "data", "octave_reference.npz"), **d)
print(f"wrote {len(d)} arrays")
