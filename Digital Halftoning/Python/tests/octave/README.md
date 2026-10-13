# Regenerating the MATLAB/C reference data

`tests/test_parity.py` compares the Python ports with outputs of the MATLAB
and C (MEX) code of this repository, stored in `../data/octave_reference.npz`.
They were produced with GNU Octave 8.4 (packages `image` and `signal`).
To regenerate them (Linux, from this folder):

```bash
sudo apt install octave octave-dev octave-image octave-signal
REPO="../../.."                                   # the "Digital Halftoning" folder
mkdir -p mex && cp "$REPO/Point Process Statistics/"*.c mex/
(cd mex && for f in *.c; do mkoctfile --mex "$f"; done)
python3 make_inputs.py
octave-cli ref.m "$REPO/Matlab" mex "$REPO/Point Process Statistics" ref_fixed.mat
python3 build_reference.py
```

With MATLAB instead of Octave: compile the C files with `mex file.c`, then run
`ref.m` after setting `argv` to return the same four arguments.
