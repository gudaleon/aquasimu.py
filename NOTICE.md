# NOTICE

This repository does not yet have a top-level licence. The packaged aquasimu
library (see pull request #1) and the prototypes in this tree, other than the
exception below, should not be treated as released under a specific licence
until one is chosen.

## BioWEN-P Python port (LGPL-3.0-or-later)

`prototypes/biowenp_py_AN_v1/` is a derivative of BioWEN-P (Zhou, Shafii,
Parsons, Passeport, Rezanezhad, Lisogorsky, Van Cappellen, 2023,
*Sci. Total Environ.* 876:162749, doi:10.1016/j.scitotenv.2023.162749;
MATLAB code FRDR doi:10.20383/103.0640) and is therefore distributed under
**LGPL-3.0-or-later**. See that subdirectory's `LICENSE`, `COPYING.LESSER`,
`COPYING`, and `NOTICE_AN_v1.md`. The original MATLAB `.m` files are not
redistributed here.

That LGPL applies only to that subdirectory (and copies of it). It does not
licence the rest of aquasimu.py.

## What is not in this repository

Eawag AQUASIM source, binaries, `.aqu` example files, and AQUASIM output
listings are not included (licence status unclear). The FRDR BioWEN-P MATLAB
snapshot is not included. Large regenerable Octave/Python result binaries
(`*.mat`, `*.npz` from the 90-day equivalence run) are gitignored; see
`prototypes/biowenp_py_AN_v1/results_AN_v1/README.md`.
