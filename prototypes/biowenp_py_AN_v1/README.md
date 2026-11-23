# biowenp.py

Python port of **BioWEN-P** (Zhou, Shafii, Parsons, Passeport, Rezanezhad, Lisogorsky, Van Cappellen, 2023)
plus a saturated flow-through augmentation for eggshell + woodchip bioreactors.

- Paper: [doi:10.1016/j.scitotenv.2023.162749](https://doi.org/10.1016/j.scitotenv.2023.162749)
- Original MATLAB code (not redistributed here): [FRDR doi:10.20383/103.0640](https://doi.org/10.20383/103.0640)
- Licence: **LGPL-3.0-or-later** (`COPYING.LESSER`, `COPYING`; credits in `NOTICE_AN_v1.md`)

The full scientific write-up, equivalence numbers, quirks reproduced from the original MATLAB, demo results
and gaps are in **[README_AN_v1.md](README_AN_v1.md)**. Every comparison in this repository is a
*consistency check* (code vs code, or code vs the published tables). Nothing here is validated against field data.
All committed runs use **synthetic** forcing.

## Quick start

```bash
pip install -r requirements_AN_v1.txt
python3 -m pytest -q tests
python3 scripts/run_demo_AN_v1.py
```

Optional (needs GNU Octave and a local copy of the FRDR model code):

```bash
python3 scripts/run_equivalence_octave_AN_v1.py 90 new_tag
```

Large regenerable binaries (`results_AN_v1/*.mat`, `*.npz`, about 28 MB and 50 MB) are gitignored.
How to regenerate them is in [`results_AN_v1/README.md`](results_AN_v1/README.md).

The original BioWEN-P `.m` files and the FRDR snapshot are **not** in this repository. The Octave driver
(`scripts/octave_AN_v1/`) is ours: it calls the unmodified original from a local snapshot path.
