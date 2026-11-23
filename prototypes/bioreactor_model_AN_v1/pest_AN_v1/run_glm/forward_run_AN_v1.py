"""PEST-compatible forward run (AN_v1).

PEST/PEST++ writes the parameter file `params_AN_v1.dat` from the template `params_AN_v1.tpl`,
calls this script, and reads `output_AN_v1.dat` through the instruction file `output_AN_v1.ins`.
The script runs the SYNTHETIC-case model; the mass balance is checked inside (the run aborts
with a non-zero exit code if the water, N or P balance does not close -> PEST sees a failed run).
"""
import sys, pathlib
import os
HERE = pathlib.Path(__file__).resolve().parent
# locate the project root (folder containing bioreactor_AN_v1/): env var, or walk up from here / cwd
for base in [os.environ.get("BIOREACTOR_AN_V1_ROOT", ""), *HERE.parents, *pathlib.Path.cwd().resolve().parents]:
    if base and (pathlib.Path(base) / "bioreactor_AN_v1").is_dir():
        sys.path.insert(0, str(base)); break
import numpy as np
from bioreactor_AN_v1 import synthetic_AN_v1 as S
from bioreactor_AN_v1 import calib_common_AN_v1 as C

vals = {}
for line in open("params_AN_v1.dat"):
    k, v = line.split(); vals[k] = float(v)
x = np.array([vals[n.lower()] for n in S.NAMES])
out = S.run_model(x, C.FORCING, tol=1e-6)          # raises MassBalanceError if residual > tol
with open("output_AN_v1.dat", "w") as f:
    for i, d in enumerate(C.D):
        f.write(f"on_{d:03d} {out['cN'][d]:.10e}\n")
    for i, d in enumerate(C.D):
        f.write(f"op_{d:03d} {out['cP'][d]:.10e}\n")
    f.write(f"fkgn {out['kgN']:.10e}\n")
    f.write(f"fkgp {out['kgP']:.10e}\n")
