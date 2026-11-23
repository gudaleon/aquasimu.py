"""PEST++ GLM then IES on the SYNTHETIC bioreactor case.

Requires pestpp-glm and pestpp-ies on PATH, or PESTPP_BIN pointing at their
directory. The binaries are not in this repository.

    python scripts/run_pest_synthetic.py

Products go to outputs/pest/run_glm and outputs/pest/run_ies. The observations
are synthetic.
"""
import shutil
import time
from pathlib import Path

from aquasimu.calib.pestpp import build_case, require_pest, run_glm, run_ies

require_pest("pestpp-glm", "pestpp-ies")
root = Path("outputs/pest")
case_dir = root / "case"
print("SYNTHETIC case. Writing the PEST control file.")
build_case(case_dir, noptmax=12, ies_num_reals=200)

glm_dir = root / "run_glm"
shutil.rmtree(glm_dir, ignore_errors=True)
shutil.copytree(case_dir, glm_dir)
t = time.time()
run_glm(glm_dir)
print(f"pestpp-glm finished in {time.time() - t:.0f} s")

template = root / "ies_template"
shutil.rmtree(template, ignore_errors=True)
shutil.copytree(case_dir, template)
ies_dir = root / "run_ies"
shutil.rmtree(ies_dir, ignore_errors=True)
t = time.time()
run_ies(ies_dir, template, noptmax=4, num_reals=200, num_workers=8, port=4004)
print(f"pestpp-ies finished in {time.time() - t:.0f} s")
shutil.rmtree(template, ignore_errors=True)
