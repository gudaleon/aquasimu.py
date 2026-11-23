#!/usr/bin/env bash
# Reproduce everything (AN_v1).  Needs: python>=3.10 with numpy scipy pandas matplotlib numba pyemu pytest,
# and PEST++ binaries (pestpp-glm, pestpp-ies) -> set PESTPP_BIN to their folder.
set -e
cd "$(dirname "$0")/.."
export PATH="$(dirname "$(which python)"):$PATH"
python -m pytest -q tests                              # unit tests (mass balance, analytical solutions, DREAM)
python scripts/run_aquasim_examples_AN_v1.py           # download Eawag AQUASIM source+tutorial, build batch version, run examples
python scripts/run_benchmarks_AN_v1.py                 # Python vs AQUASIM consistency check -> tables/, figures/
python scripts/run_dream_AN_v1.py 8000 14000           # DREAM(ZS) on the synthetic case
python scripts/run_pest_AN_v1.py                       # PEST++ GLM + IES (pyEMU)
python scripts/compare_AN_v1.py                        # recovery, R-hat, CIs, PEST vs DREAM, mass-balance table
