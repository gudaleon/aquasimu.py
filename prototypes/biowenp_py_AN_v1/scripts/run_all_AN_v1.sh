#!/usr/bin/env bash
# Part of biowenp_py_AN_v1 (LGPL-3.0-or-later). Runs tests, equivalence (needs octave), cross-check and demo.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -m pytest -q tests
TAG=${1:-run_$(date +%Y%m%d_%H%M%S)}     # new tag each run, so earlier outputs are never overwritten
python3 scripts/run_equivalence_octave_AN_v1.py 90 "$TAG"
python3 scripts/run_crosscheck_soilp_AN_v1.py
python3 scripts/run_demo_AN_v1.py
