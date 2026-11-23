"""DREAM(ZS) on the SYNTHETIC bioreactor case.

The observations, forcing and "true" parameters are seeded random numbers, not
field data. A full run is thousands of reactor simulations per chain:

    python scripts/run_dream_synthetic.py 8000 14000

The first argument is the number of generations for the fixed-sigma run
(5 parameters). The second is the heteroscedastic run (9 parameters).
Eight chains are used. Chains are written to outputs/dream/.
"""
import sys
import time
from pathlib import Path

import numpy as np

from aquasimu.calib.dream import dream_zs
from aquasimu.examples.bioreactor import likelihood as case_likelihood
from aquasimu.examples.bioreactor import synthetic as case

OUT = Path("outputs/dream")
OUT.mkdir(parents=True, exist_ok=True)
N_FIXED = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
N_HETERO = int(sys.argv[2]) if len(sys.argv) > 2 else 14000

print("SYNTHETIC case. DREAM(ZS) run 1: fixed (known) sigma, 5 parameters")
t = time.time()
fixed = dream_zs(
    case_likelihood.logpost_fixed_sigma, case.Z_LO, case.Z_HI,
    n_chains=8, n_gen=N_FIXED, seed=11, n_proc=8,
)
print(f"   {N_FIXED * 8} model runs in {time.time() - t:.0f} s, acceptance {fixed['accept_rate']:.2f}")
np.savez_compressed(OUT / "dream_fixed_sigma.npz", **fixed)

print("SYNTHETIC case. DREAM(ZS) run 2: heteroscedastic sigma estimated, 9 parameters")
t = time.time()
lo = np.r_[case.Z_LO, case_likelihood.ERR_LO]
hi = np.r_[case.Z_HI, case_likelihood.ERR_HI]
hetero = dream_zs(
    case_likelihood.logpost_hetero, lo, hi,
    n_chains=8, n_gen=N_HETERO, seed=12, n_proc=8,
)
print(f"   {N_HETERO * 8} model runs in {time.time() - t:.0f} s, acceptance {hetero['accept_rate']:.2f}")
np.savez_compressed(OUT / "dream_hetero_sigma.npz", **hetero)
