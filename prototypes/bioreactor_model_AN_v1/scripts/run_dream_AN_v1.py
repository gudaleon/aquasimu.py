"""DREAM(ZS) calibration on the SYNTHETIC data set (AN_v1).  Own numpy implementation."""
import sys, time, pathlib
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
from bioreactor_AN_v1 import synthetic_AN_v1 as S
from bioreactor_AN_v1 import calib_common_AN_v1 as C
from bioreactor_AN_v1.dream_AN_v1 import dream_zs

OUT = HERE / "results_AN_v1"; OUT.mkdir(exist_ok=True)
NG1, NG2 = int(sys.argv[1]) if len(sys.argv) > 1 else 8000, int(sys.argv[2]) if len(sys.argv) > 2 else 14000

print("== DREAM(ZS) run 1: fixed (known) sigma, 5 parameters")
t = time.time()
r1 = dream_zs(C.logpost_fixed_sigma, S.Z_LO, S.Z_HI, n_chains=8, n_gen=NG1, seed=11, n_proc=8)
print(f"   {NG1*8} model runs in {time.time()-t:.0f} s, acceptance {r1['accept_rate']:.2f}")
np.savez_compressed(OUT / "dream_fixed_sigma_AN_v1.npz", **{k: v for k, v in r1.items()})

print("== DREAM(ZS) run 2: heteroscedastic sigma estimated, 9 parameters")
t = time.time()
lo = np.r_[S.Z_LO, C.ERR_LO]; hi = np.r_[S.Z_HI, C.ERR_HI]
r2 = dream_zs(C.logpost_hetero, lo, hi, n_chains=8, n_gen=NG2, seed=12, n_proc=8)
print(f"   {NG2*8} model runs in {time.time()-t:.0f} s, acceptance {r2['accept_rate']:.2f}")
np.savez_compressed(OUT / "dream_hetero_sigma_AN_v1.npz", **{k: v for k, v in r2.items()})
