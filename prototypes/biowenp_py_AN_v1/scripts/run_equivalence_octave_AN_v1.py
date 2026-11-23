#!/usr/bin/env python3
# Part of biowenp_py_AN_v1 (LGPL-3.0-or-later). MATLAB(Octave)-vs-Python consistency check on SYNTHETIC forcing.
"""Runs the ORIGINAL BioWEN-P MATLAB code (unmodified, from the read-only FRDR snapshot) in GNU Octave and the
Python port on the same SYNTHETIC 10-minute forcing with the Table 6 parameters, then compares every output array.
This is a code-equivalence (consistency) check, not a validation against field data."""
import json, os, subprocess, sys, time
import numpy as np, pandas as pd
from scipy.io import loadmat
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from biowenp_AN_v1 import HydroParams, PParams, HydroForcing, run_hydrology, run_p_faithful
from biowenp_AN_v1.massbalance_AN_v1 import water_balance_original, water_balance_full, p_balance_original
from biowenp_AN_v1.synthetic_AN_v1 import synthetic_stormwater_10min

SNAP = "/workspace/parsons_bioretention_AN_v1/frdr_8c090e43_model_code"
n_days = int(sys.argv[1]) if len(sys.argv) > 1 else 60
tag = sys.argv[2] if len(sys.argv) > 2 else "AN_v1"
res = os.path.join(ROOT, "results_AN_v1"); tab = os.path.join(ROOT, "tables_AN_v1")
h, p = HydroParams(), PParams()
f = synthetic_stormwater_10min(n_days=n_days)
fcsv = os.path.join(res, f"synthetic_forcing_10min_{n_days}d_{tag}.csv")
pd.DataFrame(dict(datenum=f["datenum"], precip_mm_per_10min=f["precip"], temp_C=f["temp"])).to_csv(fcsv, index=False, float_format="%.12g")
mat = dict(A_BC=h.A_BC, A_PP=h.A_PP, DA=h.DA, DA_PL=h.DA_PL, f=h.f, la=h.la, Tt=h.Tt, DD=h.DD, Void_max=h.Void_max,
           n_PET=h.n_PET, R_Pref=h.R_Pref, i_t_c=0, Sp_max=h.Sp_max, hfm_u=h.hfm_u, dfm=h.dfm, n_fm=h.n_fm, rho_fm=h.rho_fm,
           Ks=h.Ks, sfc=h.sfc, s_h=h.s_h, beta=h.beta, K_1=h.K_1, u=h.u, b_Sr=h.b_Sr, n_bs=h.n_bs,
           PETMethod=h.PETMethod, FlowEstimateMethod=h.FlowEstimateMethod, InfiltrationMethod=1,
           Inflow_TP_EstimateMethod=h.Inflow_TP_EstimateMethod,
           alpha_POP=p.alpha_POP, alpha_PIP=p.alpha_PIP, alpha_TDP=p.alpha_TDP, beta_Microbial=p.beta_Microbial,
           beta_Labile=p.beta_Labile, beta_Stable=p.beta_Stable, beta_SorbedP=p.beta_SorbedP, beta_MineralP=p.beta_MineralP,
           beta_RecalcitrantP=p.beta_RecalcitrantP, beta_SRP=p.beta_SRP, C_overflow=p.C_overflow, lambda_=p.lam,
           alpha_L=p.alpha_L, Kf=p.Kf, Pu_Depth=p.Pu_Depth, I_Pu_max=p.I_Pu_max, C_SRP_Pu_min=p.C_SRP_Pu_min, Km_Pu=p.Km_Pu,
           k_LF=p.k_LF, k_ad=p.k_ad, Qmax_ad=p.Qmax_ad, kL_ad=p.kL_ad, k_mi=p.k_mi, Um=p.Um, Km_Im=p.Km_Im,
           f_debris=p.f_debris, k_Death=p.k_Death, k_Hy=p.k_Hy, lambda_s=p.lambda_s, dx=1, Macro_percent=1,
           SorptionEquilibrium=2)
mat["lambda"] = mat.pop("lambda_")
fjson = os.path.join(res, f"params_matlab_names_{tag}.json")
json.dump(mat, open(fjson, "w"), indent=1)
fout = os.path.join(res, f"octave_out_{n_days}d_{tag}.mat")
if os.path.exists(fout):
    raise SystemExit(f"{fout} exists; refusing to overwrite (use a new tag)")
t0 = time.time()
r = subprocess.run(["octave", "--no-gui", "-q", os.path.join(ROOT, "scripts/octave_AN_v1/run_biowenp_octave_AN_v1.m"),
                    SNAP, fcsv, fjson, fout], capture_output=True, text=True)
t_oct = time.time() - t0
print(r.stdout[-2000:], r.stderr[-3000:])
if r.returncode != 0:
    raise SystemExit("octave failed")
O = loadmat(fout, squeeze_me=True)
t0 = time.time()
hyd = run_hydrology(HydroForcing(f["precip"], f["temp"], f["month"]), h)
po = run_p_faithful(hyd, h, p, dx=1, macro_percent=1, sorption_mode=2)
t_py = time.time() - t0
rows = []
for k in ["s", "I", "ET", "L", "Sp", "Qo", "Qd", "Sr", "Q_Ex", "Q_in", "Q_in_d", "C_TP_in", "Ss"]:
    a, b = np.asarray(O[k], float), hyd[k]
    rows.append((k, "hydrology", a.size, np.max(np.abs(a - b)), np.max(np.abs(a)), ))
pmap = {k: k for k in ["C_macro", "C_macro_POP", "C_macro_PIP", "qe_ad", "qe_mineral", "qe_Recalcitrant", "M_POP_Labile",
                        "M_POP_Stable", "Plant_P", "TPflux", "TDPflux", "Cout_TP", "Cout_TDP", "C_Sr_TP", "Ce_TDP_FM",
                        "Ce_PP_FM", "TPExflux", "Precip_flux", "Adsorbed_flux"]}
for k in pmap:
    a, b = np.asarray(O[k], float), po[k]
    if a.shape != b.shape:
        a = a.reshape(b.shape)
    rows.append((k, "phosphorus", a.size, np.max(np.abs(a - b)), np.max(np.abs(a))))
df = pd.DataFrame(rows, columns=["variable", "module", "n_values", "max_abs_diff", "max_abs_value"])
df["max_rel_diff"] = df.max_abs_diff / df.max_abs_value.replace(0, np.nan)
out_csv = os.path.join(tab, f"equivalence_octave_vs_python_{n_days}d_{tag}.csv")
df.to_csv(out_csv, index=False)
print(df.to_string())
wb, wf = water_balance_original(hyd, h), water_balance_full(hyd, h)
pb = p_balance_original(hyd, po, h, p)
summ = dict(n_days=n_days, steps=len(f["precip"]), octave_seconds=t_oct, python_seconds=t_py,
            worst_rel_diff=float(np.nanmax(df.max_rel_diff)), water_original=wb, water_full=wf, p_original=pb)
json.dump(summ, open(os.path.join(tab, f"equivalence_summary_{n_days}d_{tag}.json"), "w"), indent=1, default=float)
print(json.dumps(summ, indent=1, default=float))
np.savez_compressed(os.path.join(res, f"python_faithful_{n_days}d_{tag}.npz"), **{k: v for k, v in hyd.items()},
                    **{"P_" + k: v for k, v in po.items() if np.ndim(v) <= 2})
