#!/usr/bin/env python3
# Part of biowenp_py_AN_v1 (LGPL-3.0-or-later).
"""Cross-implementation consistency check: the faithful port of the ORIGINAL CODE vs the Soil P thread's independent
reimplementation from the paper's Table 4/6 (worked_examples/03_parsons_zhou2023, no-uptake branch), under the
same SYNTHETIC configuration: saturated column (s = 1), I = 0.1 m/d, inflow TP 4.3 mg/L, rho_fm = 1.5 kg/L,
no ET, uptake off, 365 days, bedding drain-only. Differences measure code-vs-paper interpretation, not error."""
import os, sys, json
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from biowenp_AN_v1 import HydroParams, PParams, run_p_faithful
SOILP = "/workspace/soilp/swatplusP_sectionI_AN_v1/upstream_readonly/worked_examples/03_parsons_zhou2023/results_v1.csv"
days = 365
M = days * 144
h = HydroParams(rho_fm=1.5)
p = PParams(I_Pu_max=0.0)
I = np.full(M, 0.1 / 144)
A = h.A_BC * 1000
hyd = dict(s=np.ones(M), ET=np.zeros(M), I=I, L=I.copy(), Sr=np.zeros(M), Qd=I * A, Qo=np.zeros(M),
           Q_Ex=np.zeros(M), C_TP_in=np.full(M, 4.3))
po = run_p_faithful(hyd, h, p)
inflow = 4.3 * I.sum() * A                                   # mg
media_gain = (po["TPP"][:, -1].mean() - po["TPP"][:, 0].mean()) * h.dfm * A * h.rho_fm
d_s = p.d_s(h)
bed = np.sum(po["Ce_PP_FM"] * (1 - np.exp(-p.lambda_s * d_s)) * hyd["L"]) * A
out = np.sum(po["TPflux"])
pw = (po["C_macro"][:, -1] + po["C_macro_PP"][:, -1]).sum() * h.dfm * A * h.n_fm / po["C_macro"].shape[0]
g = lambda X: (X[:, -1].mean() - X[:, 0].mean())
gains = dict(adsorbed=g(po["qe_ad"]), mineral=g(po["qe_mineral"]), recalcitrant=g(po["qe_Recalcitrant"]),
             labile=g(po["M_POP_Labile"]), stable=g(po["M_POP_Stable"]))
tot = sum(gains.values())
port = {"media_gain_pct_of_inflow": 100 * media_gain / inflow, "bed_accum_pct_of_inflow": 100 * bed / inflow,
        "outflow_pct_of_inflow": 100 * out / inflow,
        "max_adsorbed_frac_of_Qmax": po["qe_ad"][:, -1].max() / p.Qmax_ad,
        **{f"media_share_{k}_pct": 100 * v / tot for k, v in gains.items()},
        "closure_residual_pct": 100 * (inflow - media_gain - bed - out - pw) / inflow}
ref = pd.read_csv(SOILP).set_index("quantity")["value"]
rows = []
for k, v in port.items():
    rk = "nouptake_" + k if not k.startswith("closure") else "nouptake_mass_balance_residual_pct"
    rows.append(dict(quantity=k, faithful_code_port=v, soilp_table4_reimpl=float(ref.get(rk, np.nan))))
df = pd.DataFrame(rows)
df["difference"] = df.faithful_code_port - df.soilp_table4_reimpl
fn = os.path.join(ROOT, "tables_AN_v1", "crosscheck_vs_soilp_table4_reimpl_AN_v1.csv")
df.to_csv(fn, index=False)
print(df.to_string())
