#!/usr/bin/env python3
# Part of biowenp_py_AN_v1 (LGPL-3.0-or-later).
"""Demo of the augmented flow-through model for eggshell + woodchip reactors on SYNTHETIC forcing (3 years).
Every number and figure produced here comes from invented (seeded) data and literature PRIORS. Nothing is field data,
nothing is calibrated, nothing is validated. Writes tables_AN_v1/demo_* and figures_AN_v1/demo_*."""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from biowenp_AN_v1.flowthrough_AN_v1 import FlowThroughUnit, ReactorForcing, simulate_train
from biowenp_AN_v1.synthetic_AN_v1 import synthetic_reactor_daily
from biowenp_AN_v1.priors_AN_v1 import eggshell_priors, woodchip_priors, prior_table, LEE2022

TAB, FIG = os.path.join(ROOT, "tables_AN_v1"), os.path.join(ROOT, "figures_AN_v1")
SYN = "SYNTHETIC data - not field data, not validation"
years = [synthetic_reactor_daily(365, seed=11 + k) for k in range(3)]
s = {k: np.concatenate([y[k] for y in years]) for k in years[0]}
f = ReactorForcing(s["Q"], s["T"], s["C_N"], s["C_SRP"], s["C_PP"])
f.table().to_csv(os.path.join(TAB, "demo_input_forcing_SYNTHETIC_AN_v1.csv"), index=False)
pd.DataFrame(prior_table()).to_csv(os.path.join(TAB, "demo_priors_Lee2022_Addy2016_AN_v1.csv"), index=False)
wp = woodchip_priors(0.65)
ep = eggshell_priors(0.04)
wood = FlowThroughUnit(name="woodchip_bed", n_cells=4, bed_volume=50, porosity=0.65, length=12, q_design=50,
                       alpha_L=0.3, n_kinetics="mm", vmax_N20=wp["vmax_N20"][0], K_N=wp["K_N"][0],
                       theta_N=wp["theta_N"][0])
bed_mass = 20 * 350.0                     # kg mixed bed (20 m3 at an assumed 350 kg/m3)
egg_mass = 0.04 * bed_mass                # 4 wt% calcined eggshell (proposal: 3-5 wt%)
scen = {"q_max 11.5 g/kg (Lee2022 column)": 11.5, "q_max 35 g/kg (prior median)": ep["q_max"][0],
        "q_max 108.2 g/kg (Lee2022 batch)": 108.2}
rows, results = [], {}
for label, qmax in scen.items():
    egg = FlowThroughUnit(name="eggshell_woodchip_mix", n_cells=3, bed_volume=20, porosity=0.55, length=6,
                          q_design=50, alpha_L=0.2, n_kinetics="mm", vmax_N20=0.5 * wp["vmax_N20"][0],
                          K_N=wp["K_N"][0], theta_N=wp["theta_N"][0], p_kinetics="zhou", sorbent_mass=egg_mass,
                          q_max=qmax, K_L=ep["K_L"][0], k_ad=ep["k_ad"][0], theta_P=ep["theta_P"][0],
                          lambda_pp=0.3)
    res, train = simulate_train([wood, egg], f, rtol=1e-6)        # gate on every unit + train
    results[label] = (res, train)
    for r in res:
        bt = r.balance_table(); bt.insert(0, "unit", r.unit.name); bt.insert(0, "scenario", label)
        rows.append(bt)
    if label.startswith("q_max 35"):
        wood.table().to_csv(os.path.join(TAB, "demo_input_unit_woodchip_AN_v1.csv"), index=False)
        egg.table().to_csv(os.path.join(TAB, "demo_input_unit_eggshell_mix_AN_v1.csv"), index=False)
bal = pd.concat(rows)
bal.to_csv(os.path.join(TAB, "demo_mass_balance_and_residuals_SYNTHETIC_AN_v1.csv"), index=False)
summary = []
for label, (res, train) in results.items():
    w, e = res
    summary.append(dict(scenario=label, kgN_retained_per_yr_train=(w.ledger["N_rem"] + e.ledger["N_rem"]) / 3000,
                        kgP_retained_per_yr_train=(w.ledger["P_rem"] + e.ledger["P_rem"]) / 3000,
                        P_in_kg=train["P_in"] / 1000, P_out_kg=train["P_out"] / 1000,
                        P_bypass_kg=e.ledger["P_b"] / 1000, P_sorbed_kg=e.ledger["P_sorbed"] / 1000,
                        P_filtered_kg=e.ledger["P_filtered"] / 1000,
                        final_capacity_used_frac=e.ledger["final_capacity_used_frac"],
                        max_rel_residual=max(abs(a / b) for r in res for (a, b) in r.residuals().values()),
                        train_rel_residual_N=train["residual_N"] / train["N_in"],
                        train_rel_residual_P=train["residual_P"] / train["P_in"]))
sm = pd.DataFrame(summary)
sm.to_csv(os.path.join(TAB, "demo_summary_SYNTHETIC_AN_v1.csv"), index=False)
print(sm.to_string())

# ---------- figures
lab35 = list(scen)[1]
w, e = results[lab35][0]
days = np.arange(len(f.Q))
fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
ax[0].plot(days, f.Q, lw=0.8, label="inflow Q"); ax[0].axhline(50, color="k", ls="--", lw=0.8, label="design flow (bypass above)")
ax[0].fill_between(days, 50, np.maximum(f.Q, 50), color="tab:red", alpha=0.3, label="bypass Q_b")
ax[0].set_ylabel("m$^3$/d"); ax[0].legend(fontsize=8)
ax[1].plot(days, f.C_N, lw=0.8, label="NO$_3$-N in"); ax[1].plot(days, w.daily.C_N_out, lw=0.8, label="after woodchip")
ax[1].plot(days, e.daily.C_N_out, lw=0.8, label="after eggshell-woodchip mix"); ax[1].set_ylabel("g N/m$^3$"); ax[1].legend(fontsize=8)
axb = ax[1].twinx(); axb.plot(days, f.T, color="grey", lw=0.5); axb.set_ylabel("T (°C)", color="grey")
ax[2].plot(days, f.C_SRP + f.C_PP, lw=0.8, color="tab:blue", label="TP in")
ax[2].plot(days, e.daily.C_P_out + e.daily.C_PP_out, lw=0.8, color="tab:orange", label="TP out of train")
ax[2].plot(days, e.daily.C_P_out, lw=0.8, color="tab:green", label="SRP out of train"); ax[2].set_ylabel("g P/m$^3$"); ax[2].set_xlabel("day"); ax[2].legend(fontsize=8)
fig.suptitle(f"Augmented BioWEN-P flow-through train, prior-median parameters — {SYN}", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "demo_timeseries_SYNTHETIC_AN_v2.png"), dpi=150); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for label, (res, _) in results.items():
    ax[0].plot(days, res[1].daily.q_frac_capacity, label=label)
    ax[1].plot(days, res[1].daily.C_P_out, lw=0.7, label=label)
ax[0].set_ylabel("sorbed P / q_max (eggshell)"); ax[0].set_xlabel("day"); ax[0].legend(fontsize=7)
ax[0].set_title("capacity used (Lee 2022 priors, not Zhou Qmax)", fontsize=9)
ax[1].set_ylabel("SRP leaving train (g/m$^3$)"); ax[1].set_xlabel("day"); ax[1].set_title("SRP breakthrough", fontsize=9)
fig.suptitle(SYN, fontsize=9); fig.tight_layout(); fig.savefig(os.path.join(FIG, "demo_capacity_depletion_SYNTHETIC_AN_v1.png"), dpi=150); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for i, (sp, terms) in enumerate((("water (m$^3$)", ["V_in", "V_t", "V_b"]),
                                 ("N (kg)", ["N_in", "N_out_t", "N_b", "N_rem", "dS_N"]),
                                 ("P (kg)", ["P_in", "P_out_t", "P_b", "P_sorbed", "P_mineral", "P_filtered", "dS_P"]))):
    sc = 1.0 if i == 0 else 1e-3
    x = np.arange(len(terms)); wdt = 0.38
    ax[i].bar(x - wdt / 2, [w.ledger[t] * sc for t in terms], wdt, label="woodchip bed")
    ax[i].bar(x + wdt / 2, [e.ledger[t] * sc for t in terms], wdt, label="eggshell-woodchip mix")
    ax[i].set_xticks(x); ax[i].set_xticklabels(terms, rotation=45, fontsize=8); ax[i].set_title(sp, fontsize=9)
    key = ["water", "N", "P"][i]
    rr = [abs(r.residuals()[key][0] / r.residuals()[key][1]) for r in (w, e)]
    ax[i].text(0.02, 0.95, "rel. residual: %.1e / %.1e" % tuple(rr), transform=ax[i].transAxes, fontsize=8, va="top")
ax[0].legend(fontsize=8)
fig.suptitle(f"Mass-balance gate (rtol 1e-6) — 3-yr totals, prior-median run — {SYN}", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "demo_mass_balance_SYNTHETIC_AN_v1.png"), dpi=150); plt.close(fig)
print("figures written")
