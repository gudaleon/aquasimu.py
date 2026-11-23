"""Post-process DREAM and PEST++ results on the SYNTHETIC case: recovery tables, R-hat,
credible intervals on kg N / kg P retained per year, PEST vs DREAM comparison table + figures,
and the mass-balance table.  Everything shown is a consistency check on synthetic data."""
import sys, pathlib, glob, re, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import multiprocessing as mp
HERE = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import pyemu
from bioreactor_AN_v1 import synthetic_AN_v1 as S
from bioreactor_AN_v1 import calib_common_AN_v1 as C
from bioreactor_AN_v1.dream_AN_v1 import gelman_rubin

R, FIG, TAB = HERE / "results_AN_v1", HERE / "figures_AN_v1", HERE / "tables_AN_v1"
FIG.mkdir(exist_ok=True); TAB.mkdir(exist_ok=True)
rng = np.random.default_rng(5)
NP_ = len(S.NAMES)
q = lambda a, p: np.percentile(a, p, axis=0)


def fwd(z):
    yN, yP, kN, kP = C.forward(z)
    return np.r_[yN, yP, kN, kP]


def dream_summary(file, d):
    r = np.load(R / file)
    ch = r["chains"]; n = ch.shape[0] - 1
    half = ch[n // 2:, :, :d]
    samp = half.reshape(-1, d)
    rhat = gelman_rubin(half)
    thin = samp[rng.choice(len(samp), 2000, replace=False), :]
    with mp.get_context("fork").Pool(8) as p:
        sims = np.array(p.map(fwd, list(thin[:, :NP_])))
    nobs = len(C.D)
    return dict(raw=r, samp=samp, rhat=rhat, rhat_all=gelman_rubin(ch[n // 2:]), sims=sims, thin=thin,
                yN=sims[:, :nobs], yP=sims[:, nobs:2 * nobs], kgN=sims[:, -2], kgP=sims[:, -1], accept=float(r["accept_rate"]),
                n_gen=n)


D1 = dream_summary("dream_fixed_sigma_AN_v1.npz", NP_)
D2 = dream_summary("dream_hetero_sigma_AN_v1.npz", NP_ + 4)

# ------------------------------- PEST++ GLM ------------------------------------------------------------------
g = HERE / "pest_AN_v1" / "run_glm"
pst = pyemu.Pst(str(g / "pest_AN_v1.pst")); pst.parrep(str(g / "pest_AN_v1.par"))
jco = str(g / "pest_AN_v1.jco")
sc = pyemu.Schur(jco=jco, pst=pst, forecasts=["fkgn", "fkgp"])
post = sc.posterior_parameter; names_l = [n.lower() for n in S.NAMES]
pcov = pd.DataFrame(post.x, index=post.names, columns=post.names).loc[names_l, names_l]
psd = np.sqrt(np.diag(pcov))                                          # sd in PEST estimation space (log10 for log pars)
pcorr = pcov / np.outer(psd, psd)
prior_sd = np.sqrt(np.diag(sc.parcov.to_dataframe().loc[names_l, names_l]))
best = pst.parameter_data.loc[names_l, "parval1"].values
best_z = S.to_internal(best)
glm_fc = pst.res.loc[["fkgn", "fkgp"], "modelled"].values if pst.res is not None else None
fs = sc.get_forecast_summary()
glm_sd_fc = np.sqrt(fs.loc[["fkgn", "fkgp"], "post_var"].values)
glm_prior_sd_fc = np.sqrt(fs.loc[["fkgn", "fkgp"], "prior_var"].values)
ev = pyemu.ErrVar(jco=jco, pst=pst, forecasts=["fkgn", "fkgp"])
ident5 = ev.get_identifiability_dataframe(singular_value=5)["ident"].loc[names_l].values
ident3 = ev.get_identifiability_dataframe(singular_value=3)["ident"].loc[names_l].values
sv = np.linalg.svd(sc.xtqx.x, compute_uv=False)
glm_phi = float(pst.phi) if hasattr(pst, "phi") else np.nan
rec = open(g / "pest_AN_v1.rec").read()
glm_iters = len(re.findall(r"OPTIMISATION ITERATION NUMBER", rec))

# ------------------------------- PEST++ IES ------------------------------------------------------------------
ies_dir = HERE / "pest_AN_v1" / "run_ies"
its = sorted(int(re.search(r"\.(\d+)\.par\.csv$", f).group(1)) for f in glob.glob(str(ies_dir / "pest_AN_v1.*.par.csv")))
it = its[-1]
pe = pd.read_csv(ies_dir / f"pest_AN_v1.{it}.par.csv", index_col=0)
oe = pd.read_csv(ies_dir / f"pest_AN_v1.{it}.obs.csv", index_col=0)
phi = pd.read_csv(ies_dir / "pest_AN_v1.phi.actual.csv")
ies_phi_med = float(np.median(phi.iloc[-1][[c for c in phi.columns if c.isdigit()]].astype(float)))
ies_z = S.to_internal(pe[names_l].values.T).T if False else np.column_stack([np.where(S.LOG[i], np.log10(pe[n].values), pe[n].values) for i, n in enumerate(names_l)])
ies_kgN, ies_kgP = oe["fkgn"].values, oe["fkgp"].values

# ------------------------------- tables ------------------------------------------------------------------
def ci(a): return q(a, [2.5, 50, 97.5])
rows = []
for i, n in enumerate(S.NAMES):
    tru = S.TRUE[i]; phys = lambda zz: S.to_physical(np.where(np.arange(NP_) == i, zz, 0 * S.Z_TRUE))[i] if False else (10 ** zz if S.LOG[i] else zz)
    row = dict(parameter=n, true=tru)
    for tag, Dd in (("DREAM_fixedSigma", D1), ("DREAM_heteroSigma", D2)):
        lo_, md, hi_ = ci(Dd["samp"][:, i])
        row.update({f"{tag}_median": phys(md), f"{tag}_2.5%": phys(lo_), f"{tag}_97.5%": phys(hi_),
                    f"{tag}_Rhat": Dd["rhat"][i], f"{tag}_true_in_95CI": bool(phys(lo_) <= tru <= phys(hi_))})
    zlo, zhi = best_z[i] - 1.96 * psd[i], best_z[i] + 1.96 * psd[i]
    row.update({"PESTGLM_best": best[i], "PESTGLM_lin2.5%": phys(zlo), "PESTGLM_lin97.5%": phys(zhi),
                "PESTGLM_true_in_95CI": bool(phys(zlo) <= tru <= phys(zhi))})
    lo_, md, hi_ = ci(ies_z[:, i])
    row.update({"PESTIES_median": phys(md), "PESTIES_2.5%": phys(lo_), "PESTIES_97.5%": phys(hi_),
                "PESTIES_true_in_95CI": bool(phys(lo_) <= tru <= phys(hi_))})
    rows.append(row)
rec_tab = pd.DataFrame(rows)
rec_tab.to_csv(TAB / "parameter_recovery_AN_v1.csv", index=False, float_format="%.5g")

# error-model parameters of the hetero run
er = []
for j, n in enumerate(C.ERR_NAMES):
    s = D2["samp"][:, NP_ + j]; lo_, md, hi_ = ci(s)
    f = (lambda v: 10 ** v) if n.endswith("0") else (lambda v: v)
    er.append(dict(parameter=n, true=f(C.ERR_TRUE[j]), median=f(md), **{"2.5%": f(lo_), "97.5%": f(hi_)}, Rhat=D2["rhat"][NP_ + j]))
pd.DataFrame(er).to_csv(TAB / "error_model_recovery_AN_v1.csv", index=False, float_format="%.5g")

# forecast table
fc = []
tr = {"kgN_retained_per_yr": S_true_kg for S_true_kg in [C.OBS["true_kgN"]]}
for lab, key, true_v, glm_i, ies in (("kg N retained / yr (woodchip)", "kgN", C.OBS["true_kgN"], 0, ies_kgN),
                                      ("kg P retained / yr (eggshell)", "kgP", C.OBS["true_kgP"], 1, ies_kgP)):
    row = {"quantity": lab, "truth_(synthetic)": true_v}
    for tag, Dd in (("DREAM_fixedSigma", D1), ("DREAM_heteroSigma", D2)):
        lo_, md, hi_ = ci(Dd[key]); row.update({f"{tag}_median": md, f"{tag}_2.5%": lo_, f"{tag}_97.5%": hi_, f"{tag}_95CI_width": hi_ - lo_,
                                               f"{tag}_truth_in_CI": bool(lo_ <= true_v <= hi_)})
    m_, s_ = glm_fc[glm_i], glm_sd_fc[glm_i]
    row.update({"PESTGLM_best": m_, "PESTGLM_lin2.5%": m_ - 1.96 * s_, "PESTGLM_lin97.5%": m_ + 1.96 * s_, "PESTGLM_95CI_width": 3.92 * s_,
                "PESTGLM_truth_in_CI": bool(abs(m_ - true_v) <= 1.96 * s_), "PESTGLM_prior_sd": glm_prior_sd_fc[glm_i]})
    lo_, md, hi_ = ci(ies); row.update({"PESTIES_median": md, "PESTIES_2.5%": lo_, "PESTIES_97.5%": hi_, "PESTIES_95CI_width": hi_ - lo_,
                                       "PESTIES_truth_in_CI": bool(lo_ <= true_v <= hi_)})
    fc.append(row)
fc_tab = pd.DataFrame(fc); fc_tab.to_csv(TAB / "forecast_kg_retained_AN_v1.csv", index=False, float_format="%.5g")

# one-table comparison: parameter sd (log10 units / linear for theta), correlation, identifiability
dsd = D1["samp"].std(axis=0)[:NP_]; dsd2 = D2["samp"][:, :NP_].std(axis=0)
prior_u = (S.Z_HI - S.Z_LO) / np.sqrt(12)
cmp_tab = pd.DataFrame({"parameter": S.NAMES, "sd_DREAM_fixed": dsd, "sd_DREAM_hetero": dsd2, "sd_PESTGLM_linear": psd,
                        "sd_PESTIES": ies_z.std(axis=0), "prior_sd_uniform": prior_u,
                        "DREAM_info_gain(1-sd/prior)": 1 - dsd / prior_u,
                        "PEST_identifiability_5SV": ident5, "PEST_identifiability_3SV": ident3,
                        "PEST_var_reduction_%": 100 * (1 - psd**2 / prior_sd**2)})
cmp_tab.to_csv(TAB / "pest_vs_dream_parameters_AN_v1.csv", index=False, float_format="%.5g")
dcorr = np.corrcoef(D1["samp"][:, :NP_].T); icorr = np.corrcoef(ies_z.T)
pd.DataFrame(dcorr, index=S.NAMES, columns=S.NAMES).to_csv(TAB / "correlation_DREAM_AN_v1.csv", float_format="%.3f")
pd.DataFrame(pcorr.values, index=S.NAMES, columns=S.NAMES).to_csv(TAB / "correlation_PESTGLM_AN_v1.csv", float_format="%.3f")
pd.DataFrame(icorr, index=S.NAMES, columns=S.NAMES).to_csv(TAB / "correlation_PESTIES_AN_v1.csv", float_format="%.3f")

# ------------------------------- mass-balance table (true-parameter run) ---------------------------------------
out = S.run_model(S.TRUE, C.FORCING)
mbt = pd.concat([r.balance_table() for r in out["results"]])
ann = []
for r in out["results"]:
    b = r.check_balance(1e-9)
    ann.append(dict(unit=r.unit.name, period="annual", water_in_m3=b["water_in_m3"], water_treated_out_m3=b["water_treated_out_m3"],
                    water_bypass_m3=b["water_bypass_m3"], water_residual_m3=b["water_residual_m3"],
                    N_in_kg=b["N_in_g"] / 1e3, N_out_treated_kg=b["N_out_treated_g"] / 1e3, N_bypass_kg=b["N_bypass_g"] / 1e3,
                    N_removed_kg=b["N_removed_g"] / 1e3, N_storage_change_kg=b["N_storage_change_g"] / 1e3, N_residual_kg=b["N_residual_g"] / 1e3,
                    P_in_kg=b["P_in_g"] / 1e3, P_out_treated_kg=b["P_out_treated_g"] / 1e3, P_bypass_kg=b["P_bypass_g"] / 1e3,
                    P_removed_kg=b["P_removed_g"] / 1e3, P_storage_change_kg=b["P_storage_change_g"] / 1e3, P_residual_kg=b["P_residual_g"] / 1e3))
mbt = mbt.reset_index().rename(columns={"index": "period"})
mbt["period"] = mbt["period"].astype(str).str[:7]
mbt = pd.concat([mbt, pd.DataFrame(ann)], ignore_index=True)
mbt.to_csv(TAB / "mass_balance_table_AN_v1.csv", index=False, float_format="%.6g")
# worst residual over ALL posterior-sample runs (each of which passed the 1e-6 closure test inside run_model)
worst = 0.0
for z in D1["thin"][:200, :NP_]:
    tb = S.run_model(S.to_physical(z), C.FORCING)["train_balance"]
    worst = max(worst, tb["N_rel_residual"], tb["P_rel_residual"], tb["water_rel_residual"])
summary = dict(true_run_max_rel_residual=float(max(out["train_balance"][k] for k in ("N_rel_residual", "P_rel_residual", "water_rel_residual"))),
               worst_rel_residual_over_200_posterior_runs=float(worst), tolerance=1e-6,
               dream_fixed=dict(n_gen=D1["n_gen"], chains=8, model_runs=D1["n_gen"] * 8, accept=D1["accept"], max_Rhat=float(D1["rhat"].max())),
               dream_hetero=dict(n_gen=D2["n_gen"], chains=8, model_runs=D2["n_gen"] * 8, accept=D2["accept"], max_Rhat=float(D2["rhat"].max())),
               pest_glm=dict(iterations=glm_iters, final_phi=glm_phi, n_obs_weighted=int(2 * len(C.D)),
                             singular_values=[float(x) for x in sv]),
               pest_ies=dict(final_iteration=it, n_realisations=int(len(pe)), median_phi=float(ies_phi_med)),
               truth=dict(kgN=C.OBS["true_kgN"], kgP=C.OBS["true_kgP"]))
json.dump(summary, open(TAB / "run_summary_AN_v1.json", "w"), indent=1)

# ------------------------------- figures --------------------------------------------------------------------------
lab = [f"{n}\n(log10)" if S.LOG[i] else n for i, n in enumerate(S.NAMES)]
# F1 posterior pairs
fig, ax = plt.subplots(NP_, NP_, figsize=(13, 12))
s = D1["samp"][:, :NP_][::10]
for i in range(NP_):
    for j in range(NP_):
        a = ax[i, j]
        if j > i: a.axis("off"); continue
        if i == j:
            a.hist(s[:, i], 40, color="tab:blue", alpha=.7); a.axvline(S.Z_TRUE[i], color="r", lw=2); a.axvline(best_z[i], color="k", ls="--")
        else:
            a.plot(s[:, j], s[:, i], ".", ms=1.5, alpha=.3); a.plot(S.Z_TRUE[j], S.Z_TRUE[i], "r*", ms=12)
            a.plot(best_z[j], best_z[i], "k^", ms=7)
        if i == NP_ - 1: a.set_xlabel(lab[j], fontsize=8)
        if j == 0 and i > 0: a.set_ylabel(lab[i], fontsize=8)
fig.suptitle("DREAM(ZS) posterior (SYNTHETIC data, fixed-sigma run).  red = true value, black = PEST-GLM best fit", y=0.995)
fig.tight_layout(); fig.savefig(FIG / "fig_dream_posterior_AN_v1.png", dpi=120); plt.close(fig)

# F2 convergence
fig, ax = plt.subplots(2, 3, figsize=(16, 8))
for Dd, nm, c in ((D1, "fixed sigma", "tab:blue"), (D2, "hetero sigma", "tab:orange")):
    rt = Dd["raw"]["rhat_trace"]; g_ = Dd["raw"]["rhat_gens"]
    ax[0, 0].plot(g_, rt.max(axis=1), color=c, label=f"max R-hat, {nm}")
ax[0, 0].axhline(1.2, color="grey", ls=":"); ax[0, 0].axhline(1.05, color="grey", ls="--"); ax[0, 0].set_ylim(0.98, 3); ax[0, 0].set_yscale("log")
ax[0, 0].set_xlabel("generation"); ax[0, 0].set_title("Gelman-Rubin R-hat (last 50% of chains)"); ax[0, 0].legend()
rt = D1["raw"]["rhat_trace"]
for i, n in enumerate(S.NAMES): ax[0, 1].plot(D1["raw"]["rhat_gens"], rt[:, i], label=n)
ax[0, 1].set_ylim(0.98, 2); ax[0, 1].axhline(1.05, color="grey", ls="--"); ax[0, 1].legend(fontsize=7); ax[0, 1].set_title("R-hat per parameter (fixed-sigma run)")
ax[0, 2].plot(D1["raw"]["logp"].max(axis=1)); ax[0, 2].set_title("best log-posterior (fixed-sigma run)"); ax[0, 2].set_ylim(D1["raw"]["logp"].max() - 60, D1["raw"]["logp"].max() + 2)
chn = D1["raw"]["chains"]
for k, i in enumerate((0, 3, 4)):
    for c_ in range(8): ax[1, k].plot(chn[:, c_, i], lw=.4)
    ax[1, k].axhline(S.Z_TRUE[i], color="k", ls="--"); ax[1, k].axvline(D1["n_gen"] // 2, color="r", ls=":"); ax[1, k].set_title(f"chains: {lab[i].replace(chr(10), ' ')}")
fig.suptitle("DREAM(ZS) convergence (SYNTHETIC case); dotted red = start of post-burn-in half")
fig.tight_layout(); fig.savefig(FIG / "fig_dream_convergence_AN_v1.png", dpi=120); plt.close(fig)

# F3 fit + predictive + kg
fig = plt.figure(figsize=(16, 9)); gs = fig.add_gridspec(2, 3)
for k, (nm, ob, tru, sim, key) in enumerate((("NO3-N (g/m3)", C.OBS["obsN"], C.OBS["trueN"], D1["yN"], "N"), ("PO4-P (g/m3)", C.OBS["obsP"], C.OBS["trueP"], D1["yP"], "P"))):
    a = fig.add_subplot(gs[k, :2]); lo_, md, hi_ = q(sim, [2.5, 50, 97.5])
    a.fill_between(C.D, lo_, hi_, color="tab:blue", alpha=.35, label="DREAM 95% (parameter uncertainty)")
    a.plot(C.D, md, color="tab:blue", lw=1, label="DREAM median"); a.plot(C.D, tru, "k-", lw=.8, label="truth (noise-free)")
    a.plot(C.D, ob, "r.", ms=4, label="synthetic observations")
    a.set_title(f"Train outlet {nm}  -  SYNTHETIC data"); a.set_xlabel("day of year"); a.legend(fontsize=7, ncol=4)
for k, (key, tv, gi, ies_v) in enumerate((("kgN", C.OBS["true_kgN"], 0, ies_kgN), ("kgP", C.OBS["true_kgP"], 1, ies_kgP))):
    a = fig.add_subplot(gs[k, 2])
    a.hist(D1[key], 50, density=True, alpha=.6, label="DREAM fixed-sigma"); a.hist(D2[key], 50, density=True, alpha=.5, label="DREAM hetero-sigma")
    xs = np.linspace(*a.get_xlim(), 300); m_, s_ = glm_fc[gi], glm_sd_fc[gi]
    a.plot(xs, np.exp(-.5 * ((xs - m_) / s_) ** 2) / (s_ * np.sqrt(2 * np.pi)), "g-", label="PEST-GLM linear")
    a.hist(ies_v, 25, density=True, histtype="step", color="m", lw=2, label="PEST-IES"); a.axvline(tv, color="k", ls="--", label="truth")
    a.set_title("kg N retained per year" if key == "kgN" else "kg P retained per year"); a.legend(fontsize=6)
fig.tight_layout(); fig.savefig(FIG / "fig_fit_and_forecast_AN_v1.png", dpi=120); plt.close(fig)

# F4 PEST vs DREAM comparison
fig, ax = plt.subplots(2, 3, figsize=(17, 10))
a = ax[0, 0]; x = np.arange(NP_); wd = .2
for k, (tag, mid, lo_, hi_, c) in enumerate((
        ("DREAM fixed", np.median(D1["samp"][:, :NP_], 0), q(D1["samp"][:, :NP_], 2.5), q(D1["samp"][:, :NP_], 97.5), "tab:blue"),
        ("DREAM hetero", np.median(D2["samp"][:, :NP_], 0), q(D2["samp"][:, :NP_], 2.5), q(D2["samp"][:, :NP_], 97.5), "tab:orange"),
        ("PEST-GLM (linear)", best_z, best_z - 1.96 * psd, best_z + 1.96 * psd, "tab:green"),
        ("PEST-IES", np.median(ies_z, 0), q(ies_z, 2.5), q(ies_z, 97.5), "m"))):
    a.errorbar(x + (k - 1.5) * wd, mid - S.Z_TRUE, yerr=[mid - lo_, hi_ - mid], fmt="o", color=c, capsize=3, label=tag)
a.axhline(0, color="k"); a.set_xticks(x); a.set_xticklabels(lab, fontsize=8); a.set_ylabel("estimate - truth (log10 for log params)")
a.set_title("Parameter estimates with 95% intervals"); a.legend(fontsize=7)
for a, M, t_ in ((ax[0, 1], dcorr, "DREAM posterior correlation"), (ax[0, 2], pcorr.values, "PEST-GLM linear posterior correlation")):
    im = a.imshow(M, vmin=-1, vmax=1, cmap="RdBu_r"); a.set_xticks(x); a.set_yticks(x); a.set_xticklabels(S.NAMES, rotation=45, fontsize=8); a.set_yticklabels(S.NAMES, fontsize=8)
    for i in range(NP_):
        for j in range(NP_): a.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7)
    a.set_title(t_)
a = ax[1, 0]
a.bar(x - .2, cmp_tab["DREAM_info_gain(1-sd/prior)"], .2, label="DREAM: 1 - sd_post/sd_prior")
a.bar(x, ident5, .2, label="PEST identifiability (5 SV)"); a.bar(x + .2, ident3, .2, label="PEST identifiability (3 SV)")
a.set_xticks(x); a.set_xticklabels(S.NAMES, fontsize=8); a.set_ylim(0, 1.05); a.set_title("Identifiability / information gain"); a.legend(fontsize=7)
for k, (a, key, tv, gi, ies_v, t_) in enumerate(((ax[1, 1], "kgN", C.OBS["true_kgN"], 0, ies_kgN, "kg N retained / yr"),
                                                   (ax[1, 2], "kgP", C.OBS["true_kgP"], 1, ies_kgP, "kg P retained / yr"))):
    items = [("DREAM fixed", ci(D1[key]), "tab:blue"), ("DREAM hetero", ci(D2[key]), "tab:orange"),
             ("PEST-GLM", np.array([glm_fc[gi] - 1.96 * glm_sd_fc[gi], glm_fc[gi], glm_fc[gi] + 1.96 * glm_sd_fc[gi]]), "tab:green"),
             ("PEST-IES", ci(ies_v), "m")]
    for j, (nm, (l_, m_, h_), c) in enumerate(items):
        a.errorbar(m_, j, xerr=[[m_ - l_], [h_ - m_]], fmt="o", color=c, capsize=4)
    a.axvline(tv, color="k", ls="--", label="truth"); a.set_yticks(range(4)); a.set_yticklabels([i[0] for i in items]); a.set_title(t_ + " (95% interval)"); a.legend()
fig.suptitle("PEST++ vs DREAM(ZS) on SYNTHETIC data: estimates, correlation, identifiability, forecast uncertainty (consistency check)")
fig.tight_layout(); fig.savefig(FIG / "fig_pest_vs_dream_AN_v1.png", dpi=120); plt.close(fig)

# F5 mass balance (true-parameter run)
w = mbt[(mbt.unit == "woodchip") & (mbt.period != "annual")]; e = mbt[(mbt.unit == "eggshell") & (mbt.period != "annual")]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
x_ = np.arange(len(w))
ax[0].bar(x_, w.N_out_treated_kg, label="treated outflow"); ax[0].bar(x_, w.N_bypass_kg, bottom=w.N_out_treated_kg, label="bypass")
ax[0].bar(x_, w.N_removed_kg, bottom=w.N_out_treated_kg + w.N_bypass_kg, label="removed (denitrified)"); ax[0].plot(x_, w.N_in_kg, "k_", ms=18, mew=2, label="inflow")
ax[0].set_title("Woodchip unit: monthly N balance (kg)"); ax[0].legend(fontsize=7)
ax[1].bar(x_, e.P_out_treated_kg, label="treated outflow"); ax[1].bar(x_, e.P_bypass_kg, bottom=e.P_out_treated_kg, label="bypass")
ax[1].bar(x_, e.P_removed_kg, bottom=e.P_out_treated_kg + e.P_bypass_kg, label="sorbed on eggshell"); ax[1].plot(x_, e.P_in_kg, "k_", ms=18, mew=2, label="inflow")
ax[1].set_title("Eggshell unit: monthly P balance (kg)"); ax[1].legend(fontsize=7)
ax[2].semilogy(x_, np.abs(w.N_residual_kg) + 1e-20, "o-", label="woodchip N"); ax[2].semilogy(x_, np.abs(e.P_residual_kg) + 1e-20, "s-", label="eggshell P")
ax[2].semilogy(x_, np.abs(w.water_residual_m3) + 1e-20, "^-", label="woodchip water (m3)"); ax[2].set_title("closure residuals |in - out - bypass - removed - storage|"); ax[2].legend(fontsize=7)
ax[2].set_ylim(1e-20, 1e-3)
fig.tight_layout(); fig.savefig(FIG / "fig_mass_balance_AN_v1.png", dpi=120)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
print(rec_tab.T.to_string()); print(fc_tab.T.to_string()); print(cmp_tab.to_string()); print(json.dumps(summary, indent=1))
