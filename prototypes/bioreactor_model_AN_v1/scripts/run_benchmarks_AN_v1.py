"""Consistency check of the Python equations against the real AQUASIM 2.1g (batch, built from
Eawag's public source) running Eawag's tutorial example files (see run_aquasim_examples_AN_v1.py).

This is a CONSISTENCY CHECK of equations/numerics against another implementation of the same
equations, not a validation against field data.  Nothing here is typed in from a plot: all
AQUASIM numbers are read from the listing files produced by running AQUASIM.
"""
import sys, pathlib
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
HERE = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
from bioreactor_AN_v1 import Unit, Forcing, simulate_unit
from bioreactor_AN_v1.aquasim_io_AN_v1 import read_listing

AQ = HERE / "benchmarks_AN_v1" / "aquasim_out"
rows, fig_data = [], {}


def metrics(case, what, t, y_py, y_aq, scale, bal=None, extra=None):
    err = y_py - y_aq
    d = dict(case=case, quantity=what, n_points=len(t), RMSE=float(np.sqrt(np.mean(err**2))),
             max_abs_err=float(np.max(np.abs(err))), scale=scale,
             RMSE_pct_of_scale=100 * float(np.sqrt(np.mean(err**2))) / scale,
             max_err_pct_of_scale=100 * float(np.max(np.abs(err))) / scale)
    if bal: d.update(bal)
    if extra: d.update(extra)
    rows.append(d)
    return d


def front_time(t, y, level):
    """first time y crosses `level` upward (linear interpolation)"""
    i = np.argmax(y >= level)
    if y[i] < level: return np.nan
    return float(t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))


def mb(r):
    b = r.check_balance(1e-6)
    return dict(N_rel_resid=b["N_rel_residual"], P_rel_resid=b["P_rel_residual"], water_rel_resid=b["water_rel_residual"])


# ---------------------------------------------------------------------------------------------
# B1  batch reactor, first-order decay (tutorial 2.1 part A: kA = 1/h, C0 = 1 mg/l, V = 10 l, 10 h)
# ---------------------------------------------------------------------------------------------
_, c = read_listing(AQ / "proc_a__Conc.txt"); t_aq, y_aq = c[0]; t_aq, y_aq = t_aq[1:], y_aq[1:]
f = Forcing(np.zeros(100), np.full(100, 20.0), np.zeros(100), np.zeros(100), dt=0.1)
r = simulate_unit(Unit(name="batch_first_order", volume=10, n_kinetics="first_order", k_n20=1.0, cN0=1.0), f)
y_py = r.series("cN_mean")
metrics("B1 batch first-order (proc_a)", "C_A(t) python vs AQUASIM", t_aq, y_py, y_aq, 1.0, mb(r))
metrics("B1 batch first-order (proc_a)", "C_A(t) AQUASIM vs analytic exp(-kt)", t_aq, y_aq, np.exp(-t_aq), 1.0)
fig_data["B1"] = (t_aq, y_aq, y_py, np.exp(-t_aq))

# B2  batch reactor, Michaelis-Menten (Monod) step B -> C: rmax=0.25 mg/l/h, K=0.5 mg/l, C_B(0)=1 (modified initial cond.)
_, c = read_listing(AQ / "procB_mm_mod__Conc.txt"); t_aq, y_aq = c[1][0][1:], c[1][1][1:]
r = simulate_unit(Unit(name="batch_mm", volume=10, n_kinetics="mm", vmax_n=0.25, K_n=0.5, cN0=1.0), f)
y_py = r.series("cN_mean")
metrics("B2 batch Monod (proc_b, C_B(0)=1)", "C_B(t) python vs AQUASIM", t_aq, y_py, y_aq, 1.0, mb(r))
fig_data["B2"] = (t_aq, y_aq, y_py)

# ---------------------------------------------------------------------------------------------
# B3  advective-dispersive column, tracer pulse (tutorial 5.2 part A): Q=1, A=1, L=1, D=0.1
# ---------------------------------------------------------------------------------------------
_, c = read_listing(AQ / "advdif_a__C_t.txt"); t_aq, y_aq = c[1]; t_aq, y_aq = t_aq[1:], y_aq[1:]
dt, nsm = 0.0005, 20
tm = (np.arange(int(3.0 / dt)) + 0.5) * dt
cin = np.interp(tm, [0, 0.01, 1.0, 1.01, 10], [0, 1, 1, 0, 0])
ff = Forcing(np.ones_like(tm), np.full_like(tm, 20.0), np.zeros_like(tm), cin, dt=dt)
res3 = {}
# AQUASIM uses a 'transmission' outlet BC (d2C/dx2 = 0, manual eq. 3.43b) and Danckwerts-type flux inlet BC.
for n, outlet in ((50, "transmission"), (100, "transmission"), (200, "transmission"), (200, "zero_gradient")):
    u = Unit(name=f"AD_column_n{n}_{outlet}", n_cells=n, volume=1.0, disp_ratio=0.1, scheme="central", outlet=outlet)
    r = simulate_unit(u, ff)
    cface = 1.5 * r.series("cP_last") - 0.5 * r.series("cP_last2") if outlet == "transmission" else r.series("cP_last")
    y_py = cface[nsm - 1::nsm]
    d = metrics(f"B3 AD tracer pulse (advdif_a) n={n} outlet={outlet}", "C(x=L,t) python vs AQUASIM", t_aq, y_py, y_aq, 1.0, mb(r),
                dict(mass_recovered_py=float(r.raw[:, 9].sum()) / float(r.raw[:, 8].sum())))
    res3[(n, outlet)] = y_py
fig_data["B3"] = (t_aq, y_aq, res3)

# ---------------------------------------------------------------------------------------------
# B4  saturated soil column with sorption (tutorial 5.3) -> analogue of the eggshell P-sorption bed
#     Q=0.001 m3/h, A=1e-3 m2, L=1 m, theta=0.4, rho_s=2300, no dispersion, 1 h pulse-ish input
# ---------------------------------------------------------------------------------------------
A_, L_, th, rho_s, Qin = 1e-3, 1.0, 0.4, 2300.0, 1e-3
Vpore = A_ * L_ * th; rho_th = rho_s * (1 - th) / th
dt = 0.0025; ns = 4
tm = (np.arange(int(2.0 / dt)) + 0.5) * dt
def pulse(h): return np.interp(tm, [0, 0.01, 0.5, 0.51, 10], [0, h, h, 0, 0])
soil = {}
for tag, krel in (("soil_a", 1e4), ("soil_b", 100.0)):
    _, curves = read_listing(AQ / f"{tag}__plot.txt")
    cases = {0: ("no sorption", 1.0, dict(p_kinetics="off")),
             1: ("linear", 1.0, dict(p_kinetics="relax_linear", K_d=0.00058, k_rel=krel)),
             2: ("Langmuir (pulse 1.0)", 1.0, dict(p_kinetics="relax", q_max=0.00029, K_p=0.5, k_rel=krel)),
             3: ("Langmuir (pulse 0.2)", 0.2, dict(p_kinetics="relax", q_max=0.00029, K_p=0.5, k_rel=krel))}
    for k, (lab, h, kw) in cases.items():
        t_aq, y_aq = curves[k][0][1:], curves[k][1][1:]
        for n, aqtag, ngrid in ((52, tag, "52 pts"), (502, tag + "_g502", "502 pts")):
            _, cv = read_listing(AQ / f"{aqtag}__plot.txt")
            t_aq, y_aq = cv[k][0][1:], cv[k][1][1:]
            u = Unit(name=f"{tag}_{lab}_n{n}", n_cells=n, volume=Vpore, scheme="tvd", rho_bulk_over_theta=rho_th, **kw)
            f4 = Forcing(np.full_like(tm, Qin), np.full_like(tm, 20.0), np.zeros_like(tm), pulse(h), dt=dt)
            r = simulate_unit(u, f4)
            y_py = r.series("cP_last")[ns - 1::ns]
            ft_py, ft_aq = front_time(t_aq, y_py, 0.5 * h), front_time(t_aq, y_aq, 0.5 * h)
            metrics(f"B4 {tag} sorption column", f"{lab}, python n={n} vs AQUASIM ({ngrid})", t_aq, y_py, y_aq, h, mb(r),
                    dict(t50_front_py=ft_py, t50_front_aquasim=ft_aq))
            soil[(tag, k, n)] = (t_aq, y_aq, y_py, lab)
fig_data["B4"] = soil

# analytical retardation check for the linear isotherm (arrival of the 50% front, pore-volume time tau = 0.4 h)
R = 1 + rho_th * 0.00058
rows.append(dict(case="B4 soil_a linear isotherm", quantity="analytic 50%%-front time tau*R (tau=0.4 h, R=%.3f)" % R,
                 t50_front_py=0.4 * R, t50_front_aquasim=np.nan))

# ---------------------------------------------------------------------------------------------
# B5  biofilm: quasi-steady diffusion-reaction (Monod) profile in AQUASIM's biofilm (tutorial 5.1 part A)
#     Given AQUASIM's own biomass profile X(z), thickness LF and bulk S, solve D S'' = (mu/Y) X S/(K+S)
#     (S'(0)=0 at substratum, S(LF)=S_surface) and compare with AQUASIM's pore-water profile S(z).
# ---------------------------------------------------------------------------------------------
_, lf = read_listing(AQ / "film_a__LF_t.txt"); tL, LF = lf[0]
_, sz = read_listing(AQ / "film_a__S_z.txt"); _, xz = read_listing(AQ / "film_a__X_z.txt")
par = {"S1": (0.4, 0.01, 5.0, 2e-5), "S2": (0.1, 0.5, 5.0, 2e-5)}
rho = 25000.0
film = {}
# curve order in S_z: for each of S1(10,25,50), S2(10,25,50): [pore water, bulk]
order = [("S1", 10), ("S1", 25), ("S1", 50), ("S2", 10), ("S2", 25), ("S2", 50)]
for j, (sp, day) in enumerate(order):
    z, S = sz[2 * j]; zb, Sb = sz[2 * j + 1]
    Sbulk = float(Sb[0]); Lf = float(np.interp(day, tL, LF))
    xi = {"S1": 0, "S2": 3}[sp] + {10: 0, 25: 1, 50: 2}[day]
    zx, eps = xz[xi]
    mu, Y, K_, D = par[sp]
    Xf = lambda zz: np.interp(zz, zx, eps) * rho
    # pore-water fraction  eps_w = 1 - eps_X1 - eps_X2  (AQUASIM eq. 3.16: flux = eps_w * D * dS/dz, rate per total volume)
    d_ = {10: 0, 25: 1, 50: 2}[day]
    z1, e1 = xz[d_]; z2, e2 = xz[3 + d_]
    ew = lambda zz: 1.0 - np.interp(zz, z1, e1) - np.interp(zz, z2, e2)
    def fun(zz, y):
        S_ = np.maximum(y[0], 0); w_ = ew(zz)
        dw = (ew(zz + 1e-9) - ew(zz - 1e-9)) / 2e-9                 # d(eps_w)/dz
        r = mu / Y * Xf(zz) * S_ / (K_ + S_)
        # y[1] = eps_w * D * dS/dz  (flux variable);  d(flux)/dz = r
        return np.vstack([y[1] / (w_ * D), r])
    def bc(ya, yb): return np.array([ya[1], yb[0] - S[-1]])
    zg = np.linspace(0, zx[-1], 400)
    sol = solve_bvp(fun, bc, zg, np.vstack([np.full_like(zg, S[-1] / 2), np.zeros_like(zg)]), tol=1e-8, max_nodes=100000)
    y_py = sol.sol(z)[0]
    metrics(f"B5 biofilm profile (film_a)", f"{sp} at day {day}: BVP vs AQUASIM pore-water profile", z, y_py, S, max(float(S.max()), 1e-9),
            extra=dict(bvp_converged=bool(sol.success)))
    film[(sp, day)] = (z, S, y_py)
fig_data["B5"] = film

# ---------------------------------------------------------------------------------------------
tab = pd.DataFrame(rows)
(HERE / "tables_AN_v1").mkdir(exist_ok=True)
tab.to_csv(HERE / "tables_AN_v1" / "aquasim_benchmark_errors_AN_v1.csv", index=False, float_format="%.6g")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(tab[["case", "quantity", "RMSE_pct_of_scale", "max_err_pct_of_scale"]].to_string())

# ------------------------------- figure --------------------------------------------------------
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
t, a, p, an = fig_data["B1"]
ax[0, 0].plot(t, a, "k-", lw=2, label="AQUASIM 2.1g"); ax[0, 0].plot(t, p, "r--", label="Python"); ax[0, 0].plot(t, an, "b:", label="analytic")
ax[0, 0].set_title("B1 batch first-order decay"); ax[0, 0].set_xlabel("t (h)"); ax[0, 0].set_ylabel("C (mg/L)"); ax[0, 0].legend()
t, a, p = fig_data["B2"]
ax[0, 1].plot(t, a, "k-", lw=2, label="AQUASIM"); ax[0, 1].plot(t, p, "r--", label="Python (MM kinetics)")
ax[0, 1].set_title("B2 batch Monod (Michaelis-Menten) decay"); ax[0, 1].set_xlabel("t (h)"); ax[0, 1].legend()
t, a, d = fig_data["B3"]
ax[0, 2].plot(t, a, "k-", lw=2, label="AQUASIM (outlet)")
for (n, ol), y in d.items(): ax[0, 2].plot(t, y, "--", label=f"Python n={n}, {ol} outlet")
ax[0, 2].set_title("B3 advective-dispersive tracer pulse (D=0.1)"); ax[0, 2].set_xlabel("t (arb.)"); ax[0, 2].legend()
for col, tag in ((0, "soil_a"), (1, "soil_b")):
    a_ = ax[1, col]
    for k, cl in zip(range(4), ("k", "tab:red", "tab:green", "tab:blue")):
        t, y_aq, y_py, lab = fig_data["B4"][(tag, k, 502)]
        a_.plot(t, y_aq, "-", color=cl, lw=2, label=f"AQUASIM {lab}")
        a_.plot(t, y_py, "--", color="orange" if k == 0 else cl, lw=1.2, alpha=.9)
    a_.set_title(f"B4 sorption column, k_rel={'1e4' if tag=='soil_a' else '100'}/h\nsolid=AQUASIM(502 pts), dashed=Python (TVD n=502)", fontsize=9)
    a_.set_xlabel("t (h)"); a_.set_ylabel("C (mg/m3)"); a_.legend(fontsize=7)
a_ = ax[1, 2]
for (sp, day), (z, S, y_py) in fig_data["B5"].items():
    a_.plot(z * 1e3, S, "-", lw=2, label=f"AQUASIM {sp} d{day}"); a_.plot(z * 1e3, y_py, "k--", lw=0.8)
a_.set_title("B5 biofilm quasi-steady S(z)\nAQUASIM (solid) vs diffusion-reaction BVP (dashed)", fontsize=9); a_.set_xlabel("z (mm)"); a_.set_ylabel("S (gCOD/m3)"); a_.legend(fontsize=6, ncol=2)
fig.suptitle("Consistency check of the Python equations vs AQUASIM 2.1g tutorial examples (not validation against field data)")
fig.tight_layout(); (HERE / "figures_AN_v1").mkdir(exist_ok=True)
fig.savefig(HERE / "figures_AN_v1" / "fig_aquasim_consistency_AN_v1.png", dpi=130)
