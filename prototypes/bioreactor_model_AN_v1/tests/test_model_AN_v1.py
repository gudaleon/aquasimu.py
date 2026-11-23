"""Unit tests (AN_v1).  Run:  python -m pytest -q tests/"""
import numpy as np, pytest
from scipy.integrate import solve_ivp, solve_bvp
from bioreactor_AN_v1 import (Unit, Forcing, simulate_unit, simulate_train, train_balance,
                              MassBalanceError)
from bioreactor_AN_v1 import kernel_AN_v1 as K


def const_forcing(days=400, Q=20.0, T=20.0, cN=10.0, cP=0.5, dt=1.0):
    return Forcing(np.full(days, Q), np.full(days, T), np.full(days, cN), np.full(days, cP), dt)


def test_mass_balance_closes_all_processes_with_bypass():
    rng = np.random.default_rng(1)
    d = 365
    Q = 15 + 10 * rng.random(d); Q[::17] = 150.0                 # storms -> bypass
    f = Forcing(Q, 12 + 8 * np.sin(np.arange(d) / 58), 8 + 4 * rng.random(d), 0.3 + 0.4 * rng.random(d))
    u = Unit(n_cells=5, volume=30, q_design=40, disp_ratio=0.05, n_kinetics="mm", vmax_n=6, K_n=2,
             theta_n=1.07, p_kinetics="langmuir", k_ads=0.005, k_des=0.01, q_max=1.5, theta_p=1.02,
             rho_bulk_over_theta=2400, bio_area=50, bio_kv20=20, bio_kL=0.5, bio_Df=1e-4, bio_Lf=2e-4)
    r = simulate_unit(u, f)
    b = r.check_balance(1e-9)
    assert b["water_bypass_m3"] > 0
    assert b["N_rel_residual"] < 1e-9 and b["P_rel_residual"] < 1e-9 and b["water_rel_residual"] < 1e-12
    assert b["kgN_retained_per_yr"] > 0 and b["kgP_retained_per_yr"] > 0
    tbl = r.balance_table()
    assert np.allclose(tbl[["N_residual_kg", "P_residual_kg", "water_residual_m3"]].to_numpy(), 0, atol=1e-8)


def test_balance_fails_loudly():
    r = simulate_unit(Unit(volume=10, n_kinetics="first_order", k_n20=1.0), const_forcing(30))
    r.raw[5, K.R_NOUT] *= 1.5                                      # corrupt a flux
    with pytest.raises(MassBalanceError):
        r.check_balance(1e-6)


def test_cstr_first_order_steady_state_matches_analytical():
    V, Q, k, c0 = 20.0, 20.0, 1.0, 10.0
    r = simulate_unit(Unit(n_cells=1, volume=V, n_kinetics="first_order", k_n20=k), const_forcing(Q=Q, cN=c0))
    tau = V / Q
    assert r.cout_N[-1] == pytest.approx(c0 / (1 + k * tau), rel=1e-8)
    assert r.balance()["kgN_retained_per_yr"] == pytest.approx(Q * (c0 - c0 / (1 + k * tau)) * 365 / 1000, rel=2e-3)


def test_cstr_michaelis_menten_steady_state_matches_analytical():
    V, Q, vmax, Kn, c0 = 20.0, 20.0, 8.0, 2.0, 10.0
    r = simulate_unit(Unit(n_cells=1, volume=V, n_kinetics="mm", vmax_n=vmax, K_n=Kn), const_forcing(Q=Q, cN=c0))
    # steady state: Q(c0 - c) = V vmax c/(K+c) -> c^2 + (K + V vmax/Q - c0) c - c0 K = 0
    b = Kn + V * vmax / Q - c0
    c = (-b + np.sqrt(b * b + 4 * c0 * Kn)) / 2
    assert r.cout_N[-1] == pytest.approx(c, rel=1e-8)


def test_tanks_in_series_steady_state_matches_analytical():
    V, Q, k, c0, n = 30.0, 20.0, 0.8, 10.0, 6
    r = simulate_unit(Unit(n_cells=n, volume=V, n_kinetics="first_order", k_n20=k), const_forcing(Q=Q, cN=c0, days=200))
    assert r.cout_N[-1] == pytest.approx(c0 / (1 + k * V / Q / n) ** n, rel=1e-8)


def test_arrhenius_temperature_dependence():
    V, Q, k, th = 20.0, 20.0, 1.0, 1.1
    for T in (5.0, 12.0, 25.0):
        r = simulate_unit(Unit(volume=V, n_kinetics="first_order", k_n20=k, theta_n=th), const_forcing(Q=Q, T=T))
        kT = k * th ** (T - 20)
        assert r.cout_N[-1] == pytest.approx(10 / (1 + kT * V / Q), rel=1e-8)


def test_langmuir_sorption_depletes_and_conserves():
    # batch (Q=0): equilibrium of the kinetic Langmuir model must satisfy mass conservation and the isotherm
    V, rho_th, kads, kdes, qmax, c0 = 10.0, 1000.0, 0.01, 0.05, 2.0, 1.0
    f = const_forcing(days=4000, Q=0.0, cP=0.0)
    u = Unit(volume=V, p_kinetics="langmuir", k_ads=kads, k_des=kdes, q_max=qmax, rho_bulk_over_theta=rho_th, cP0=c0)
    r = simulate_unit(u, f, nsub=40)
    mass_total = c0 * V
    sorbed = r.raw[-1, K.R_PQMASS]
    dissolved = mass_total - sorbed
    ce = dissolved / V
    q_eq = sorbed / (rho_th * V)
    assert q_eq == pytest.approx(kads * ce * qmax / (kads * ce + kdes), rel=1e-4)   # Langmuir isotherm
    assert q_eq < qmax                                                              # finite capacity
    # a flow-through column saturates: removal efficiency falls with time as capacity depletes
    f2 = const_forcing(days=300, Q=10.0, cP=1.0)
    u2 = u.with_(cP0=0.0, q_max=0.02, k_des=0.0)
    r2 = simulate_unit(u2, f2, nsub=40)
    assert r2.cout_P[-1] > r2.cout_P[1] + 0.2
    assert r2.raw[-1, K.R_PQMASS] <= u2.q_max * rho_th * V * (1 + 1e-6)


def test_matches_independent_scipy_solution_advective_dispersive():
    n, V, Q = 8, 12.0, 6.0
    u = Unit(n_cells=n, volume=V, disp_ratio=0.2, scheme="central", n_kinetics="mm", vmax_n=5, K_n=3, theta_n=1.0,
             p_kinetics="relax", k_rel=2.0, q_max=0.05, K_p=0.4, rho_bulk_over_theta=300)
    f = const_forcing(days=30, Q=Q, T=20, cN=9, cP=0.6)
    r = simulate_unit(u, f)
    Vc = V / n; E = 0.2 * n * Q

    def rhs(t, y):
        cN, cP, q = y[:n], y[n:2 * n], y[2 * n:3 * n]
        dN = np.zeros(n); dP = np.zeros(n)
        for c, cin, d in ((cN, 9, dN), (cP, 0.6, dP)):
            for i in range(n):
                fin = Q * cin if i == 0 else Q * 0.5 * (c[i - 1] + c[i]) - E * (c[i] - c[i - 1])
                fout = Q * c[i] if i == n - 1 else Q * 0.5 * (c[i] + c[i + 1]) - E * (c[i + 1] - c[i])
                d[i] = (fin - fout) / Vc
        dq = 2.0 * (0.05 * cP / (0.4 + cP) - q)
        return np.r_[dN - 5 * cN / (3 + cN), dP - 300 * dq, dq, Q * cN[-1], Q * cP[-1]]
    sol = solve_ivp(rhs, (0, 30), np.zeros(3 * n + 2), method="LSODA", rtol=1e-10, atol=1e-12)
    # compare cumulative treated-outflow load over the 30 d (independent integration of the same equations)
    outN_ref, outP_ref = sol.y[3 * n, -1], sol.y[3 * n + 1, -1]
    assert r.raw[:, K.R_NOUT].sum() == pytest.approx(outN_ref, rel=1e-6)
    assert r.raw[:, K.R_POUT].sum() == pytest.approx(outP_ref, rel=1e-6)
    assert 1.0 < outP_ref < 0.6 * 6 * 30


def test_biofilm_effectiveness_matches_numerical_diffusion_reaction():
    kv, Df, Lf, kL, C = 30.0, 1e-4, 4e-4, 0.4, 5.0                # 1/d, m2/d, m, m/d, g/m3
    coef = K._biofilm_flux_coeff(kv, Df, Lf, kL)                  # J/C
    # solve Df C'' = kv C, C'(0)=0 (substratum), Df C'(L) = kL (Cb - C(L)) numerically (bulk C fixed)
    def f(x, y): return np.vstack([y[1], kv / Df * y[0]])
    def bc(ya, yb): return np.array([ya[1], Df * yb[1] - kL * (C - yb[0])])
    x = np.linspace(0, Lf, 200)
    sol = solve_bvp(f, bc, x, np.vstack([np.full_like(x, C / 2), np.zeros_like(x)]), tol=1e-10)
    J_num = Df * sol.y[1, -1]
    assert J_num == pytest.approx(coef * C, rel=1e-5)


def test_bypass_water_and_load_accounting():
    f = Forcing(np.array([10.0, 100.0, 10.0]), np.full(3, 20.0), np.full(3, 10.0), np.full(3, 1.0))
    r = simulate_unit(Unit(volume=10, q_design=30, n_kinetics="first_order", k_n20=1.0), f)
    b = r.balance()
    assert b["water_bypass_m3"] == pytest.approx(70.0)
    assert b["N_bypass_g"] == pytest.approx(70.0 * 10.0)


def test_train_balance():
    f = const_forcing(days=200, Q=25, cN=9, cP=0.8)
    wood = Unit(name="woodchip", n_cells=4, volume=30, q_design=40, n_kinetics="mm", vmax_n=6, K_n=2, theta_n=1.07)
    egg = Unit(name="eggshell", n_cells=2, volume=2, q_design=40, p_kinetics="langmuir", k_ads=0.005, k_des=0.01,
               q_max=1.5, rho_bulk_over_theta=2400)
    rs = simulate_train([wood, egg], f)
    tb = train_balance(rs, 1e-9)
    assert tb["kgN_retained_per_yr"] > 0 and tb["kgP_retained_per_yr"] > 0


def test_tvd_scheme_conserves_mass_and_is_less_diffusive_than_upwind():
    """A square N pulse through a plug-flow column: TVD stays non-negative/bounded, closes the balance, and keeps a
    sharper front than first-order upwind (smaller numerical diffusion)."""
    days = 60
    cN = np.zeros(days); cN[5:8] = 10.0          # 3-day pulse, residence time 5 d (daily steps; peak stays below 10)
    f = Forcing(np.full(days, 10.0), np.full(days, 20.0), cN, np.zeros(days))
    # pore volume 50 m3, Q=10 -> tau=5 d; dispersion-free column, so a sharper scheme keeps a higher peak
    out = {}
    for sch in ("upwind", "tvd"):
        r = simulate_unit(Unit(n_cells=100, volume=50, scheme=sch), f)
        r.check_balance(1e-9)
        out[sch] = r.cout_N
        assert out[sch].min() > -1e-9 and out[sch].max() < 10.0 + 1e-6      # no under/overshoot
    # the sharper front reaches a higher peak concentration
    assert out["tvd"].max() > out["upwind"].max()


def test_transmission_outlet_requires_central_and_closes_balance():
    with pytest.raises(ValueError):
        Unit(scheme="upwind", outlet="transmission")
    f = const_forcing(days=40)
    r = simulate_unit(Unit(n_cells=40, volume=20, disp_ratio=0.1, scheme="central", outlet="transmission",
                           n_kinetics="first_order", k_n20=0.3), f)
    b = r.check_balance(1e-9)
    assert b["N_rel_residual"] < 1e-9
