# Part of biowenp_py_AN_v1 (LGPL-3.0-or-later). Tests: rate laws vs Zhou et al. (2023) Table 4, CN assembly vs the
# MATLAB layout, Octave equivalence (if the stored run exists), mass-balance gates and analytic reactor solutions.
import os, sys, json, glob
import numpy as np
import pytest
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from biowenp_AN_v1 import HydroParams, PParams, HydroForcing, run_hydrology, run_p_faithful, MassBalanceError
from biowenp_AN_v1 import p_pools_AN_v1 as pp
from biowenp_AN_v1.transport_AN_v1 import cn_matrices, cn_step, dense_matrix
from biowenp_AN_v1.massbalance_AN_v1 import water_balance_full
from biowenp_AN_v1.synthetic_AN_v1 import synthetic_stormwater_10min, synthetic_reactor_daily
from biowenp_AN_v1.flowthrough_AN_v1 import FlowThroughUnit, ReactorForcing, simulate_flowthrough, simulate_train


# ---------------------------------------------------------------- rate laws (Table 4)
def test_langmuir_inverse_roundtrip():
    C = np.array([0.0, 0.01, 0.5, 3.0, 40.0])
    q = pp.langmuir_q(C, 683.5, 0.2)
    assert np.allclose(pp.langmuir_ceq(q, 683.5, 0.2), C, rtol=1e-12, atol=1e-14)
    assert np.all(q < 683.5)


def test_adsorption_rate_sign_and_zero_at_equilibrium():
    q = pp.langmuir_q(1.0, 683.5, 0.2)
    assert abs(pp.adsorption_rate(1.0, q, 0.15, 683.5, 0.2)) < 1e-12
    assert pp.adsorption_rate(2.0, q, 0.15, 683.5, 0.2) > 0 > pp.adsorption_rate(0.5, q, 0.15, 683.5, 0.2)


def test_barber_uptake_threshold_and_saturation():
    Cmin, Km, Imax = 0.02, 1.2, 3.0
    assert pp.barber_uptake(Cmin, Imax, Km, Cmin) == 0
    assert pp.barber_uptake(0.0, Imax, Km, Cmin) == 0
    assert np.isclose(pp.barber_uptake(Cmin + Km, Imax, Km, Cmin), Imax / 2)
    assert pp.barber_uptake(1e6, Imax, Km, Cmin) == pytest.approx(Imax, rel=1e-5)


def test_filtration_hydrolysis_precip_bedding():
    assert pp.filtration_rate(2.0, 3.0, 0.025) == pytest.approx(0.15)
    assert pp.hydrolysis_rate(10.0, 5e-4) == pytest.approx(5e-3)
    assert np.allclose(pp.precipitation_zero_order(np.array([0.0, 1e-5, 1.0]), 3.6e-4), [0.0, 1e-5, 3.6e-4])
    assert pp.bedding_filtration(1.0, 0.06, 72.9) == pytest.approx(np.exp(-0.06 * 72.9))


# ---------------------------------------------------------------- CN assembly vs the MATLAB statements
def _matlab_dense(C, C0t, C0t1, q, s, s1, alpha, p, Rn, eps=0.5):
    N = C.size
    alpha = np.broadcast_to(alpha, (N,)).astype(float); p = np.broadcast_to(p, (N,)).astype(float)
    Rn = np.broadcast_to(Rn, (N,)).astype(float)
    A = np.diag(alpha + 2 * eps * s - p) + np.diag((-eps * s - q / 4) * np.ones(N - 1), -1) + np.diag(
        (-eps * s + q / 4) * np.ones(N - 1), 1)
    B = np.zeros(N)
    B[1:N - 1] = ((1 - eps) * s + q / 4) * C[0:N - 2] + (alpha[1:N - 1] - 2 * (1 - eps) * s + p[1:N - 1]) * C[1:N - 1] + (
        (1 - eps) * s - q / 4) * C[2:N] + Rn[1:N - 1]
    A[0, 0:2] = [alpha[0] + 2 * eps * s1 - p[0], -eps * s1 + q / 4]
    B[0] = ((1 - eps) * s1 + q / 4) * (C0t + C0t1) + (alpha[0] - 2 * (1 - eps) * s1 + p[0]) * C[0] + (
        (1 - eps) * s1 - q / 4) * C[1] + Rn[0]
    A[N - 1, N - 2:N] = [-2 * eps * s, alpha[N - 1] + 2 * eps * s - p[N - 1]]
    B[N - 1] = (2 * (1 - eps) * s) * C[N - 2] + (alpha[N - 1] - 2 * (1 - eps) * s + p[N - 1]) * C[N - 1] + Rn[N - 1]
    return A, B


def test_cn_matches_matlab_layout():
    rng = np.random.default_rng(0)
    C = rng.random(45); alpha = 1 + rng.random(45); p = -0.01 * rng.random(45); Rn = 1e-3 * rng.random(45)
    ab, B = cn_matrices(C, 1.2, 1.3, 0.7, 0.4, 0.4, alpha, p, Rn)
    A_m, B_m = _matlab_dense(C, 1.2, 1.3, 0.7, 0.4, 0.4, alpha, p, Rn)
    assert np.allclose(dense_matrix(ab), A_m, atol=0) and np.allclose(B, B_m, rtol=1e-15)
    assert np.allclose(cn_step(C, 1.2, 1.3, 0.7, 0.4, 0.4, alpha, p, Rn), np.linalg.solve(A_m, B_m), rtol=1e-12)


def test_cn_converges_to_inflow_concentration():
    C = np.zeros(45)
    for _ in range(3000):
        C = cn_step(C, 2.0, 2.0, 0.5, 0.5, 0.5)
    assert np.allclose(C, 2.0, rtol=1e-6)


# ---------------------------------------------------------------- equivalence with the original (Octave)
def test_python_matches_octave_run_if_available():
    files = sorted(glob.glob(os.path.join(ROOT, "tables_AN_v1", "equivalence_octave_vs_python_*d_AN_v1.csv")))
    if not files:
        pytest.skip("no stored Octave equivalence run")
    import pandas as pd
    df = pd.read_csv(files[-1])
    assert np.nanmax(df.max_rel_diff.values) < 1e-8


# ---------------------------------------------------------------- hydrology water ledger
def test_hydrology_extended_water_ledger_closes():
    f = synthetic_stormwater_10min(n_days=30)
    h = HydroParams()
    hyd = run_hydrology(HydroForcing(f["precip"], f["temp"], f["month"]), h)
    wb = water_balance_full(hyd, h)
    assert abs(wb["residual_rel"]) < 1e-3


# ---------------------------------------------------------------- augmented flow-through model
def _const_forcing(n=200, Q=10.0, T=20.0, CN=10.0, CP=0.5, CPP=0.0):
    o = np.ones(n)
    return ReactorForcing(Q * o, T * o, CN * o, CP * o, CPP * o)


def test_cstr_first_order_analytic():
    u = FlowThroughUnit(n_cells=1, bed_volume=20, porosity=0.5, length=5, q_design=100, n_kinetics="first", k_N20=0.8)
    r = simulate_flowthrough(u, _const_forcing(), nsub=400)
    tau = 10.0 / 10.0
    assert r.daily.C_N_out.values[-1] == pytest.approx(10.0 / (1 + 0.8 * tau), rel=2e-3)


def test_tanks_in_series_first_order_analytic():
    u = FlowThroughUnit(n_cells=5, bed_volume=20, porosity=0.5, length=5, q_design=100, alpha_L=0.0,
                        n_kinetics="first", k_N20=0.8)
    r = simulate_flowthrough(u, _const_forcing(), nsub=800)
    tau = 1.0
    assert r.daily.C_N_out.values[-1] == pytest.approx(10.0 / (1 + 0.8 * tau / 5) ** 5, rel=5e-3)


def test_theta_temperature_factor():
    u = FlowThroughUnit(n_cells=1, bed_volume=20, porosity=0.5, length=5, q_design=100, n_kinetics="first",
                        k_N20=0.8, theta_N=1.08)
    r = simulate_flowthrough(u, _const_forcing(T=10.0), nsub=400)
    k = 0.8 * 1.08 ** (-10)
    assert r.daily.C_N_out.values[-1] == pytest.approx(10.0 / (1 + k), rel=2e-3)


def test_bypass_split_and_remix():
    u = FlowThroughUnit(n_cells=2, bed_volume=20, porosity=0.5, length=5, q_design=6.0, n_kinetics="first", k_N20=1.0)
    r = simulate_flowthrough(u, _const_forcing(Q=10.0))
    assert np.allclose(r.daily.Q_t, 6.0) and np.allclose(r.daily.Q_b, 4.0)
    assert r.ledger["V_b"] == pytest.approx(4.0 * 200)
    assert r.daily.C_N_out.values[-1] > 4.0 / 10.0 * 10.0    # bypassed water is untreated


@pytest.mark.parametrize("kin", ["zhou", "maria"])
def test_langmuir_capacity_depletes_to_equilibrium(kin):
    kw = dict(K_L=0.5, k_ad=20.0) if kin == "zhou" else dict(k_ads=10.0, k_des=20.0)
    u = FlowThroughUnit(n_cells=1, bed_volume=1, porosity=0.5, length=1, q_design=100, p_kinetics=kin,
                        sorbent_mass=10.0, q_max=5.0, **kw)
    r = simulate_flowthrough(u, _const_forcing(n=400, Q=5.0, CP=2.0))
    q_eq = 5.0 * 0.5 * 2.0 / (1 + 0.5 * 2.0)
    assert r.daily.q_mean.values[-1] == pytest.approx(q_eq, rel=1e-3)
    assert r.daily.C_P_out.values[-1] == pytest.approx(2.0, rel=1e-3)        # breakthrough after depletion
    assert np.all(np.diff(r.daily.q_mean.values) >= -1e-12)


def test_mass_gate_on_every_run_and_raises_when_broken():
    s = synthetic_reactor_daily(n_days=120)
    f = ReactorForcing(s["Q"], s["T"], s["C_N"], s["C_SRP"], s["C_PP"])
    wood = FlowThroughUnit(name="w", n_kinetics="mm", vmax_N20=12.0, K_N=2.0)
    egg = FlowThroughUnit(name="e", n_cells=3, bed_volume=3, porosity=0.45, length=3, p_kinetics="zhou",
                          sorbent_mass=150, q_max=35, K_L=0.5, k_ad=2.0, k_mi=0.01, lambda_pp=0.5)
    res, train = simulate_train([wood, egg], f)
    for r in res:
        for k, (a, b) in r.residuals().items():
            assert abs(a / b) < 1e-10
    assert abs(train["residual_N"] / train["N_in"]) < 1e-10
    res[1].ledger["P_filtered"] *= 1.01
    with pytest.raises(MassBalanceError):
        res[1].check_balance(1e-6)


def test_faithful_p_model_runs_and_reports_original_closure():
    f = synthetic_stormwater_10min(n_days=15)
    h, p = HydroParams(), PParams()
    hyd = run_hydrology(HydroForcing(f["precip"], f["temp"], f["month"]), h)
    po = run_p_faithful(hyd, h, p)
    assert np.all(np.isfinite(po["C_macro"])) and po["qe_ad"].max() < p.Qmax_ad
