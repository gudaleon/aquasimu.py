"""Wrappers added around the ported reactor, and the readable rate functions."""
import numpy as np
import pytest

from aquasimu import AdvDiffReactor, Biofilm, Forcing, MixedReactor, simulate_unit
from aquasimu.processes.kinetics import langmuir_dqdt, nitrate_rate


def test_mixed_reactor_is_one_cell():
    unit = MixedReactor(n_cells=9, volume=10, n_kinetics="first_order", k_n20=0.2)
    assert unit.n_cells == 1
    forcing = Forcing(np.full(20, 5.0), np.full(20, 20.0), np.full(20, 4.0), np.full(20, 0.2))
    result = simulate_unit(unit, forcing)
    assert result.check_balance(1e-9)["N_rel_residual"] < 1e-9


def test_advdiff_needs_two_cells():
    with pytest.raises(ValueError):
        AdvDiffReactor(n_cells=1)


def test_biofilm_term_applies_and_closes():
    film = Biofilm(area=40.0, kv20=15.0, Df=1e-4, Lf=2e-4, kL=0.3, theta=1.05)
    unit = film.apply(AdvDiffReactor(n_cells=3, volume=12, n_kinetics="off"))
    assert unit.bio_area == 40.0
    forcing = Forcing(np.full(30, 8.0), np.full(30, 15.0), np.full(30, 6.0), np.zeros(30))
    result = simulate_unit(unit, forcing)
    balance = result.check_balance(1e-9)
    assert balance["N_removed_g"] > 0
    assert film.flux_coefficient(20.0) > 0


def test_documented_rates_match_the_equations():
    assert nitrate_rate(10.0, 20.0, mode="first_order", k20=1.0) == pytest.approx(10.0)
    assert nitrate_rate(10.0, 10.0, mode="first_order", k20=1.0, theta=1.1) == pytest.approx(10.0 * 1.1 ** -10)
    assert nitrate_rate(10.0, 20.0, mode="mm", vmax=8.0, half_saturation=2.0) == pytest.approx(8.0 * 10.0 / 12.0)
    assert langmuir_dqdt(1.0, 0.2, 20.0, k_ads=0.01, k_des=0.05, q_max=2.0) == pytest.approx(
        0.01 * 1.0 * (2.0 - 0.2) - 0.05 * 0.2
    )
