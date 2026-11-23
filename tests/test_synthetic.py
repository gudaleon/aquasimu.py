"""The shipped calibration case is synthetic and closes its mass balance.

The retained-mass numbers are the seeded case this port was taken from.
They are not field measurements.
"""
import pytest

from aquasimu.examples.bioreactor.synthetic import TRUE, make_forcing, run_model


def test_synthetic_truth_closes_and_matches_the_seeded_case():
    out = run_model(TRUE, make_forcing())
    assert out["kgN"] == pytest.approx(25.9129, rel=1e-4)
    assert out["kgP"] == pytest.approx(1.57006, rel=1e-4)
    balance = out["train_balance"]
    assert balance["N_rel_residual"] < 1e-9
    assert balance["P_rel_residual"] < 1e-9
    assert balance["water_rel_residual"] < 1e-12
