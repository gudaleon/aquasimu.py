"""Water, nitrogen and phosphorus closure.

Every simulation checks

    inflow = treated outflow + bypass + removed + change in storage

for water (no removal or storage change), NO3-N and PO4-P. The relative
residual is the absolute residual divided by the largest term of that balance.
`simulate_unit` and `simulate_train` raise `MassBalanceError` when any relative
residual exceeds the tolerance (default 1e-6). That check is the hard gate:
calibration forward runs use it as well, so a run that does not close is a
failed run.

The cumulative outflow and process sinks are extra states of the same ODE
system as the concentrations. A Runge-Kutta step conserves linear invariants,
so a closed balance is expected to round-off and a large residual means a
flux was dropped or a step blew up. Closure does not show that the physics
are the right physics; the analytical and independent-solver tests exist for
that.
"""
from __future__ import annotations

import numpy as np

from aquasimu.compartments import kernel as K
from aquasimu.core.errors import MassBalanceError


def _scale(*terms: float) -> float:
    return max(max(abs(x) for x in terms), 1e-300)


def component_balance(raw: np.ndarray, period_days: float) -> dict:
    """Integrate one unit's kernel output into a water / N / P balance."""
    summed = lambda col: float(np.sum(raw[:, col]))
    water_in = summed(K.R_VIN)
    water_treated = summed(K.R_VTR)
    water_bypass = summed(K.R_VBY)
    out = dict(
        water_in_m3=water_in,
        water_treated_out_m3=water_treated,
        water_bypass_m3=water_bypass,
        water_storage_change_m3=0.0,
        water_residual_m3=water_in - water_treated - water_bypass,
        N_in_g=summed(K.R_NIN),
        N_out_treated_g=summed(K.R_NOUT),
        N_bypass_g=summed(K.R_NBY),
        N_removed_g=summed(K.R_NREM),
        N_storage_change_g=summed(K.R_NSTO),
        P_in_g=summed(K.R_PIN),
        P_out_treated_g=summed(K.R_POUT),
        P_bypass_g=summed(K.R_PBY),
        P_removed_g=summed(K.R_PREM),
        P_sorbed_g=summed(K.R_PSORB),
        P_storage_change_g=summed(K.R_PSTO),
    )
    out["N_residual_g"] = (
        out["N_in_g"] - out["N_out_treated_g"] - out["N_bypass_g"]
        - out["N_removed_g"] - out["N_storage_change_g"]
    )
    out["P_residual_g"] = (
        out["P_in_g"] - out["P_out_treated_g"] - out["P_bypass_g"]
        - out["P_removed_g"] - out["P_storage_change_g"]
    )
    out["water_rel_residual"] = abs(out["water_residual_m3"]) / _scale(
        water_in, water_treated, water_bypass
    )
    out["N_rel_residual"] = abs(out["N_residual_g"]) / _scale(
        out["N_in_g"], out["N_out_treated_g"], out["N_bypass_g"],
        out["N_removed_g"], out["N_storage_change_g"],
    )
    out["P_rel_residual"] = abs(out["P_residual_g"]) / _scale(
        out["P_in_g"], out["P_out_treated_g"], out["P_bypass_g"],
        out["P_removed_g"], out["P_storage_change_g"],
    )
    out["P_sorbed_end_g"] = float(raw[-1, K.R_PQMASS])
    out["period_days"] = period_days
    per_year = 365.0 / period_days / 1000.0
    out["kgN_retained_per_yr"] = out["N_removed_g"] * per_year
    out["kgP_retained_per_yr"] = out["P_removed_g"] * per_year
    out["kgN_in_per_yr"] = out["N_in_g"] * per_year
    out["kgP_in_per_yr"] = out["P_in_g"] * per_year
    out["kgN_bypass_per_yr"] = out["N_bypass_g"] * per_year
    out["kgP_bypass_per_yr"] = out["P_bypass_g"] * per_year
    return out


def assert_unit_balance(raw: np.ndarray, period_days: float, name: str, tol: float = 1e-6) -> dict:
    """Return the balance, or raise `MassBalanceError` if a relative residual exceeds `tol`."""
    balance = component_balance(raw, period_days)
    bad = [key for key in ("water_rel_residual", "N_rel_residual", "P_rel_residual") if balance[key] > tol]
    if bad:
        detail = ", ".join(f"{key}={balance[key]:.3e}" for key in bad)
        raise MassBalanceError(f"[{name}] mass balance failed (tol={tol:g}): {detail}")
    return balance


def train_balance(results: list, tol: float = 1e-6) -> dict:
    """Overall balance of units in series (inflow to the first vs outlet of the last plus every sink).

    Raises `MassBalanceError` when water, N or P fails to close.
    """
    first = results[0].balance()
    last = results[-1]
    forcing = results[0].forcing
    water_in = float(np.sum(forcing.Q) * forcing.dt)
    out_n = float(np.sum(last.raw[:, K.R_NOUT] + last.raw[:, K.R_NBY]))
    out_p = float(np.sum(last.raw[:, K.R_POUT] + last.raw[:, K.R_PBY]))
    removed_n = sum(result.balance()["N_removed_g"] for result in results)
    removed_p = sum(result.balance()["P_removed_g"] for result in results)
    storage_n = sum(result.balance()["N_storage_change_g"] for result in results)
    storage_p = sum(result.balance()["P_storage_change_g"] for result in results)
    out = dict(
        water_in_m3=water_in,
        water_out_m3=float(np.sum(last.raw[:, K.R_VIN])),
        N_in_g=first["N_in_g"],
        N_out_total_g=out_n,
        N_removed_g=removed_n,
        N_storage_change_g=storage_n,
        P_in_g=first["P_in_g"],
        P_out_total_g=out_p,
        P_removed_g=removed_p,
        P_storage_change_g=storage_p,
    )
    out["N_residual_g"] = first["N_in_g"] - out_n - removed_n - storage_n
    out["P_residual_g"] = first["P_in_g"] - out_p - removed_p - storage_p
    # Denominator is the inflow. A zero inflow (a batch) uses a tiny floor so the
    # gate still runs; with a nonzero inflow this is the original ratio.
    out["N_rel_residual"] = abs(out["N_residual_g"]) / max(abs(first["N_in_g"]), 1e-300)
    out["P_rel_residual"] = abs(out["P_residual_g"]) / max(abs(first["P_in_g"]), 1e-300)
    out["water_rel_residual"] = abs(water_in - out["water_out_m3"]) / max(abs(water_in), 1e-300)
    bad = [key for key in ("water_rel_residual", "N_rel_residual", "P_rel_residual") if out[key] > tol]
    if bad:
        detail = ", ".join(f"{key}={out[key]:.3e}" for key in bad)
        raise MassBalanceError(f"train mass balance failed: {detail}")
    per_year = 365.0 / (len(forcing) * forcing.dt) / 1000.0
    out["kgN_retained_per_yr"] = removed_n * per_year
    out["kgP_retained_per_yr"] = removed_p * per_year
    return out
