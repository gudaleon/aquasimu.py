"""Consistency check against user-supplied AQUASIM listing files.

This compares the reactor equations with another implementation of published
AQUASIM tutorial cases. It is a consistency check of the numerics, not a
validation against measurements, and not a soil, river or lake compartment.

No listing is shipped. Set the directory yourself (see ``LISTING_ENV``). The
expected files are the batch-program listings from those tutorial runs:

    proc_a__Conc.txt
    procB_mm_mod__Conc.txt
    advdif_a__C_t.txt
    soil_a__plot.txt, soil_a_g502__plot.txt
    soil_b__plot.txt, soil_b_g502__plot.txt
    film_a__LF_t.txt, film_a__S_z.txt, film_a__X_z.txt

The soil-column cases are run with the advective-diffusive reactor as a
sorption analogue. They do not add a soil compartment.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_bvp

from aquasimu import Forcing, Unit, simulate_unit
from aquasimu.compartments import kernel as K
from aquasimu.io.listing import read_listing

LISTING_ENV = "AQUASIMU_AQUASIM_LISTINGS"


def listing_dir_from_env() -> Path | None:
    value = os.environ.get(LISTING_ENV)
    if not value:
        return None
    return Path(value)


def _metrics(rows, case, what, time, y_py, y_aq, scale, balance=None, extra=None):
    err = y_py - y_aq
    record = dict(
        case=case,
        quantity=what,
        n_points=len(time),
        RMSE=float(np.sqrt(np.mean(err ** 2))),
        max_abs_err=float(np.max(np.abs(err))),
        scale=scale,
        RMSE_pct_of_scale=100 * float(np.sqrt(np.mean(err ** 2))) / scale,
        max_err_pct_of_scale=100 * float(np.max(np.abs(err))) / scale,
    )
    if balance:
        record.update(balance)
    if extra:
        record.update(extra)
    rows.append(record)
    return record


def _front_time(time, values, level):
    index = int(np.argmax(values >= level))
    if values[index] < level:
        return np.nan
    return float(
        time[index - 1]
        + (level - values[index - 1]) * (time[index] - time[index - 1]) / (values[index] - values[index - 1])
    )


def _balance(result):
    balance = result.check_balance(1e-6)
    return dict(
        N_rel_resid=balance["N_rel_residual"],
        P_rel_resid=balance["P_rel_residual"],
        water_rel_resid=balance["water_rel_residual"],
    )


def run_consistency(listing_dir: Path | str) -> pd.DataFrame:
    """Run the tutorial comparisons and return an error table. Raises if a balance fails."""
    listing_dir = Path(listing_dir)
    rows: list[dict] = []

    # B1 batch, first-order decay. k = 1/h, C0 = 1, V = 10, 10 h at dt = 0.1.
    _, curves = read_listing(listing_dir / "proc_a__Conc.txt")
    time_aq, aquasim = curves[0]
    time_aq, aquasim = time_aq[1:], aquasim[1:]
    forcing = Forcing(np.zeros(100), np.full(100, 20.0), np.zeros(100), np.zeros(100), dt=0.1)
    result = simulate_unit(
        Unit(name="batch_first_order", volume=10, n_kinetics="first_order", k_n20=1.0, cN0=1.0),
        forcing,
    )
    python = result.series("cN_mean")
    _metrics(rows, "B1 batch first-order (proc_a)", "C_A(t) python vs AQUASIM", time_aq, python, aquasim, 1.0, _balance(result))
    _metrics(
        rows, "B1 batch first-order (proc_a)", "C_A(t) AQUASIM vs analytic exp(-kt)",
        time_aq, aquasim, np.exp(-time_aq), 1.0,
    )

    # B2 batch Michaelis-Menten. The listing is a modified initial condition (C_B(0) = 1).
    _, curves = read_listing(listing_dir / "procB_mm_mod__Conc.txt")
    time_aq, aquasim = curves[1][0][1:], curves[1][1][1:]
    result = simulate_unit(
        Unit(name="batch_mm", volume=10, n_kinetics="mm", vmax_n=0.25, K_n=0.5, cN0=1.0),
        forcing,
    )
    _metrics(
        rows, "B2 batch Monod (proc_b, C_B(0)=1)", "C_B(t) python vs AQUASIM",
        time_aq, result.series("cN_mean"), aquasim, 1.0, _balance(result),
    )

    # B3 advective-diffusive tracer. Q = 1, A = 1, L = 1, D = 0.1. The pulse is carried as PO4-P.
    _, curves = read_listing(listing_dir / "advdif_a__C_t.txt")
    time_aq, aquasim = curves[1]
    time_aq, aquasim = time_aq[1:], aquasim[1:]
    dt, stride = 0.0005, 20
    tm = (np.arange(int(3.0 / dt)) + 0.5) * dt
    cin = np.interp(tm, [0, 0.01, 1.0, 1.01, 10], [0, 1, 1, 0, 0])
    pulse = Forcing(np.ones_like(tm), np.full_like(tm, 20.0), np.zeros_like(tm), cin, dt=dt)
    for n_cells, outlet in ((50, "transmission"), (100, "transmission"), (200, "transmission"), (200, "zero_gradient")):
        unit = Unit(
            name=f"AD_column_n{n_cells}_{outlet}",
            n_cells=n_cells, volume=1.0, disp_ratio=0.1, scheme="central", outlet=outlet,
        )
        result = simulate_unit(unit, pulse)
        if outlet == "transmission":
            face = 1.5 * result.series("cP_last") - 0.5 * result.series("cP_last2")
        else:
            face = result.series("cP_last")
        recovered = float(result.raw[:, K.R_POUT].sum()) / float(result.raw[:, K.R_PIN].sum())
        _metrics(
            rows,
            f"B3 AD tracer pulse (advdif_a) n={n_cells} outlet={outlet}",
            "C(x=L,t) python vs AQUASIM",
            time_aq, face[stride - 1::stride], aquasim, 1.0, _balance(result),
            dict(mass_recovered_py=recovered),
        )

    # B4 sorption columns, integrated with the advective-diffusive reactor (not a soil compartment).
    area, length, theta, rho_s, q_in = 1e-3, 1.0, 0.4, 2300.0, 1e-3
    pore_volume = area * length * theta
    rho_over_theta = rho_s * (1 - theta) / theta
    dt, stride = 0.0025, 4
    tm = (np.arange(int(2.0 / dt)) + 0.5) * dt

    def pulse_height(height):
        return np.interp(tm, [0, 0.01, 0.5, 0.51, 10], [0, height, height, 0, 0])

    for tag, k_rel in (("soil_a", 1e4), ("soil_b", 100.0)):
        cases = {
            0: ("no sorption", 1.0, dict(p_kinetics="off")),
            1: ("linear", 1.0, dict(p_kinetics="relax_linear", K_d=0.00058, k_rel=k_rel)),
            2: ("Langmuir (pulse 1.0)", 1.0, dict(p_kinetics="relax", q_max=0.00029, K_p=0.5, k_rel=k_rel)),
            3: ("Langmuir (pulse 0.2)", 0.2, dict(p_kinetics="relax", q_max=0.00029, K_p=0.5, k_rel=k_rel)),
        }
        for index, (label, height, kwargs) in cases.items():
            for n_cells, file_tag, grid in ((52, tag, "52 pts"), (502, tag + "_g502", "502 pts")):
                _, curves = read_listing(listing_dir / f"{file_tag}__plot.txt")
                time_aq, aquasim = curves[index][0][1:], curves[index][1][1:]
                unit = Unit(
                    name=f"{tag}_{label}_n{n_cells}",
                    n_cells=n_cells, volume=pore_volume, scheme="tvd",
                    rho_bulk_over_theta=rho_over_theta, **kwargs,
                )
                column_forcing = Forcing(
                    np.full_like(tm, q_in), np.full_like(tm, 20.0), np.zeros_like(tm), pulse_height(height), dt=dt,
                )
                result = simulate_unit(unit, column_forcing)
                python = result.series("cP_last")[stride - 1::stride]
                _metrics(
                    rows, f"B4 {tag} sorption column",
                    f"{label}, python n={n_cells} vs AQUASIM ({grid})",
                    time_aq, python, aquasim, height, _balance(result),
                    dict(
                        t50_front_py=_front_time(time_aq, python, 0.5 * height),
                        t50_front_aquasim=_front_time(time_aq, aquasim, 0.5 * height),
                    ),
                )

    retardation = 1 + rho_over_theta * 0.00058
    rows.append(dict(
        case="B4 soil_a linear isotherm",
        quantity=f"analytic 50%-front time tau*R (tau=0.4 h, R={retardation:.3f})",
        t50_front_py=0.4 * retardation,
        t50_front_aquasim=np.nan,
    ))

    # B5 quasi-steady diffusion-reaction profile, given AQUASIM's own biomass profile.
    # The package kernel's biofilm term is first-order only; this check is the BVP.
    _, thickness_curves = read_listing(listing_dir / "film_a__LF_t.txt")
    time_thickness, thickness = thickness_curves[0]
    _, substrate = read_listing(listing_dir / "film_a__S_z.txt")
    _, biomass = read_listing(listing_dir / "film_a__X_z.txt")
    parameters = {"S1": (0.4, 0.01, 5.0, 2e-5), "S2": (0.1, 0.5, 5.0, 2e-5)}
    rho = 25000.0
    order = [("S1", 10), ("S1", 25), ("S1", 50), ("S2", 10), ("S2", 25), ("S2", 50)]
    for j, (species, day) in enumerate(order):
        z, profile = substrate[2 * j]
        _zb, bulk = substrate[2 * j + 1]
        film_thickness = float(np.interp(day, time_thickness, thickness))
        xi = {"S1": 0, "S2": 3}[species] + {10: 0, 25: 1, 50: 2}[day]
        zx, eps = biomass[xi]
        mu, yield_coef, half_sat, diffusivity = parameters[species]
        biomass_density = lambda zz, zx=zx, eps=eps: np.interp(zz, zx, eps) * rho
        day_index = {10: 0, 25: 1, 50: 2}[day]
        z1, e1 = biomass[day_index]
        z2, e2 = biomass[3 + day_index]
        water_fraction = lambda zz, z1=z1, e1=e1, z2=z2, e2=e2: (
            1.0 - np.interp(zz, z1, e1) - np.interp(zz, z2, e2)
        )

        def ode(zz, y, mu=mu, yield_coef=yield_coef, half_sat=half_sat, diffusivity=diffusivity):
            substrate_c = np.maximum(y[0], 0)
            fraction = water_fraction(zz)
            rate = mu / yield_coef * biomass_density(zz) * substrate_c / (half_sat + substrate_c)
            return np.vstack([y[1] / (fraction * diffusivity), rate])

        surface = profile[-1]

        def boundary(ya, yb, surface=surface):
            return np.array([ya[1], yb[0] - surface])

        grid = np.linspace(0, zx[-1], 400)
        solution = solve_bvp(
            ode, boundary, grid,
            np.vstack([np.full_like(grid, profile[-1] / 2), np.zeros_like(grid)]),
            tol=1e-8, max_nodes=100000,
        )
        _metrics(
            rows, "B5 biofilm profile (film_a)",
            f"{species} at day {day}: BVP vs AQUASIM pore-water profile",
            z, solution.sol(z)[0], profile, max(float(profile.max()), 1e-9),
            extra=dict(bvp_converged=bool(solution.success), Lf_m=film_thickness),
        )

    return pd.DataFrame(rows)
