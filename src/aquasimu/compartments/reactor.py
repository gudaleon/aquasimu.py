"""Shared reactor: one mixed cell, or many cells with advection and dispersion.

Units: time d, volume m3, concentration g/m3 (= mg/L), mass g (reported in kg).

A `Unit` is one reactor compartment:

* `n_cells = 1` is a completely mixed reactor (CSTR). `MixedReactor` fixes this.
* `n_cells > 1` is an advective-dispersive plug-flow compartment (finite volume,
  tanks in series, optional physical dispersion `disp_ratio = alpha_L / L`).
  `AdvDiffReactor` requires at least two cells.

High-flow bypass: inflow above `q_design` bypasses the unit and is mixed back
into the outlet in proportion to flow.

`simulate_unit` and `simulate_train` always run the mass-balance gate unless
`check=False` is passed explicitly for debugging.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from aquasimu.checks.massbalance import assert_unit_balance, component_balance, train_balance
from aquasimu.compartments import kernel as K
from aquasimu.variables.forcing import Forcing

__all__ = [
    "Forcing",
    "Result",
    "Unit",
    "choose_nsub",
    "simulate_train",
    "simulate_unit",
    "train_balance",
]


@dataclass
class Unit:
    name: str = "unit"
    n_cells: int = 1
    volume: float = 10.0            # pore-water volume (m3)
    q_design: float = 1e9           # design flow (m3/d); inflow above this bypasses
    disp_ratio: float = 0.0         # alpha_L / L  (dimensionless)
    scheme: str = "upwind"          # 'upwind' | 'central' | 'tvd' (van Leer flux limiter)
    outlet: str = "zero_gradient"   # 'zero_gradient' (Danckwerts) | 'transmission' (d2C/dx2 = 0; central only)
    n_kinetics: str = "off"         # 'off' | 'first_order' | 'mm'
    k_n20: float = 0.0              # 1/d
    vmax_n: float = 0.0             # g N/m3/d (per m3 pore water)
    K_n: float = 1.0                # g/m3
    theta_n: float = 1.0
    p_kinetics: str = "off"         # 'off' | 'first_order' | 'langmuir' | 'relax' | 'relax_linear'
    k_p20: float = 0.0              # 1/d
    k_ads: float = 0.0              # m3/g/d
    k_des: float = 0.0              # 1/d
    q_max: float = 0.0              # g P / kg media
    K_p: float = 1.0                # g/m3 (relax mode half-saturation)
    k_rel: float = 0.0              # 1/d (relax mode)
    theta_p: float = 1.0
    rho_bulk_over_theta: float = 0.0  # kg media per m3 pore water
    K_d: float = 0.0                # linear isotherm Kd (g/kg per g/m3), 'relax_linear' mode
    q0: float = 0.0                 # initial sorbed P (g/kg)
    cN0: float = 0.0                # initial pore-water NO3-N (g/m3)
    cP0: float = 0.0
    # optional biofilm (applies to NO3-N, first-order biofilm kinetics)
    bio_area: float = 0.0           # biofilm area per m3 pore water (m2/m3); 0 = off
    bio_kv20: float = 0.0           # intrinsic first-order rate in biofilm (1/d)
    bio_Df: float = 1e-4            # biofilm diffusivity (m2/d)
    bio_Lf: float = 1e-3            # biofilm thickness (m)
    bio_kL: float = 0.0             # external mass-transfer coefficient (m/d), 0 = no resistance
    bio_theta: float = 1.0
    t_ref: float = 20.0

    def __post_init__(self):
        if self.outlet == "transmission" and self.scheme != "central":
            raise ValueError("outlet='transmission' requires scheme='central'")

    def with_(self, **kwargs) -> Unit:
        return replace(self, **kwargs)

    def _pack(self):
        p = np.zeros(K.NP)
        p[K.P_V] = self.volume
        p[K.P_DISP] = self.disp_ratio
        p[K.P_QDES] = self.q_design
        p[K.P_KN20] = self.k_n20
        p[K.P_VMAXN] = self.vmax_n
        p[K.P_KN] = self.K_n
        p[K.P_THN] = self.theta_n
        p[K.P_KP20] = self.k_p20
        p[K.P_KADS] = self.k_ads
        p[K.P_KDES] = self.k_des
        p[K.P_QMAX] = self.q_max
        p[K.P_KP] = self.K_p
        p[K.P_KREL] = self.k_rel
        p[K.P_THP] = self.theta_p
        p[K.P_RHOB_TH] = self.rho_bulk_over_theta
        p[K.P_BIO_A] = self.bio_area
        p[K.P_BIO_KV] = self.bio_kv20
        p[K.P_BIO_DF] = self.bio_Df
        p[K.P_BIO_LF] = self.bio_Lf
        p[K.P_BIO_KL] = self.bio_kL
        p[K.P_THBIO] = self.bio_theta
        p[K.P_KD] = self.K_d
        p[K.P_OUTMODE] = 1.0 if self.outlet == "transmission" else 0.0
        p[K.P_Q0] = self.q0
        p[K.P_CN0] = self.cN0
        p[K.P_CP0] = self.cP0
        p[K.P_TREF] = self.t_ref
        return p

    def _modes(self):
        n_mode = {"off": 0, "first_order": 1, "mm": 2}[self.n_kinetics]
        p_mode = {"off": 0, "first_order": 1, "langmuir": 2, "relax": 3, "relax_linear": 4}[self.p_kinetics]
        scheme = {"upwind": 0, "central": 1, "tvd": 2}[self.scheme]
        return n_mode, p_mode, scheme

    def stiffness_rate(self, forcing: Forcing) -> float:
        """Upper estimate of the fastest rate (1/d) used to choose RK4 sub-steps."""
        n = self.n_cells
        cell_volume = self.volume / n
        qmax_flow = min(float(np.max(forcing.Q)), self.q_design)
        tmax = float(np.max(forcing.T))
        lam = qmax_flow * (2 + 4 * self.disp_ratio * n) / cell_volume
        if self.n_kinetics == "first_order":
            lam += self.k_n20 * max(self.theta_n ** (tmax - self.t_ref), 1.0)
        elif self.n_kinetics == "mm":
            lam += self.vmax_n * max(self.theta_n ** (tmax - self.t_ref), 1.0) / max(self.K_n, 1e-9)
        if self.bio_area > 0:
            kv = self.bio_kv20 * max(self.bio_theta ** (tmax - self.t_ref), 1.0)
            lam += self.bio_area * K._biofilm_flux_coeff(kv, self.bio_Df, self.bio_Lf, self.bio_kL)
        fp = max(self.theta_p ** (tmax - self.t_ref), 1.0)
        if self.p_kinetics == "first_order":
            lam += self.k_p20 * fp
        elif self.p_kinetics == "langmuir":
            cmax = float(max(np.max(forcing.cinP), self.cP0, 1e-9))
            lam += fp * (self.k_ads * (self.rho_bulk_over_theta * self.q_max + cmax) + self.k_des)
        elif self.p_kinetics == "relax":
            lam += fp * self.k_rel * (1 + self.rho_bulk_over_theta * self.q_max / max(self.K_p, 1e-9))
        elif self.p_kinetics == "relax_linear":
            lam += fp * self.k_rel * (1 + self.rho_bulk_over_theta * self.K_d)
        return lam


@dataclass
class Result:
    unit: Unit
    forcing: Forcing
    raw: np.ndarray
    nsub: int

    def series(self, name: str) -> np.ndarray:
        """State at the END of each interval.

        Names: 'cN_last', 'cP_last', 'cN_first', 'cP_first', 'cN_mean', 'cP_mean',
        'cN_last2', 'cP_last2'.
        """
        col = {
            "cN_last": K.R_CNL, "cP_last": K.R_CPL, "cN_first": K.R_CNF, "cP_first": K.R_CPF,
            "cN_mean": K.R_CNM, "cP_mean": K.R_CPM, "cN_last2": K.R_CNL2, "cP_last2": K.R_CPL2,
        }[name]
        return self.raw[:, col]

    @property
    def cout_N(self):
        """Flow-weighted outlet NO3-N, including bypass (g/m3)."""
        volume = self.raw[:, K.R_VIN]
        with np.errstate(invalid="ignore", divide="ignore"):
            concentration = (self.raw[:, K.R_NOUT] + self.raw[:, K.R_NBY]) / volume
        return np.where(volume > 0, concentration, self.forcing.cinN)

    @property
    def cout_P(self):
        """Flow-weighted outlet PO4-P, including bypass (g/m3)."""
        volume = self.raw[:, K.R_VIN]
        with np.errstate(invalid="ignore", divide="ignore"):
            concentration = (self.raw[:, K.R_POUT] + self.raw[:, K.R_PBY]) / volume
        return np.where(volume > 0, concentration, self.forcing.cinP)

    @property
    def cout_treated_N(self):
        volume = self.raw[:, K.R_VTR]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(volume > 0, self.raw[:, K.R_NOUT] / volume, np.nan)

    @property
    def cout_treated_P(self):
        volume = self.raw[:, K.R_VTR]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(volume > 0, self.raw[:, K.R_POUT] / volume, np.nan)

    def balance(self) -> dict:
        days = len(self.forcing) * self.forcing.dt
        return component_balance(self.raw, days)

    def check_balance(self, tol: float = 1e-6) -> dict:
        days = len(self.forcing) * self.forcing.dt
        return assert_unit_balance(self.raw, days, self.unit.name, tol)

    def balance_table(self, freq: str = "MS") -> pd.DataFrame:
        """Balance per calendar period (default month) in kg and m3, with closure residuals."""
        index = pd.date_range(
            self.forcing.start, periods=len(self.forcing), freq=f"{self.forcing.dt:g}D"
        )
        frame = pd.DataFrame(self.raw, index=index)
        grouped = frame.groupby(pd.Grouper(freq=freq)).sum()
        table = pd.DataFrame(index=grouped.index)
        table["water_in_m3"] = grouped[K.R_VIN]
        table["water_treated_out_m3"] = grouped[K.R_VTR]
        table["water_bypass_m3"] = grouped[K.R_VBY]
        table["water_residual_m3"] = table.water_in_m3 - table.water_treated_out_m3 - table.water_bypass_m3
        columns = {
            "N": (K.R_NIN, K.R_NOUT, K.R_NBY, K.R_NREM, K.R_NSTO),
            "P": (K.R_PIN, K.R_POUT, K.R_PBY, K.R_PREM, K.R_PSTO),
        }
        for symbol, (inflow, outflow, bypass, removed, storage) in columns.items():
            table[f"{symbol}_in_kg"] = grouped[inflow] / 1e3
            table[f"{symbol}_out_treated_kg"] = grouped[outflow] / 1e3
            table[f"{symbol}_bypass_kg"] = grouped[bypass] / 1e3
            table[f"{symbol}_removed_kg"] = grouped[removed] / 1e3
            table[f"{symbol}_storage_change_kg"] = grouped[storage] / 1e3
            table[f"{symbol}_residual_kg"] = (
                grouped[inflow] - grouped[outflow] - grouped[bypass] - grouped[removed] - grouped[storage]
            ) / 1e3
        table.insert(0, "unit", self.unit.name)
        return table


def choose_nsub(unit: Unit, forcing: Forcing, nsub_min: int = 4, safety: float = 2.0) -> int:
    lam = unit.stiffness_rate(forcing)
    return int(max(nsub_min, np.ceil(forcing.dt * lam / safety)))


def simulate_unit(
    unit: Unit,
    forcing: Forcing,
    nsub: int | None = None,
    tol: float = 1e-6,
    check: bool = True,
) -> Result:
    """Integrate one compartment. Mass balance is checked unless `check` is False."""
    n_mode, p_mode, scheme = unit._modes()
    steps = nsub or choose_nsub(unit, forcing)
    raw = K.simulate_kernel(
        forcing.Q, forcing.T, forcing.cinN, forcing.cinP, float(forcing.dt), int(steps),
        unit._pack(), int(unit.n_cells), scheme, n_mode, p_mode,
    )
    result = Result(unit, forcing, raw, steps)
    if check:
        result.check_balance(tol)
    return result


def simulate_train(
    units: list[Unit],
    forcing: Forcing,
    nsub: int | None = None,
    tol: float = 1e-6,
    check: bool = True,
) -> list[Result]:
    """Units in series. The mixed (treated + bypass) outlet of one unit feeds the next.

    Each unit is checked on its own. The train balance is checked at the end.
    """
    results = []
    current = forcing
    for unit in units:
        result = simulate_unit(unit, current, nsub=nsub, tol=tol, check=check)
        results.append(result)
        current = Forcing(current.Q, current.T, result.cout_N, result.cout_P, current.dt, current.start)
    if check:
        train_balance(results, tol)
    return results
