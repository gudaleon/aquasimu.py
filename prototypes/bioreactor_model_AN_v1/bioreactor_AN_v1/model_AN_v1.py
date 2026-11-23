"""AQUASIM-style bioreactor model, Python port (AN_v1).

Units: time d, volume m3, concentration g/m3 (= mg/L), mass g (reported in kg).

A `Unit` is one reactor compartment:
  (a) n_cells = 1  -> completely mixed reactor (CSTR)
  (b) n_cells > 1  -> advective-dispersive plug-flow compartment (finite volume,
                      tanks in series + optional physical dispersion `disp_ratio` = alpha_L / L)
Processes
  NO3-N  : 'off' | 'first_order'  r = k20 * th^(T-20) * C
                 | 'mm'           r = vmax * th^(T-20) * C / (K + C)      (Michaelis-Menten)
           optional biofilm term (first-order biofilm with diffusion limitation,
           effectiveness tanh(phi)/phi, and external mass-transfer resistance)
  PO4-P  : 'off' | 'first_order'  permanent removal r = k20 * th^(T-20) * C
                 | 'langmuir'     kinetic Langmuir sorption with FINITE capacity qmax that
                                  depletes:  dq/dt = th^(T-20) (kads*C*(qmax-q) - kdes*q)
                 | 'relax'        relaxation to the Langmuir isotherm (AQUASIM-style):
                                  dq/dt = krel (qmax*C/(K+C) - q)
High-flow bypass: inflow above `q_design` bypasses the unit and re-joins the outlet.

MASS BALANCE FIRST: every call to `simulate_unit` / `simulate_train` computes the water,
N and P balances and raises `MassBalanceError` if a closure residual exceeds `tol`.
"""
from __future__ import annotations
from dataclasses import dataclass, field, replace
import numpy as np
import pandas as pd
from . import kernel_AN_v1 as K


class MassBalanceError(RuntimeError):
    """Raised when a water, N or P balance does not close within tolerance."""


@dataclass
class Forcing:
    """Piecewise-constant forcing, one value per interval of length `dt` (days)."""
    Q: np.ndarray        # inflow m3/d
    T: np.ndarray        # temperature degC
    cinN: np.ndarray     # NO3-N g/m3
    cinP: np.ndarray     # PO4-P g/m3
    dt: float = 1.0
    start: str = "2025-01-01"

    def __post_init__(self):
        self.Q, self.T, self.cinN, self.cinP = [np.ascontiguousarray(np.asarray(a, float))
                                                for a in (self.Q, self.T, self.cinN, self.cinP)]

    def __len__(self):
        return len(self.Q)


@dataclass
class Unit:
    name: str = "unit"
    n_cells: int = 1
    volume: float = 10.0            # pore-water volume (m3)
    q_design: float = 1e9           # design flow (m3/d); inflow above this bypasses
    disp_ratio: float = 0.0         # alpha_L / L  (dimensionless)
    scheme: str = "upwind"          # 'upwind' | 'central' | 'tvd' (van Leer flux limiter)
    outlet: str = "zero_gradient"   # 'zero_gradient' (Danckwerts) | 'transmission' (AQUASIM, d2C/dx2=0; central only)
    n_kinetics: str = "off"         # 'off' | 'first_order' | 'mm'
    k_n20: float = 0.0              # 1/d
    vmax_n: float = 0.0             # g N/m3/d (per m3 pore water)
    K_n: float = 1.0                # g/m3
    theta_n: float = 1.0
    p_kinetics: str = "off"         # 'off' | 'first_order' | 'langmuir' | 'relax'
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
        if self.outlet == 'transmission' and self.scheme != 'central':
            raise ValueError("outlet='transmission' requires scheme='central'")

    def with_(self, **kw) -> "Unit":
        return replace(self, **kw)

    # -- packing for the kernel ---------------------------------------------------
    def _pack(self):
        p = np.zeros(K.NP)
        p[K.P_V] = self.volume; p[K.P_DISP] = self.disp_ratio; p[K.P_QDES] = self.q_design
        p[K.P_KN20] = self.k_n20; p[K.P_VMAXN] = self.vmax_n; p[K.P_KN] = self.K_n
        p[K.P_THN] = self.theta_n
        p[K.P_KP20] = self.k_p20; p[K.P_KADS] = self.k_ads; p[K.P_KDES] = self.k_des
        p[K.P_QMAX] = self.q_max; p[K.P_KP] = self.K_p; p[K.P_KREL] = self.k_rel
        p[K.P_THP] = self.theta_p; p[K.P_RHOB_TH] = self.rho_bulk_over_theta
        p[K.P_BIO_A] = self.bio_area; p[K.P_BIO_KV] = self.bio_kv20; p[K.P_BIO_DF] = self.bio_Df
        p[K.P_BIO_LF] = self.bio_Lf; p[K.P_BIO_KL] = self.bio_kL; p[K.P_THBIO] = self.bio_theta
        p[K.P_KD] = self.K_d; p[K.P_OUTMODE] = 1.0 if self.outlet == 'transmission' else 0.0; p[K.P_Q0] = self.q0; p[K.P_CN0] = self.cN0; p[K.P_CP0] = self.cP0; p[K.P_TREF] = self.t_ref
        return p

    def _modes(self):
        nm = {"off": 0, "first_order": 1, "mm": 2}[self.n_kinetics]
        pm = {"off": 0, "first_order": 1, "langmuir": 2, "relax": 3, "relax_linear": 4}[self.p_kinetics]
        return nm, pm, {"upwind": 0, "central": 1, "tvd": 2}[self.scheme]

    def stiffness_rate(self, f: Forcing) -> float:
        """Upper estimate of the fastest rate (1/d) used to choose RK4 sub-steps."""
        n = self.n_cells
        Vc = self.volume / n
        qmax_flow = min(float(np.max(f.Q)), self.q_design)
        tmax = float(np.max(f.T))
        lam = qmax_flow * (2 + 4 * self.disp_ratio * n) / Vc      # upwind eigenvalue ~2Q/Vc
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
            cmax = float(max(np.max(f.cinP), self.cP0, 1e-9))
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
    raw: np.ndarray                 # kernel output per interval
    nsub: int
    qdt: np.ndarray = field(default=None)

    # ---- per-interval series ---------------------------------------------------------
    def series(self, name: str) -> np.ndarray:
        """State at the END of each interval: 'cN_last','cP_last','cN_first','cP_first','cN_mean','cP_mean'."""
        col = {"cN_last": K.R_CNL, "cP_last": K.R_CPL, "cN_first": K.R_CNF, "cP_first": K.R_CPF,
               "cN_mean": K.R_CNM, "cP_mean": K.R_CPM, "cN_last2": K.R_CNL2, "cP_last2": K.R_CPL2}[name]
        return self.raw[:, col]

    @property
    def cout_N(self):               # flow-weighted outlet concentration incl. bypass (g/m3)
        v = self.raw[:, K.R_VIN]
        with np.errstate(invalid="ignore", divide="ignore"):
            c = (self.raw[:, K.R_NOUT] + self.raw[:, K.R_NBY]) / v
        return np.where(v > 0, c, self.forcing.cinN)

    @property
    def cout_P(self):
        v = self.raw[:, K.R_VIN]
        with np.errstate(invalid="ignore", divide="ignore"):
            c = (self.raw[:, K.R_POUT] + self.raw[:, K.R_PBY]) / v
        return np.where(v > 0, c, self.forcing.cinP)

    @property
    def cout_treated_N(self):
        v = self.raw[:, K.R_VTR]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(v > 0, self.raw[:, K.R_NOUT] / v, np.nan)

    @property
    def cout_treated_P(self):
        v = self.raw[:, K.R_VTR]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(v > 0, self.raw[:, K.R_POUT] / v, np.nan)

    # ---- mass balance --------------------------------------------------------------------
    def balance(self) -> dict:
        r = self.raw
        S = lambda c: float(np.sum(r[:, c]))
        w_in, w_tr, w_by = S(K.R_VIN), S(K.R_VTR), S(K.R_VBY)
        out = dict(
            water_in_m3=w_in, water_treated_out_m3=w_tr, water_bypass_m3=w_by,
            water_storage_change_m3=0.0,
            water_residual_m3=w_in - w_tr - w_by - 0.0,
            N_in_g=S(K.R_NIN), N_out_treated_g=S(K.R_NOUT), N_bypass_g=S(K.R_NBY),
            N_removed_g=S(K.R_NREM), N_storage_change_g=S(K.R_NSTO),
            P_in_g=S(K.R_PIN), P_out_treated_g=S(K.R_POUT), P_bypass_g=S(K.R_PBY),
            P_removed_g=S(K.R_PREM), P_sorbed_g=S(K.R_PSORB), P_storage_change_g=S(K.R_PSTO),
        )
        out["N_residual_g"] = out["N_in_g"] - out["N_out_treated_g"] - out["N_bypass_g"] \
            - out["N_removed_g"] - out["N_storage_change_g"]
        out["P_residual_g"] = out["P_in_g"] - out["P_out_treated_g"] - out["P_bypass_g"] \
            - out["P_removed_g"] - out["P_storage_change_g"]
        # residuals are normalised by the largest term of the balance (robust when inflow load is 0)
        sc = lambda *a: max(max(abs(x) for x in a), 1e-300)
        out["water_rel_residual"] = abs(out["water_residual_m3"]) / sc(w_in, w_tr, w_by)
        out["N_rel_residual"] = abs(out["N_residual_g"]) / sc(
            out["N_in_g"], out["N_out_treated_g"], out["N_bypass_g"], out["N_removed_g"], out["N_storage_change_g"])
        out["P_rel_residual"] = abs(out["P_residual_g"]) / sc(
            out["P_in_g"], out["P_out_treated_g"], out["P_bypass_g"], out["P_removed_g"], out["P_storage_change_g"])
        out["P_sorbed_end_g"] = float(self.raw[-1, K.R_PQMASS])
        days = len(self.forcing) * self.forcing.dt
        out["period_days"] = days
        scale = 365.0 / days / 1000.0                     # g over run -> kg per year
        out["kgN_retained_per_yr"] = out["N_removed_g"] * scale
        out["kgP_retained_per_yr"] = out["P_removed_g"] * scale
        out["kgN_in_per_yr"] = out["N_in_g"] * scale
        out["kgP_in_per_yr"] = out["P_in_g"] * scale
        out["kgN_bypass_per_yr"] = out["N_bypass_g"] * scale
        out["kgP_bypass_per_yr"] = out["P_bypass_g"] * scale
        return out

    def check_balance(self, tol: float = 1e-6) -> dict:
        b = self.balance()
        bad = [k for k in ("water_rel_residual", "N_rel_residual", "P_rel_residual") if b[k] > tol]
        if bad:
            raise MassBalanceError(
                f"[{self.unit.name}] mass balance failed (tol={tol:g}): "
                + ", ".join(f"{k}={b[k]:.3e}" for k in bad))
        return b

    def balance_table(self, freq: str = "MS") -> pd.DataFrame:
        """Balance per calendar period (default month) in kg / m3 plus closure residuals."""
        idx = pd.date_range(self.forcing.start, periods=len(self.forcing),
                            freq=f"{self.forcing.dt:g}D")
        d = pd.DataFrame(self.raw, index=idx)
        g = d.groupby(pd.Grouper(freq=freq)).sum()
        t = pd.DataFrame(index=g.index)
        t["water_in_m3"] = g[K.R_VIN]; t["water_treated_out_m3"] = g[K.R_VTR]
        t["water_bypass_m3"] = g[K.R_VBY]
        t["water_residual_m3"] = t.water_in_m3 - t.water_treated_out_m3 - t.water_bypass_m3
        for s, (i, o, b, r, st) in {"N": (K.R_NIN, K.R_NOUT, K.R_NBY, K.R_NREM, K.R_NSTO),
                                   "P": (K.R_PIN, K.R_POUT, K.R_PBY, K.R_PREM, K.R_PSTO)}.items():
            t[f"{s}_in_kg"] = g[i] / 1e3; t[f"{s}_out_treated_kg"] = g[o] / 1e3
            t[f"{s}_bypass_kg"] = g[b] / 1e3; t[f"{s}_removed_kg"] = g[r] / 1e3
            t[f"{s}_storage_change_kg"] = g[st] / 1e3
            t[f"{s}_residual_kg"] = (g[i] - g[o] - g[b] - g[r] - g[st]) / 1e3
        t.insert(0, "unit", self.unit.name)
        return t


def choose_nsub(unit: Unit, f: Forcing, nsub_min: int = 4, safety: float = 2.0) -> int:
    lam = unit.stiffness_rate(f)
    return int(max(nsub_min, np.ceil(f.dt * lam / safety)))


def simulate_unit(unit: Unit, f: Forcing, nsub: int | None = None, tol: float = 1e-6,
                  check: bool = True) -> Result:
    nm, pm, sc = unit._modes()
    ns = nsub or choose_nsub(unit, f)
    raw = K.simulate_kernel(f.Q, f.T, f.cinN, f.cinP, float(f.dt), int(ns), unit._pack(),
                            int(unit.n_cells), sc, nm, pm)
    res = Result(unit, f, raw, ns)
    if check:
        res.check_balance(tol)
    return res


def simulate_train(units: list[Unit], f: Forcing, nsub: int | None = None, tol: float = 1e-6):
    """Units in series: the (mixed treated + bypass) outlet of one unit feeds the next.
    Returns list of Results; the train outlet is that of the last unit."""
    results = []
    cur = f
    for u in units:
        r = simulate_unit(u, cur, nsub=nsub, tol=tol)
        results.append(r)
        cur = Forcing(cur.Q, cur.T, r.cout_N, r.cout_P, cur.dt, cur.start)
    return results


def train_balance(results: list[Result], tol: float = 1e-6) -> dict:
    """Overall balance of a train (inflow to unit 1 vs outlet of last unit + all sinks)."""
    first, last = results[0].balance(), results[-1]
    f0 = results[0].forcing
    vin = float(np.sum(f0.Q) * f0.dt)
    outN = float(np.sum(last.raw[:, K.R_NOUT] + last.raw[:, K.R_NBY]))
    outP = float(np.sum(last.raw[:, K.R_POUT] + last.raw[:, K.R_PBY]))
    inN, inP = first["N_in_g"], first["P_in_g"]
    remN = sum(r.balance()["N_removed_g"] for r in results)
    remP = sum(r.balance()["P_removed_g"] for r in results)
    stoN = sum(r.balance()["N_storage_change_g"] for r in results)
    stoP = sum(r.balance()["P_storage_change_g"] for r in results)
    out = dict(
        water_in_m3=vin,
        water_out_m3=float(np.sum(last.raw[:, K.R_VIN])),
        N_in_g=inN, N_out_total_g=outN, N_removed_g=remN, N_storage_change_g=stoN,
        P_in_g=inP, P_out_total_g=outP, P_removed_g=remP, P_storage_change_g=stoP)
    out["N_residual_g"] = inN - outN - remN - stoN
    out["P_residual_g"] = inP - outP - remP - stoP
    out["N_rel_residual"] = abs(out["N_residual_g"]) / inN
    out["P_rel_residual"] = abs(out["P_residual_g"]) / inP
    out["water_rel_residual"] = abs(vin - out["water_out_m3"]) / vin
    bad = [k for k in ("water_rel_residual", "N_rel_residual", "P_rel_residual") if out[k] > tol]
    if bad:
        raise MassBalanceError("train mass balance failed: " + ", ".join(f"{k}={out[k]:.3e}" for k in bad))
    days = len(f0) * f0.dt
    sc = 365.0 / days / 1000.0
    out["kgN_retained_per_yr"] = remN * sc
    out["kgP_retained_per_yr"] = remP * sc
    return out
