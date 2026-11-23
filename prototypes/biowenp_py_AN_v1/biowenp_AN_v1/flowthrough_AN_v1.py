# Part of biowenp_py_AN_v1 (LGPL-3.0-or-later). AUGMENTATION of BioWEN-P (Zhou et al. 2023) for saturated
# flow-through eggshell + woodchip reactors (Maria Dittrich's CWA FEI proposal, modelling section v4).
"""Saturated flow-through reactor with design-flow bypass, N and P modules and a hard mass-balance gate.

Equations (numbers refer to the modelling section v4):
  bypass (4):         Q_t = min(Q, Q_design), Q_b = Q - Q_t
  outlet (5):         C_out = (Q_t C_t + Q_b C_in) / Q
  transport (6)-(8):  n completely mixed cells, V_c dC_i/dt = Q_t (C_{i-1}-C_i) + E (C_{i+1}-2C_i+C_{i-1}) - V_c r_i,
                      E = (alpha_L/L) n Q_t ; Danckwerts (total-flux) inlet, zero-gradient outlet; n = 1 -> CSTR
  temperature (9):    f(T) = theta^(T - 20)
  nitrate (10)/(11):  r_N = k_N20 f_N C_N   or   r_N = vmax_N20 f_N C_N/(K_N + C_N)      (woodchip denitrification)
  SRP sorption:       'zhou'  : R_ad = k_ad f_P (C - C_eq(q)),  C_eq = q/(K_L (q_max - q))   (BioWEN-P T4.13)
                      'maria' : dq/dt = f_P (k_ads C (q_max - q) - k_des q)                  (v4 eq. 12-13)
                      capacity depletes as q -> q_max; optional zero-order sorbed -> mineral transfer k_mi
                      (BioWEN-P precipitation pathway; frees sorption sites).
  particulate P:      colloid filtration R_filt = lambda v C_PP (BioWEN-P T4.12), v = Q_t/(A phi)

Numerics: explicit conservative finite-volume sub-steps. Every reaction is applied as a mass TRANSFER between
two booked pools, limited so that it can neither overshoot the closed-cell equilibrium nor drive a pool negative.
Hence every transfer appears once with + and once with -, and the ledger closes to round-off.

Mass-balance gate: every simulate_* call checks water, N and P closure (default rtol 1e-6) and raises
MassBalanceError otherwise. Residuals and the full input table are attached to the Result.
Units: d, m3, g/m3 (= mg/L), kg media, q in g P per kg sorbing media.
"""
from dataclasses import dataclass, asdict, field
import numpy as np
import pandas as pd
from .massbalance_AN_v1 import MassBalanceError, gate


@dataclass
class FlowThroughUnit:
    name: str = "unit"
    n_cells: int = 4
    bed_volume: float = 20.0          # m3 bulk bed
    porosity: float = 0.6
    length: float = 10.0              # m (flow path)
    q_design: float = 50.0            # m3/d
    alpha_L: float = 0.1              # m dispersivity
    # nitrate (woodchip)
    n_kinetics: str = "none"          # 'none' | 'first' | 'mm'
    k_N20: float = 0.0                # 1/d
    vmax_N20: float = 0.0             # g N/m3 pore/d
    K_N: float = 1.0                  # g N/m3
    theta_N: float = 1.08
    # SRP sorption (eggshell)
    p_kinetics: str = "none"          # 'none' | 'zhou' | 'maria'
    sorbent_mass: float = 0.0         # kg sorbing media (e.g. eggshell) in the unit
    q_max: float = 0.0                # g P/kg sorbent
    K_L: float = 1.0                  # m3/g (= L/mg)  ('zhou')
    k_ad: float = 0.0                 # 1/d            ('zhou')
    k_ads: float = 0.0                # m3/g/d         ('maria')
    k_des: float = 0.0                # 1/d            ('maria')
    k_mi: float = 0.0                 # g P/kg/d zero-order sorbed -> mineral
    theta_P: float = 1.0
    q0: float = 0.0                   # initial sorbed P, g/kg
    # particulate P
    lambda_pp: float = 0.0            # 1/m filtration coefficient
    T_ref: float = 20.0

    @property
    def pore_volume(self):
        return self.bed_volume * self.porosity

    @property
    def area(self):
        return self.bed_volume / self.length

    def table(self):
        return pd.DataFrame(list(asdict(self).items()), columns=["parameter", "value"])


@dataclass
class ReactorForcing:
    Q: np.ndarray        # m3/d
    T: np.ndarray        # degC
    C_N: np.ndarray      # nitrate-N g/m3
    C_SRP: np.ndarray    # g/m3
    C_PP: np.ndarray     # particulate P g/m3

    def table(self):
        return pd.DataFrame(dict(day=np.arange(len(self.Q)), Q_m3d=self.Q, T_C=self.T, C_N=self.C_N,
                                 C_SRP=self.C_SRP, C_PP=self.C_PP))


def _langmuir_ceq_closed(total, Vc, m, qmax, K):
    """Closed-cell equilibrium C* for total P = Vc C + m qmax K C/(1+K C) (positive root)."""
    a = Vc * K
    b = Vc + m * qmax * K - total * K
    c = -total
    disc = np.sqrt(b * b - 4 * a * c)
    return np.where(a > 0, (-b + disc) / (2 * np.where(a > 0, a, 1)), total / Vc)


@dataclass
class Result:
    unit: FlowThroughUnit
    forcing: ReactorForcing
    daily: pd.DataFrame
    ledger: dict
    profiles: dict

    def balance(self):
        return self.ledger

    def residuals(self):
        L = self.ledger
        return {
            "water": (L["V_in"] - L["V_t"] - L["V_b"], L["V_in"]),
            "N": (L["N_in"] - L["N_out_t"] - L["N_b"] - L["N_rem"] - L["dS_N"], L["N_in"]),
            "P": (L["P_in"] - L["P_out_t"] - L["P_b"] - L["P_sorbed"] - L["P_mineral"] - L["P_filtered"] - L["dS_P"],
                  L["P_in"]),
        }

    def check_balance(self, rtol=1e-6):
        return gate(self.residuals(), rtol)

    def balance_table(self):
        rows = [dict(term=k, value=v) for k, v in self.ledger.items()]
        for k, (r, s) in self.residuals().items():
            rows.append(dict(term=f"residual_{k}", value=r))
            rows.append(dict(term=f"residual_{k}_rel", value=r / s if s else 0.0))
        return pd.DataFrame(rows)


def simulate_flowthrough(u: FlowThroughUnit, f: ReactorForcing, nsub=None, rtol=1e-6, keep_profiles=True):
    n = u.n_cells
    Vc = u.pore_volume / n
    mc = u.sorbent_mass / n
    D = len(f.Q)
    qmax_flow = min(float(np.max(f.Q)), u.q_design)
    E_max = (u.alpha_L / u.length) * n * qmax_flow
    if nsub is None:
        nsub = int(np.ceil((qmax_flow + 2 * E_max) / Vc / 0.5)) + 1
        nsub = max(nsub, 24)
    dt = 1.0 / nsub
    CN = np.zeros(n); CP = np.zeros(n); CPP = np.zeros(n)
    q = np.full(n, u.q0); qm = np.zeros(n); filt = np.zeros(n)
    S0_N = CN.sum() * Vc
    S0_P = (CP.sum() + CPP.sum()) * Vc
    q_init = q.copy()
    L = dict(V_in=0.0, V_t=0.0, V_b=0.0, N_in=0.0, N_out_t=0.0, N_b=0.0, N_rem=0.0, P_in=0.0, P_out_t=0.0, P_b=0.0,
             P_ads_gross=0.0, P_des_gross=0.0, P_min_transfer=0.0, P_filt_transfer=0.0)
    out = {k: np.zeros(D) for k in ["Q", "Q_t", "Q_b", "C_N_in", "C_N_out", "C_P_in", "C_P_out", "C_PP_out",
                                    "N_rem_g", "P_sorb_g", "q_mean", "q_frac_capacity", "T"]}
    prof_q = np.zeros((D, n)); prof_CP = np.zeros((D, n)); prof_CN = np.zeros((D, n))
    v_cells = np.arange(n)
    for d in range(D):
        Q = f.Q[d]; T = f.T[d]
        Qt = min(Q, u.q_design); Qb = Q - Qt
        E = (u.alpha_L / u.length) * n * Qt
        fN = u.theta_N ** (T - u.T_ref); fP = u.theta_P ** (T - u.T_ref)
        vel = Qt / (u.area * u.porosity)  # m/d
        nrem_day = 0.0; psorb_day = 0.0
        day_out = {"N": 0.0, "P": 0.0, "PP": 0.0}           # treated outflow mass this day (g)
        for _ in range(nsub):
            # ---- transport (conservative): mass changes = dt*(fluxes)
            for C, cin, key_in, key_out, sp in ((CN, f.C_N[d], "N_in", "N_out_t", "N"),
                                                (CP, f.C_SRP[d], "P_in", "P_out_t", "P"),
                                                (CPP, f.C_PP[d], "P_in", "P_out_t", "PP")):
                up = np.concatenate(([cin], C[:-1]))
                adv = Qt * (up - C)                      # upwind; column sum telescopes to Qt (cin - C_last)
                disp = np.zeros(n)
                if n > 1 and E > 0:
                    fl = E * (C[1:] - C[:-1])           # dispersive exchange between i and i+1 (sums to 0)
                    disp[:-1] += fl
                    disp[1:] -= fl
                L[key_in] += Qt * cin * dt
                L[key_out] += Qt * C[-1] * dt           # advective outflow uses the pre-update outlet cell
                day_out[sp] += Qt * C[-1] * dt
                C += dt * (adv + disp) / Vc
            # ---- nitrate removal (woodchip)
            if u.n_kinetics != "none":
                if u.n_kinetics == "first":
                    r = u.k_N20 * fN * CN
                else:
                    r = u.vmax_N20 * fN * CN / (u.K_N + CN)
                tr = np.minimum(r * dt, CN)              # g/m3 removed this sub-step
                CN -= tr
                g = tr.sum() * Vc
                L["N_rem"] += g; nrem_day += g
            # ---- SRP sorption (eggshell)
            if u.p_kinetics != "none" and mc > 0:
                total = CP * Vc + mc * q
                if u.p_kinetics == "zhou":
                    K = u.K_L
                    ceq_q = np.where(q < u.q_max, q / (u.K_L * np.maximum(u.q_max - q, 1e-300)), np.inf)
                    rate = u.k_ad * fP * (CP - ceq_q)                       # g/m3/d (+ = adsorption)
                else:
                    K = u.k_ads / u.k_des if u.k_des > 0 else np.inf
                    rate = (u.k_ads * CP * (u.q_max - q) - u.k_des * q) * fP * mc / Vc
                want = rate * dt * Vc                                     # g
                if np.isfinite(K):
                    cstar = _langmuir_ceq_closed(total, Vc, mc, u.q_max, K)
                    lim = (CP - cstar) * Vc                               # transfer that reaches equilibrium
                else:
                    lim = np.minimum(CP * Vc, mc * (u.q_max - q))
                want = np.where(np.isfinite(want), want, np.sign(lim) * np.abs(lim))
                tr = np.where(lim >= 0, np.clip(want, 0, lim), np.clip(want, lim, 0))
                tr = np.minimum(tr, CP * Vc)                              # cannot take more than dissolved
                tr = np.maximum(tr, -mc * q)                              # cannot release more than sorbed
                CP -= tr / Vc
                q += tr / mc
                L["P_ads_gross"] += tr[tr > 0].sum(); L["P_des_gross"] += -tr[tr < 0].sum()
                psorb_day += tr.sum()
                if u.k_mi > 0:
                    pm = np.minimum(u.k_mi * dt, q)
                    q -= pm; qm += pm
                    L["P_min_transfer"] += pm.sum() * mc
            # ---- particulate P filtration
            if u.lambda_pp > 0 and vel > 0:
                frac = 1 - np.exp(-u.lambda_pp * vel * dt)
                tr = CPP * frac
                CPP -= tr
                filt += tr * Vc
                L["P_filt_transfer"] += tr.sum() * Vc
        # daily record
        L["V_in"] += Q; L["V_t"] += Qt; L["V_b"] += Qb
        L["N_b"] += Qb * f.C_N[d]; L["P_b"] += Qb * (f.C_SRP[d] + f.C_PP[d])
        out["Q"][d], out["Q_t"][d], out["Q_b"][d], out["T"][d] = Q, Qt, Qb, T
        out["C_N_in"][d], out["C_P_in"][d] = f.C_N[d], f.C_SRP[d] + f.C_PP[d]
        # eq. (5) with the day-integrated treated outflow (flow-weighted daily mean), so that loads passed to a
        # downstream unit equal the loads that left this one exactly
        out["C_N_out"][d] = (day_out["N"] + Qb * f.C_N[d]) / Q if Q > 0 else 0.0
        out["C_P_out"][d] = (day_out["P"] + Qb * f.C_SRP[d]) / Q if Q > 0 else 0.0
        out["C_PP_out"][d] = (day_out["PP"] + Qb * f.C_PP[d]) / Q if Q > 0 else 0.0
        out["N_rem_g"][d] = nrem_day; out["P_sorb_g"][d] = psorb_day
        out["q_mean"][d] = q.mean()
        out["q_frac_capacity"][d] = q.mean() / u.q_max if u.q_max > 0 else 0.0
        prof_q[d] = q; prof_CP[d] = CP; prof_CN[d] = CN
    return _finish(u, f, L, out, CN, CP, CPP, q, q_init, qm, filt, Vc, mc, S0_N, S0_P,
                   dict(q=prof_q, C_SRP=prof_CP, C_N=prof_CN) if keep_profiles else {}, rtol)


def _finish(u, f, L, out, CN, CP, CPP, q, q_init, qm, filt, Vc, mc, S0_N, S0_P, prof, rtol):
    L["N_in_t"], L["P_in_t"] = L["N_in"], L["P_in"]          # loads entering the treated path
    L["N_in"] = L["N_in_t"] + L["N_b"]                         # total inflow load (treated + bypassed)
    L["P_in"] = L["P_in_t"] + L["P_b"]
    L["dS_N"] = CN.sum() * Vc - S0_N
    L["dS_P"] = (CP.sum() + CPP.sum()) * Vc - S0_P
    L["P_sorbed"] = float(((q - q_init) * mc).sum())          # net change of the sorbed pool (g)
    L["P_mineral"] = float((qm * mc).sum())                   # mineral pool (g)
    L["P_filtered"] = float(filt.sum())                       # filtered particulate pool (g)
    L["P_rem"] = L["P_sorbed"] + L["P_mineral"] + L["P_filtered"]
    yrs = len(f.Q) / 365.0
    L["kgN_retained_per_yr"] = L["N_rem"] / 1000 / yrs
    L["kgP_retained_per_yr"] = L["P_rem"] / 1000 / yrs
    L["final_capacity_used_frac"] = float(q.mean() / u.q_max) if u.q_max > 0 else 0.0
    daily = pd.DataFrame(out)
    r = Result(u, f, daily, L, prof)
    r.check_balance(rtol)                                     # HARD GATE on every run
    return r


def forcing_from_result(r: Result):
    """Outlet of one unit (bypass already remixed, eq. 5) as the forcing of the next unit in a train."""
    d = r.daily
    return ReactorForcing(Q=d.Q.values, T=d["T"].values, C_N=d.C_N_out.values, C_SRP=d.C_P_out.values,
                          C_PP=d.C_PP_out.values)


def simulate_train(units, f: ReactorForcing, rtol=1e-6, **kw):
    """Units in series (e.g. woodchip -> eggshell). Each unit is gated; the train ledger is gated too."""
    results = []
    fi = f
    for u in units:
        r = simulate_flowthrough(u, fi, rtol=rtol, **kw)
        results.append(r)
        fi = forcing_from_result(r)
    # train-level closure: inflow = final outlet + all retained + storage change
    first, last = results[0], results[-1]
    N_in = first.ledger["N_in"]
    P_in = first.ledger["P_in"]
    N_out = last.ledger["N_out_t"] + last.ledger["N_b"]
    P_out = last.ledger["P_out_t"] + last.ledger["P_b"]
    N_ret = sum(r.ledger["N_rem"] + r.ledger["dS_N"] for r in results)
    P_ret = sum(r.ledger["P_rem"] + r.ledger["dS_P"] for r in results)
    train = dict(N_in=N_in, N_out=N_out, N_retained_and_stored=N_ret, P_in=P_in, P_out=P_out,
                 P_retained_and_stored=P_ret, residual_N=N_in - N_out - N_ret, residual_P=P_in - P_out - P_ret)
    gate({"train_N": (train["residual_N"], N_in), "train_P": (train["residual_P"], P_in)}, rtol)
    return results, train
