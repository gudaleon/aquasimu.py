# Copyright 2020-2023 Bowen Zhou (original MATLAB BioWEN-P); Python port 2026 Alex Neumann / EnviroManagement.
# Part of biowenp_py_AN_v1, a Python port of BioWEN-P (Zhou et al. 2023, doi:10.1016/j.scitotenv.2023.162749;
# code: FRDR doi:10.20383/103.0640). GNU LGPL v3.0 or later (see COPYING.LESSER, COPYING). NO WARRANTY.
"""Hydrologic module: line-by-line port of the main loop of ``ModelMainloop.m`` (ElmDrive == 1 branch),
``InflowEstimator.m`` and ``Discharge.m`` (method 1).

Explicit Euler, 10-minute step (dt = 1).  Supported switches: InfiltrationMethod 1, PETMethod 0/1,
FlowEstimateMethod 1, Inflow_TP_EstimateMethod 1/2/3.  Not ported: Green-Ampt (unused in the paper),
Priestley-Taylor (PETMethod 2), the neural-network rating curve (FlowEstimateMethod 2, stored as a MATLAB
network object in Qd_Sr_Net.mat), the SWMM-input branch (ElmDrive == 0).

Quirks of the original are reproduced deliberately (the purpose is equivalence) and flagged with "QUIRK".
"""
import numpy as np
from dataclasses import dataclass


@dataclass
class HydroForcing:
    precip: np.ndarray     # mm per 10 min
    temp: np.ndarray       # air temperature, degC
    month: np.ndarray      # calendar month 1..12 of each step

    def __post_init__(self):
        self.precip = np.nan_to_num(np.asarray(self.precip, float), nan=0.0)      # as ModelMainloop.m
        self.temp = np.asarray(self.temp, float).copy()
        self.temp[np.isnan(self.temp)] = 3.0
        self.month = np.asarray(self.month, int)


def discharge(Sr, Sr_crit, method=1):
    """Discharge.m, method 1 (nonlinear rating curve). Returns underdrain discharge in L per 10 min."""
    if method != 1:
        raise NotImplementedError("FlowEstimateMethod 2 (MATLAB neural net) is not portable; use method 1")
    if Sr >= Sr_crit:
        head = Sr - Sr_crit
        if head <= 0.2:
            return 77.721 * np.exp(28.418 * head)
        return 20245.0 * np.log(head) + 42646.0
    return 0.0


def inflow_estimator(m_B, DA, la, Q_in, dt, method, c_const=4.3):
    """InflowEstimator.m. Q_in is a volume (m3 per step). Returns (m_B, C_TP_in mg/L)."""
    if method == 1:
        C_B = 0.2 / 24 / 6 / 1000 * 170
        K_w, N_w = 7.71, 1.8
        if Q_in == 0:
            return m_B + C_B * dt * DA * la, 0.0
        k = K_w * Q_in ** N_w
        w = -k * m_B
        return max(0.0, m_B + w * dt), -w / Q_in / 1000
    if method == 2:
        C_1, C_2, C_3, C_4 = 3e6 / 1e4, 0.06 / 24 / 6, 500.0, 2.35
        if Q_in == 0:
            return max(0.0, m_B + (C_1 * C_2 - C_2 * m_B) * dt), 0.0
        w = -C_3 * Q_in ** C_4 * m_B
        return max(0.0, m_B + w * dt), -w / Q_in / 1000 + 0.35
    return m_B, (0.0 if Q_in == 0 else c_const)


def monthly_pet(temp, month, matlab_quirk=True):
    """Thornthwaite monthly PET (m/month) and monthly mean temperature as in ModelMainloop.m lines 32-45.

    QUIRK (found by the Octave equivalence run): in the original, Temp_MonthlyAve is a 12x1 column and b_PET a
    1x12 row, so ``b_PET.*(...)`` broadcasts to a 12x12 matrix and ``PET_Month(m)`` (linear index) picks column 1,
    i.e. the daylength factor of JANUARY (0.82) is used for every month.  ``matlab_quirk=True`` reproduces this;
    False applies the intended month-specific factor."""
    tavg = np.full(12, np.nan)
    for m in range(1, 13):
        sel = temp[month == m]
        if sel.size:
            tavg[m - 1] = np.nanmean(sel)
    tpos = np.fmax(0.0, tavg)                      # MATLAB max(0,NaN) = 0
    I_PET = np.nansum((tpos / 5.0) ** 1.51)
    a = 67.5e-8 * I_PET ** 3 - 77.1e-6 * I_PET ** 2 + 0.0179 * I_PET + 0.492
    b = np.array([0.82, 0.82, 1.025, 1.12, 1.26, 1.27, 1.29, 1.20, 1.04, 0.95, 0.81, 0.88])
    bb = np.full(12, b[0]) if matlab_quirk else b
    pet = 0.01 * 1.62 * bb * (10.0 * tpos / I_PET) ** a
    return pet, tavg


def run_hydrology(forc: HydroForcing, h, dt=1.0, matlab_quirks=True):
    """Port of the ElmDrive hydrologic loop. Returns a dict of arrays of length stime (as Model_Output)."""
    P, T, mon = forc.precip, forc.temp, forc.month
    stime = len(P)
    PET_Month, Tavg = monthly_pet(T, mon, matlab_quirks)
    Sr_crit = h.Sr_crit
    z = lambda: np.zeros(stime + 1)
    Ss, Qm, Q_in, Q_in_d, Q_in_c, s, I, ET, L, Sp, Qo, Qd, Sr, Q_Ex, m_B, C_TP_in = (z() for _ in range(16))
    Void, Retain = z(), z()
    s[0] = h.sfc
    for i in range(stime):
        n = i
        Rv = 0.1 if (T[n] > 22 and P[n] < 15) else 0.05 + h.f * h.la
        if T[n] <= h.Tt:
            Qm[i] = 0.0
            Ss[i + 1] = Ss[i] + P[n] / 1000 * dt
            Q_in[i] = (h.DA * Rv + h.A_BMP) * Qm[i]
            Q_in_c[i] = (h.DA * Rv) * Qm[i] / h.A_BC
        else:
            Qm[i] = min(Ss[i], h.DD * (T[n] - h.Tt) * dt)
            Ss[i + 1] = Ss[i] - Qm[i]
            Q_in[i] = h.DA_PL * Rv * Qm[i] + (h.DA * Rv + h.A_BMP) * (Qm[i] + P[n] / 1000 * dt)
            Q_in_c[i] = (h.DA * Rv) * (Qm[i] + P[n] / 1000 * dt) / h.A_BC
        Q_in_d[i] = Q_in[i] / h.A_BC
        Retain[i] = max(0.0, min(h.Void_max - Void[i], Q_in_c[i]))
        Q_in_d[i] -= Retain[i]
        m_B[i + 1], C_TP_in[i] = inflow_estimator(m_B[i], h.DA, h.la, Q_in_d[i] * h.A_BC, dt,
                                                  h.Inflow_TP_EstimateMethod, h.C_TP_in_const)
        # infiltration (method 1)
        if s[i] != 1:
            I[i] = min(Q_in_d[i] + Sp[i], h.n_fm * h.dfm * (1 - s[i]))
        else:
            I[i] = h.Ks * (Q_in_d[i] + Sp[i]) / h.dfm
        m = mon[n] - 1
        if h.PETMethod == 0:
            PET_now = PET_Month[m] / 30 / 24 / 6   # QUIRK: PET_now undefined for method 0 in the original
            ET[i] = max(0.0, s[i] - h.s_h) * PET_now if T[n] > 0 else 0.0
        elif h.PETMethod == 1:
            PET_now = max(0.0, PET_Month[m] * min(1 + 0.5 * (T[n] - Tavg[m]), 2) / 30 / 24 / 6)
            ET[i] = max(0.0, s[i] - h.s_h) * PET_now if T[n] > 0 else 0.0
        else:
            raise NotImplementedError("PETMethod 2 (Priestley-Taylor) not ported")
        Void[i + 1] = Retain[i] + (max(0.0, Void[i] - h.n_PET * PET_now) if T[n] > 0 else Void[i])
        # leakage
        if s[i] != 1:
            L[i] = max(0.0, h.Ks * (np.exp(h.beta * (s[i] - h.sfc)) - 1) / (np.exp(h.beta * (1 - h.sfc)) - 1))
        else:
            L[i] = I[i]
        Pref_L = L[i] * h.R_Pref if s[i] > 0.5 else 0.0
        L[i] -= Pref_L
        s[i + 1] = max(0.0, s[i] + (I[i] - ET[i] - L[i]) / h.n_fm / h.dfm)   # QUIRK: s is not capped at 1
        Sp[i + 1] = min(h.Sp_max, max(0.0, Sp[i] + Q_in_d[i] - I[i]))
        Qo[i] = 1000 * h.A_BC * max(0.0, Sp[i] + Q_in_d[i] - I[i] - h.Sp_max)
        Qd[i] = discharge(Sr[i], Sr_crit, h.FlowEstimateMethod) + Pref_L * h.A_BC * 1000
        Sr[i + 1] = max(0.0, Sr[i] + ((L[i] - Qd[i] / 1000 / h.A_BC) / h.n_bs - 1 / h.K_1 * Sr[i] * dt))
        Q_Ex[i] = 1 / h.K_1 * Sr[i] * h.A_BC * 1000 * h.n_bs * dt
        Qd[i] += Qo[i]
    out = dict(Ss=Ss, Qm=Qm, Q_in=Q_in, Q_in_d=Q_in_d, s=s, I=I, ET=ET, L=L, Sp=Sp, Qo=Qo, Qd=Qd, Sr=Sr,
               Q_Ex=Q_Ex, m_B=m_B, C_TP_in=C_TP_in, Void=Void, Retain=Retain)
    # "Remove the last row of simulated values" (all arrays of length stime+1 are truncated)
    out = {k: v[:stime] for k, v in out.items()}
    out["precip"], out["temp"] = P, T
    return out
