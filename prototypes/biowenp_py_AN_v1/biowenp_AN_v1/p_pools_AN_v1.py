# Copyright 2020-2023 Bowen Zhou (original MATLAB BioWEN-P); Python port 2026 Alex Neumann / EnviroManagement.
# Part of biowenp_py_AN_v1, a Python port of BioWEN-P (Zhou et al. 2023, doi:10.1016/j.scitotenv.2023.162749;
# code: FRDR doi:10.20383/103.0640). GNU LGPL v3.0 or later (see COPYING.LESSER, COPYING). NO WARRANTY.
"""Phosphorus pools and rate laws: port of ``BioC_P_RTM.m`` (BioWEN-P P reactive-transport model).

Pure rate-law functions (Zhou et al. 2023, Table 4) are exposed separately so they can be unit-tested and
re-used by the augmented flow-through model.  ``run_p_faithful`` reproduces the MATLAB time loop statement
by statement (including its operator-splitting inconsistencies, flagged "QUIRK"), so that MATLAB/Octave and
Python outputs can be compared array by array.
"""
import numpy as np
from .transport_AN_v1 import cn_step

# --------------------------------------------------------------------------- rate laws (Table 4)

def langmuir_q(C, Qmax, KL):
    """T4.13 (equilibrium part): adsorbed content (mg/kg) in equilibrium with C (mg/L)."""
    return Qmax * KL * C / (1.0 + KL * C)


def langmuir_ceq(q, Qmax, KL):
    """Inverse Langmuir: equilibrium porewater concentration for adsorbed content q (as coded: q/KL/(Qmax-q))."""
    return q / KL / (Qmax - q)


def adsorption_rate(C, q, k_ad, Qmax, KL):
    """T4.13 (kinetic part): R_Ad = k_ad (C - C_eq(q))  [mg/L per time]."""
    return k_ad * (C - langmuir_ceq(q, Qmax, KL))


def filtration_rate(C_PP, vx, lam):
    """T4.12: R_Filt = lambda v_x C_PP (colloid filtration, 1st order in depth)."""
    return lam * vx * C_PP


def bedding_filtration(C_PP_FM, lambda_s, d_s):
    """T4.16: C_PP leaving the bedding layer = C_PP_FM exp(-lambda_s d_s)."""
    return C_PP_FM * np.exp(-lambda_s * d_s)


def barber_uptake(C, I_max, Km, Cmin):
    """T4.14 (Barber 1984): I_max (C - Cmin)/(Km + C - Cmin) for C > Cmin, else 0 (per root-area basis)."""
    C = np.asarray(C, float)
    out = I_max * (C - Cmin) / (Km + C - Cmin)
    return np.where(C > Cmin, out, 0.0)


def hydrolysis_rate(M_labile, k_hy):
    """T4.15: R_Hy = k_Hy M_labile,POP (mg/kg per time)."""
    return k_hy * M_labile


def precipitation_zero_order(q_ad, k_mi_dt):
    """Zero-order transfer adsorbed -> mineral, limited by the adsorbed pool: min(k_mi dt, q_ad)."""
    return np.minimum(k_mi_dt, q_ad)


# --------------------------------------------------------------------------- faithful time loop

def run_p_faithful(hyd, h, p, dx=1.0, macro_percent=1.0, sorption_mode=2, C_init=0.0, dt=1.0, eps=0.5):
    """Port of BioC_P_RTM(soilmoisture, C0, C_init, I, dt, stime, dx, Macro_percent, SorptionEquilibrium, Model_Output).

    hyd : dict returned by hydrology_AN_v1.run_hydrology (needs s, ET, I, L, Sr, Qd, Qo, Q_Ex, C_TP_in).
    Returns dict with all pools (N x M arrays) and outflow series, same names as Model_Output fields.
    """
    soil = hyd["s"]; I = hyd["I"]; C0 = hyd["C_TP_in"]
    M = len(soil)
    N = int(round(h.dfm * 100 / dx))
    phi_macro = h.n_fm * macro_percent
    phi_micro = h.n_fm * (1 - macro_percent)
    with np.errstate(divide="ignore", invalid="ignore"):
        vx = I * 100 / phi_macro / soil
        ET_factor = soil / (soil - hyd["ET"] / h.n_fm / h.dfm)
    Dx = p.alpha_L * vx
    rho = h.rho_fm
    d_s = p.d_s(h)
    Kf = p.Kf if macro_percent != 1 else 0.0
    pu_idx = int(round(p.Pu_Depth / dx)) - 1        # MATLAB I(Pu_Depth/dx:end) = 0  -> python [pu_idx:]
    Z = lambda: np.zeros((N, M))
    C_PP, C_DP, C_mi, C_POP, C_PIP = Z(), Z(), Z(), Z(), Z()
    qa_ma, qa_mi, q_min, q_rec = Z(), Z(), Z(), Z()
    th_POP, th_PIP = Z(), Z()
    P_mic, P_lab, P_sta = Z(), Z(), Z()
    I_ad_ma, I_ad_mi, Ads_flux, Prec_flux, Hyd_flux = Z(), Z(), Z(), Z(), Z()
    I_Pu_ma, I_Pu_mi = Z(), Z()
    PlantP = np.zeros(M)
    Ce_TDP_FM, Ce_PP_FM, Ce_PP_S = np.zeros(M), np.zeros(M), np.zeros(M)
    C_Sr_DP, C_Sr_PP, C_Sr_TP = np.zeros(M), np.zeros(M), np.zeros(M)
    TPflux, TDPflux, Cout_TP, Cout_TDP = np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M)
    Rn_Pu_all = Z()   # liquid-side uptake actually applied (diagnostic of the QUIRK below)

    C_PP[:, 0] = C_init; C_DP[:, 0] = C_init
    q_min[:, 0] = p.qe_mineral0; q_rec[:, 0] = p.qe_recalcitrant0
    P_mic[:, 0] = p.pop_microbial0
    PlantP[0] = p.plantP0
    C_0 = C0 * p.alpha_TDP
    C_s0 = C0 * (1 - p.alpha_TDP)
    fpop = p.alpha_POP / (p.alpha_POP + p.alpha_PIP)
    fpip = p.alpha_PIP / (p.alpha_POP + p.alpha_PIP)
    BaseStorage = h.b_Sr - h.b_bs
    L, Sr, Qd, Qo = hyd["L"], hyd["Sr"], hyd["Qd"], hyd["Qo"]

    for t in range(M - 1):
        C_PP[:, t] *= ET_factor[t]
        C_DP[:, t] *= ET_factor[t]
        q = vx[t] / dx
        s = Dx[t] / dx ** 2
        p_df = -Kf * dt / 2 if macro_percent == 1 else -Kf * dx * dt / 2
        UCa = phi_macro * soil[t] / rho
        UCi = phi_micro * soil[t] / rho
        C_POP[:, t] = np.fmax(0, C_PP[:, t] * fpop)
        C_PIP[:, t] = np.fmax(0, C_PP[:, t] * fpip)
        k_dep = vx[t] * p.lam
        th_POP[:, t] = k_dep * C_POP[:, t] * UCa * dt
        th_PIP[:, t] = k_dep * C_PIP[:, t] * UCa * dt
        # fast sorption, macropores
        Ceq_ma = qa_ma[:, t] / p.kL_ad / (p.Qmax_ad - qa_ma[:, t])
        I_ad_ma[:, t] = p.k_ad * (C_DP[:, t] - Ceq_ma) * dt
        qa_ma[:, t + 1] = np.fmax(0, qa_ma[:, t] - p.k_mi * dt) + th_PIP[:, t] * p.beta_SorbedP + I_ad_ma[:, t] * UCa
        if macro_percent != 1:
            Ceq_mi = qa_mi[:, t] / p.kL_ad / (p.Qmax_ad - qa_mi[:, t])
            I_ad_mi[:, t] = p.k_ad * (C_mi[:, t] - Ceq_mi) * dt
            qa_mi[:, t + 1] = (1 - p.k_mi * dt) * qa_mi[:, t] + I_ad_mi[:, t] * UCi
        Ads_flux[:, t] = I_ad_ma[:, t] * UCa + I_ad_mi[:, t] * UCi
        qa = qa_ma[:, t] + qa_mi[:, t]
        prec = np.minimum(p.k_mi * dt, qa)
        q_min[:, t + 1] = q_min[:, t] + th_PIP[:, t] * p.beta_MineralP + prec
        q_rec[:, t + 1] = q_rec[:, t] + th_PIP[:, t] * p.beta_RecalcitrantP
        Prec_flux[:, t] = prec
        # QUIRK: in the macropore branch the precipitated mass is taken from qa_ma via max(0, qa_ma - k_mi dt) and
        # added to q_min via min(k_mi dt, qa) -- identical only when qa_ma >= k_mi dt.
        # plant uptake
        with np.errstate(divide="ignore", invalid="ignore"):
            Rn_Pu = -p.I_Pu_max / phi_macro / soil[t] / p.Pu_Depth * (C_DP[:, t] - p.C_SRP_Pu_min) / (
                p.Km_Pu + C_DP[:, t] - p.C_SRP_Pu_min) * dt
        Rn_Pu[C_DP[:, t] <= p.C_SRP_Pu_min] = 0.0
        Rn_Pu_all[:, t] = Rn_Pu
        k_Pu = 0.0
        I_Pu_ma[:, t] = -Rn_Pu
        I_Pu_ma[pu_idx:, t] = 0.0
        # QUIRK: Rn_Pu (all depths) is removed from the porewater in the CN solve, but only I_Pu_ma (top Pu_Depth-dx)
        # is credited to the plant pool -> P sink below Pu_Depth is not booked anywhere.
        if macro_percent != 1:
            I_Pu_mi[:, t] = np.fmax(0, C_mi[:, t] - p.C_SRP_Pu_min) * p.I_Pu_max / phi_micro / soil[t] / p.Pu_Depth / (
                p.Km_Pu + C_mi[:, t] - p.C_SRP_Pu_min) * dt
            I_Pu_mi[pu_idx:, t] = 0.0
        # POP pool
        with np.errstate(divide="ignore", invalid="ignore"):
            k_Im = np.fmax(0, p.Um / (p.Km_Im + C_DP[:, t]) * P_mic[:, t] / UCa)
        k_Im = np.nan_to_num(k_Im, nan=0.0)
        I_Im_ma = k_Im * C_DP[:, t] * dt
        if macro_percent != 1:
            with np.errstate(divide="ignore", invalid="ignore"):
                I_Im_mi = np.nan_to_num(np.fmax(0, p.Um * C_mi[:, t] / (p.Km_Im + C_mi[:, t]) * P_mic[:, t] / UCi) * dt)
        else:
            I_Im_mi = np.zeros(N)
        R_Death = -p.k_Death * P_mic[:, t] * dt
        R_Hy = -p.k_Hy * P_lab[:, t] * dt
        with np.errstate(divide="ignore", invalid="ignore"):
            I_Hy_ma = np.nan_to_num(np.fmax(0, -R_Hy / UCa), nan=0.0)
            I_Hy_mi = np.nan_to_num(np.fmax(0, -R_Hy / UCi), nan=0.0) if macro_percent != 1 else np.zeros(N)
        Hyd_flux[:, t] = -R_Hy
        P_mic[:, t + 1] = P_mic[:, t] + I_Im_ma * UCa + I_Im_mi * UCi + R_Death + p.beta_Microbial * th_POP[:, t]
        P_lab[:, t + 1] = P_lab[:, t] - R_Death * (1 - p.f_debris) + R_Hy + p.beta_Labile * th_POP[:, t]
        P_sta[:, t + 1] = P_sta[:, t] - R_Death * p.f_debris + p.beta_Stable * th_POP[:, t]
        if macro_percent != 1:
            C_mi[:, t + 1] = np.fmax(0, C_mi[:, t] + 2 * p_df * (C_mi[:, t] - C_DP[:, t]) - I_ad_mi[:, t] - I_Im_mi
                                     - I_Pu_mi[:, t] + I_Hy_mi)
        # macropore CN coefficients
        if sorption_mode == 1:
            p_filt = -k_dep * dt / 2
            p_Im = -k_Im * dt / 2
            p_Pu = -k_Pu * dt / 2
            pr = -p.k_mi / UCa * p.Qmax_ad / (1 / p.kL_ad + C_DP[:, t]) * dt / 2 + p_df + p_Im + p_Pu
            alpha = 1 + p.Qmax_ad / p.kL_ad / (1 / p.kL_ad + C_DP[:, t]) ** 2 / UCa
            Rn = -2 * p_df * C_mi[:, t] + I_Hy_ma + Rn_Pu
        elif sorption_mode == 2:
            p_ad = -p.k_ad * dt / 2
            p_filt = -k_dep * dt / 2
            p_Im = -k_Im * dt / 2
            p_Pu = -k_Pu * dt / 2
            pr = p_ad + p_df + p_Pu + p_Im
            Rn = -2 * p_ad * Ceq_ma - 2 * p_df * C_mi[:, t] + I_Hy_ma + Rn_Pu
            alpha = np.ones(N)
            # QUIRK: porewater adsorption is CN-averaged (k_ad dt (C^n+C^{n+1})/2) while the solid gains the explicit
            # k_ad dt C^n; filtration likewise (implicit in liquid, explicit theta in solid). Not exactly conservative.
        else:
            pr = 0.0; p_filt = 0.0; Rn = np.zeros(N); alpha = np.ones(N)
            C_DP[:, t] = C_DP[:, t] - I_ad_ma[:, t] - I_Pu_ma[:, t] - I_Im_ma + I_Hy_ma
            if vx[t] != 0:
                C_POP[:, t] = (1 - k_dep * dt) * C_POP[:, t]
                C_PIP[:, t] = (1 - k_dep * dt) * C_PIP[:, t]
                C_PP[:, t] = C_POP[:, t] + C_PIP[:, t]
        s1 = 0.0 if I[t] == 0 else s
        C_DP[:, t + 1] = cn_step(C_DP[:, t], C_0[t], C_0[t + 1], q, s, s1, alpha, pr, Rn, eps)
        C_PP[:, t + 1] = cn_step(C_PP[:, t], C_s0[t], C_s0[t + 1], q, s, s1, 1.0, p_filt, 0.0, eps)
        PlantP[t + 1] = ((1 - p.k_LF * dt) * PlantP[t]
                         + max(0.0, np.sum(I_Pu_mi[:, t]) * h.A_BC * h.dfm * 1000 / N * phi_micro * soil[t])
                         + max(0.0, np.sum(I_Pu_ma[:, t]) * h.A_BC * h.dfm * 1000 / N * phi_macro * soil[t]))
        P_lab[0, t + 1] += p.k_LF * dt * PlantP[t]
        if sorption_mode == 1:
            qa_ma[:, t + 1] = langmuir_q(C_DP[:, t + 1], p.Qmax_ad, p.kL_ad)
        # lower zone
        Ce_TDP_FM[t + 1] = C_DP[-1, t + 1]
        Ce_PP_FM[t + 1] = C_PP[-1, t + 1]
        Ce_PP_S[t + 1] = bedding_filtration(Ce_PP_FM[t + 1], p.lambda_s, d_s)
        Vst = h.n_bs * (Sr[t] + BaseStorage)
        C_Sr_DP[t + 1] = (L[t] * Ce_TDP_FM[t] + Vst * C_Sr_DP[t]) / (L[t] + Vst)
        C_Sr_PP[t + 1] = (L[t] * Ce_PP_S[t] + Vst * C_Sr_PP[t]) / (L[t] + Vst)
        C_Sr_TP[t + 1] = C_Sr_DP[t + 1] + C_Sr_PP[t + 1]
        TPflux[t + 1] = C_Sr_TP[t + 1] * Qd[t + 1] + p.C_overflow * Qo[t + 1]
        TDPflux[t + 1] = C_Sr_DP[t + 1] * Qd[t + 1] + p.C_overflow * p.alpha_TDP * Qo[t + 1]
        if Qd[t] + Qo[t] > 0:
            Cout_TP[t] = TPflux[t] / (Qd[t] + Qo[t])
            Cout_TDP[t] = TDPflux[t] / (Qd[t] + Qo[t])

    qe_ad = qa_ma + qa_mi
    PIP = qe_ad + q_min + q_rec
    POP = P_mic + P_lab + P_sta
    out = dict(C_macro=C_DP, C_micro=C_mi, C_macro_PP=C_PP, C_macro_POP=C_POP, C_macro_PIP=C_PIP, POP=POP,
               M_POP_Microbial=P_mic, M_POP_Labile=P_lab, M_POP_Stable=P_sta, PIP=PIP, qe_ad=qe_ad,
               qe_mineral=q_min, qe_Recalcitrant=q_rec, Plant_P=PlantP, vx=vx, Adsorbed_flux=Ads_flux,
               Precip_flux=Prec_flux, Filtrated_sorbed_flux=th_PIP * p.beta_SorbedP,
               Filtrated_mineral_flux=th_PIP * p.beta_MineralP, hydrolysis_flux=Hyd_flux,
               Ce_TDP_FM=Ce_TDP_FM, Ce_PP_FM=Ce_PP_FM, Ce_PP_Sand=Ce_PP_S, C_Sr_DP=C_Sr_DP, C_Sr_PP=C_Sr_PP,
               C_Sr_TP=C_Sr_TP, TPflux=TPflux, TDPflux=TDPflux, SRPflux=p.beta_SRP * TDPflux,
               TPExflux=C_Sr_TP * hyd["Q_Ex"], TDPExflux=C_Sr_DP * hyd["Q_Ex"], Cout_TP=Cout_TP, Cout_TDP=Cout_TDP,
               theta_POP=th_POP, theta_PIP=th_PIP, Rn_Pu=Rn_Pu_all, I_Pu_macro=I_Pu_ma, I_ad_macro=I_ad_ma)
    out["TPP"] = PIP + POP
    out["Soil_P_content"] = (PIP.mean(axis=0) + POP.mean(axis=0))
    out["Ce_TP_FM"] = Ce_TDP_FM + Ce_PP_FM
    out["Cout_SRP"] = p.beta_SRP * Cout_TDP
    out["TPflux_in"] = C0 * (I * 1000 * h.A_BC + Qo)
    return out
