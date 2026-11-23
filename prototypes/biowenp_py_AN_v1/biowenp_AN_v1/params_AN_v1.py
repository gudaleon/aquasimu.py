# Copyright 2020-2023 Bowen Zhou (original MATLAB BioWEN-P); Python port 2026 Alex Neumann / EnviroManagement.
# This file is part of biowenp_py_AN_v1, a Python port of BioWEN-P (Zhou et al. 2023, doi:10.1016/j.scitotenv.2023.162749;
# code: FRDR doi:10.20383/103.0640). Licensed under the GNU Lesser General Public License v3.0 or later
# (see COPYING.LESSER and COPYING). Distributed WITHOUT ANY WARRANTY.
"""Parameter sets for the BioWEN-P port.

Two sets are provided and kept separate on purpose:

* ``HydroParams`` / ``PParams`` defaults = the **calibrated values printed in Table 6 of Zhou et al. (2023)**,
  converted to the internal units of the MATLAB code (10-minute time step, cm for the P column, mg/L, mg/kg).
  The original code read these from ``Inputs/Calibration Results/*/BioC_sbest.out`` files that are NOT in the
  FRDR deposit, so Table 6 is the only public source.
* ``setparameters_prior()`` = the hard-coded values in ``Setparameters.m`` *before* the calibration files are read
  (useful only to reproduce the code literally).

Unit conventions of the original (kept in the faithful port):
  time step dt = 1 means 10 min; hydrologic depths in m per 10 min; water volumes Qo/Qd/Q_Ex in L per 10 min;
  P column in cm; concentrations mg/L; solid pools mg P / kg soil; rho_fm in kg/L.
"""
from dataclasses import dataclass, asdict, replace, field

PO4_MG_PER_UMOL = 95e-3          # the MATLAB code converts umol -> mg with 95 mg/mmol (PO4 basis), kept for fidelity
MIN10_PER_DAY = 24 * 6


@dataclass
class HydroParams:
    # inflow (Table 6A)
    DA: float = 5781.0            # drainage area, m2
    DA_PL: float = 0.0            # parking-lot area, m2
    f: float = 0.62               # fraction of rainfall producing runoff (calibrated)
    la: float = 0.446             # impervious fraction
    A_BC: float = 145.0           # bioretention cell area, m2
    A_PP: float = 530.0           # permeable pavement area, m2
    Tt: float = -0.23             # snowmelt threshold, degC
    DD: float = 0.57e-5           # degree-day factor, m/degC/(10 min)
    Void_max: float = 0.09        # catchment depression storage Sc,max, m
    n_PET: float = 10.0
    R_Pref: float = 0.0           # preferential-flow ratio
    # ponding
    Sp_max: float = 0.15          # m
    # filter media
    hfm_u: float = -0.271         # elevation of top of filter media, m
    dfm: float = 0.45             # filter media depth, m
    n_fm: float = 0.40
    rho_fm: float = 1.6           # kg/L (Setparameters.m; not printed in Table 6)
    Ks: float = 0.0213            # m/(10 min)
    sfc: float = 0.35
    s_h: float = 0.083
    beta: float = 12.5
    # bedding storage
    K_1: float = 189.0            # exfiltration coefficient ke, (10 min)
    u: float = -1.216             # underdrain bottom, m
    b_Sr: float = -1.45           # bottom of IWS water level, m
    n_bs: float = 0.48
    # switches (published run: InfiltrationMethod=1, PETMethod=1, FlowEstimateMethod=2 (neural net); NN not portable,
    # so the port uses FlowEstimateMethod=1 (the code's nonlinear rating curve) — see README gaps)
    PETMethod: int = 1
    FlowEstimateMethod: int = 1
    Inflow_TP_EstimateMethod: int = 3
    C_TP_in_const: float = 4.3    # mg/L used by Inflow_TP_EstimateMethod == 3

    @property
    def A_BMP(self):
        return self.A_BC + self.A_PP

    @property
    def b_bs(self):
        return self.hfm_u - 1.6

    @property
    def Sr_crit(self):
        return self.u - self.b_Sr


@dataclass
class PParams:
    # inflow partitions (Table 6B)
    alpha_POP: float = 0.42
    alpha_PIP: float = 0.18
    alpha_TDP: float = 0.40
    beta_Microbial: float = 0.0
    beta_Labile: float = 0.27
    beta_Stable: float = 0.73
    beta_SorbedP: float = 0.36
    beta_MineralP: float = 0.02
    beta_RecalcitrantP: float = 0.62
    beta_SRP: float = 0.5
    C_overflow: float = 0.33      # mg/L
    # filter media
    lam: float = 0.025            # filtration coefficient, 1/cm
    alpha_L: float = 1.0          # dispersivity, cm
    Kf: float = 0.0               # macro-micro exchange (only if Macro_percent != 1)
    Pu_Depth: float = 15.0        # cm
    I_Pu_max: float = 1e4 * 600 * PO4_MG_PER_UMOL * 3.33e-7   # code units (see Setparameters.m)
    C_SRP_Pu_min: float = PO4_MG_PER_UMOL * 0.25             # mg/L
    Km_Pu: float = PO4_MG_PER_UMOL * 12.7                    # mg/L
    k_LF: float = 0.0
    k_ad: float = 0.15            # (10 min)^-1 (Table 6; the code's manual override is 0.153)
    Qmax_ad: float = 683.5        # mg P/kg
    kL_ad: float = 0.2            # L/mg
    k_mi: float = 0.052 / MIN10_PER_DAY   # zero-order precipitation, mg/kg/(10 min)
    Um: float = 0.0               # immobilisation switched off in published run
    Km_Im: float = 0.0
    f_debris: float = 0.2
    k_Death: float = 0.0
    k_Hy: float = 5.0e-4 / MIN10_PER_DAY  # (10 min)^-1
    # bedding
    lambda_s: float = 0.06        # 1/cm
    # initial solid pools (Section 2.3 / BioC_P_RTM.m)
    qe_mineral0: float = 90.0
    qe_recalcitrant0: float = 190.0
    pop_microbial0: float = 10.0
    plantP0: float = 0.1          # mg

    def d_s(self, h: HydroParams):
        return (h.hfm_u - h.b_Sr - h.dfm) * 100.0


def setparameters_prior():
    """Literal hard-coded values of Setparameters.m before calibration files are loaded (for reference only)."""
    h = HydroParams(f=0.4, Tt=-0.12, DD=0.96e-5, Void_max=0.05, n_PET=5, Ks=138 / 6 / 1000, s_h=0.08, beta=12.1,
                    K_1=173, b_Sr=-1.423, n_bs=0.4)
    p = PParams(alpha_POP=0.4, alpha_PIP=0.2, beta_Labile=0.1, beta_Stable=0.9, beta_SorbedP=0.3125,
                beta_MineralP=0.0625, beta_RecalcitrantP=0.625, C_overflow=0.15, alpha_L=0.1,
                I_Pu_max=1e4 * 600 * PO4_MG_PER_UMOL * 2.68e-7, C_SRP_Pu_min=PO4_MG_PER_UMOL * 0.6,
                Km_Pu=PO4_MG_PER_UMOL * 16, k_ad=0.13, Qmax_ad=400, kL_ad=0.22, k_mi=0.055 / MIN10_PER_DAY,
                k_Hy=0.0003 / MIN10_PER_DAY, lambda_s=0.1384)
    return h, p


def as_table(obj):
    """Return a list of (name, value) rows, used to publish the full input table alongside every run."""
    return list(asdict(obj).items())
