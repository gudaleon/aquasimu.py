# Copyright 2020-2023 Bowen Zhou (original MATLAB BioWEN-P); Python port 2026 Alex Neumann / EnviroManagement.
# Part of biowenp_py_AN_v1, a Python port of BioWEN-P (Zhou et al. 2023, doi:10.1016/j.scitotenv.2023.162749;
# code: FRDR doi:10.20383/103.0640). GNU LGPL v3.0 or later (see COPYING.LESSER, COPYING). NO WARRANTY.
"""Mass-balance checks.

* ``water_balance_original`` / ``p_balance_original``: ports of PostsimulationAnalysis/Waterbalancechecking.m and
  Pmassbalancechecking.m (the original authors' own diagnostics, units L and mg).
* ``water_balance_full``: an extended water ledger that also books ponding storage and catchment depression
  retention, so the residual is the integration error alone.
* ``MassBalanceError`` and ``gate``: the explicit closure gate used by every run of the augmented model.
"""
import numpy as np


class MassBalanceError(RuntimeError):
    pass


def gate(residuals: dict, rtol: float):
    """Raise MassBalanceError if any |relative residual| exceeds rtol. residuals: name -> (residual, scale)."""
    bad = {k: r / s for k, (r, s) in residuals.items() if s > 0 and abs(r / s) > rtol}
    if bad:
        raise MassBalanceError("mass balance not closed (rtol=%g): %s" % (rtol, bad))
    return True


def water_balance_original(hyd, h):
    """Waterbalancechecking.m (L): [in, out(Qd incl. overflow), overflow, underdrain, exfil, ET, storage, residual]."""
    tin = np.nansum(hyd["Q_in"]) * 1000
    tout = np.nansum(hyd["Qd"])
    tov = np.nansum(hyd["Qo"])
    tex = np.nansum(hyd["Q_Ex"])
    tet = np.nansum(hyd["ET"]) * h.A_BC * 1000
    sto = ((hyd["s"][-1] - hyd["s"][0]) * h.A_BC * h.dfm * h.n_fm * 1000
           + (hyd["Sr"][-1] - hyd["Sr"][0]) * h.A_BC * h.n_bs * 1000)
    res = tin - tout - tex - tet - sto
    return dict(inflow_L=tin, outflow_L=tout, overflow_L=tov, underdrain_L=tout - tov, exfiltration_L=tex,
                ET_L=tet, storage_change_L=sto, residual_L=res, residual_rel=res / tin if tin else 0.0)


def water_balance_full(hyd, h):
    """Extended ledger (L) on the BMP inflow volume actually delivered (Q_in_d, after catchment retention),
    including ponding storage. Remaining residual = clamps (max(0,.)) and the s>1 overshoot of the original."""
    A = h.A_BC * 1000
    tin = np.sum(hyd["Q_in_d"]) * A
    tout = np.sum(hyd["Qd"])
    tex = np.sum(hyd["Q_Ex"])
    tet = np.sum(hyd["ET"]) * A
    sto = ((hyd["s"][-1] - hyd["s"][0]) * h.dfm * h.n_fm * A + (hyd["Sr"][-1] - hyd["Sr"][0]) * h.n_bs * A
           + (hyd["Sp"][-1] - hyd["Sp"][0]) * A)
    retained = np.sum(hyd["Retain"]) * A
    res = tin - tout - tex - tet - sto
    return dict(inflow_after_catchment_L=tin, catchment_retention_L=retained, outflow_L=tout, exfiltration_L=tex,
                ET_L=tet, storage_change_L=sto, residual_L=res, residual_rel=res / tin if tin else 0.0)


def p_balance_original(hyd, po, h, p):
    """Pmassbalancechecking.m (mg). Column names follow the comments in the original script."""
    A = h.A_BC * 1000
    rho = h.rho_fm
    m = {}
    m["inflow_mg"] = np.nansum(hyd["C_TP_in"] * hyd["I"]) * A
    m["outflow_mg"] = np.sum(po["TPflux"])
    m["leak_to_deeper_mg"] = A * np.nansum(po["Ce_TP_FM"] * hyd["L"])
    m["exfil_plus_deep_storage_mg"] = m["leak_to_deeper_mg"] - m["outflow_mg"] + np.nansum(p.C_overflow * hyd["Qo"])
    m["exfiltration_mg"] = np.sum(po["TPExflux"])
    m["soil_accum_mg"] = (po["TPP"][:, -1].mean() - po["TPP"][:, 0].mean()) * h.dfm * A * rho
    m["plant_accum_mg"] = po["Plant_P"][-1] - po["Plant_P"][0]
    N = po["C_macro"].shape[0]
    Mpw = (np.sum(po["C_macro"][:, -1] + po["C_macro_POP"][:, -1] + po["C_macro_PIP"][:, -1])
           * h.dfm * A * rho * hyd["s"][-1] / N)
    m["porewater_mg_as_coded"] = Mpw
    m["closure_residual_mg"] = (m["inflow_mg"] - m["leak_to_deeper_mg"] - m["soil_accum_mg"]
                                - m["plant_accum_mg"] - Mpw)
    m["closure_residual_rel"] = m["closure_residual_mg"] / m["inflow_mg"] if m["inflow_mg"] else 0.0
    # diagnostic of the uptake booking quirk: porewater uptake below Pu_Depth that is not credited to the plant pool
    pu_idx = int(round(p.Pu_Depth)) - 1
    vol = h.A_BC * h.dfm * 1000 / N * h.n_fm * hyd["s"][None, :]
    m["unbooked_uptake_mg"] = float(np.sum(-po["Rn_Pu"][pu_idx:, :] * vol))
    return m
