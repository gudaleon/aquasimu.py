"""Likelihoods shared by DREAM and PEST++ on the SYNTHETIC bioreactor case.

Importing this module builds the synthetic forcing and the synthetic
observations (one model run). Nothing in those arrays is a measurement.
"""
from __future__ import annotations

import numpy as np

from aquasimu.examples.bioreactor import synthetic as case

FORCING = case.make_forcing()
OBS = case.make_observations(FORCING)
SIG_N, SIG_P = case.obs_sigma(OBS)
D = OBS["days"]

# Heteroscedastic error model for the second DREAM run: sigma = s0 + s1 * y_sim.
# s0 is sampled in log10; s1 is linear.
ERR_NAMES = ["sN0", "sN1", "sP0", "sP1"]
ERR_LO = np.array([np.log10(0.01), 0.0, np.log10(0.001), 0.0])
ERR_HI = np.array([np.log10(1.0), 0.3, np.log10(0.1), 0.3])
ERR_TRUE = np.array([
    np.log10(case.NOISE_N[0]), case.NOISE_N[1],
    np.log10(case.NOISE_P[0]), case.NOISE_P[1],
])


def forward(z):
    """Sampling-space parameters (length 5) -> outlet observations and kg retained.

    The mass balance is checked inside the model and raises MassBalanceError
    when it does not close.
    """
    out = case.run_model(case.to_physical(z), FORCING)
    return out["cN"][D], out["cP"][D], out["kgN"], out["kgP"]


def logpost_fixed_sigma(z):
    """Gaussian log posterior with the known observation standard deviations. Uniform prior."""
    z = np.asarray(z)
    if np.any(z < case.Z_LO) or np.any(z > case.Z_HI):
        return -np.inf
    y_n, y_p, _, _ = forward(z)
    return -0.5 * (
        np.sum(((OBS["obsN"] - y_n) / SIG_N) ** 2) + np.sum(((OBS["obsP"] - y_p) / SIG_P) ** 2)
    )


def logpost_hetero(zz):
    """Heteroscedastic Gaussian log posterior. sigma_i = s0 + s1 * y_sim,i is estimated."""
    zz = np.asarray(zz)
    z, error = zz[:5], zz[5:]
    lo = np.r_[case.Z_LO, ERR_LO]
    hi = np.r_[case.Z_HI, ERR_HI]
    if np.any(zz < lo) or np.any(zz > hi):
        return -np.inf
    y_n, y_p, _, _ = forward(z)
    sigma_n = 10 ** error[0] + error[1] * y_n
    sigma_p = 10 ** error[2] + error[3] * y_p
    return -(np.sum(np.log(sigma_n)) + np.sum(np.log(sigma_p))) - 0.5 * (
        np.sum(((OBS["obsN"] - y_n) / sigma_n) ** 2) + np.sum(((OBS["obsP"] - y_p) / sigma_p) ** 2)
    )
