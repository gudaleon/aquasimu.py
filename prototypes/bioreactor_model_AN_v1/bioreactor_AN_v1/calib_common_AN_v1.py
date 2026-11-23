"""Shared calibration set-up used by DREAM and PEST++ (AN_v1): synthetic data, forward model
wrapper on the sampling (log10) space, and the two likelihoods.  ALL DATA ARE SYNTHETIC."""
import numpy as np
from . import synthetic_AN_v1 as S

FORCING = S.make_forcing()
OBS = S.make_observations(FORCING)
SIG_N, SIG_P = S.obs_sigma(OBS)
D = OBS["days"]

# extra error-model parameters for the heteroscedastic DREAM run: sigma = s0 + s1 * y_sim
ERR_NAMES = ["sN0", "sN1", "sP0", "sP1"]
ERR_LO = np.array([np.log10(0.01), 0.0, np.log10(0.001), 0.0])      # s0 in log10, s1 linear
ERR_HI = np.array([np.log10(1.0), 0.3, np.log10(0.1), 0.3])
ERR_TRUE = np.array([np.log10(S.NOISE_N[0]), S.NOISE_N[1], np.log10(S.NOISE_P[0]), S.NOISE_P[1]])


def forward(z):
    """z = sampling-space parameter vector (len 5) -> (yN_obs, yP_obs, kgN, kgP).  Mass balance is
    checked inside (raises MassBalanceError if it does not close)."""
    out = S.run_model(S.to_physical(z), FORCING)
    return out["cN"][D], out["cP"][D], out["kgN"], out["kgP"]


def logpost_fixed_sigma(z):
    """Gaussian likelihood with the (known) observation sigmas used as PEST weights; uniform prior."""
    z = np.asarray(z)
    if np.any(z < S.Z_LO) or np.any(z > S.Z_HI):
        return -np.inf
    yN, yP, _, _ = forward(z)
    return -0.5 * (np.sum(((OBS["obsN"] - yN) / SIG_N) ** 2) + np.sum(((OBS["obsP"] - yP) / SIG_P) ** 2))


def logpost_hetero(zz):
    """Heteroscedastic Gaussian likelihood; sigma_i = s0 + s1 * y_sim,i estimated with the model parameters."""
    zz = np.asarray(zz); z, e = zz[:5], zz[5:]
    lo = np.r_[S.Z_LO, ERR_LO]; hi = np.r_[S.Z_HI, ERR_HI]
    if np.any(zz < lo) or np.any(zz > hi):
        return -np.inf
    yN, yP, _, _ = forward(z)
    sN = 10 ** e[0] + e[1] * yN; sP = 10 ** e[2] + e[3] * yP
    return -(np.sum(np.log(sN)) + np.sum(np.log(sP))) - 0.5 * (
        np.sum(((OBS["obsN"] - yN) / sN) ** 2) + np.sum(((OBS["obsP"] - yP) / sP) ** 2))
