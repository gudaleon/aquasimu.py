"""DREAM(ZS) sanity test on a known 4-d correlated Gaussian target (recovers mean/cov, R-hat < 1.1)."""
import numpy as np
from bioreactor_AN_v1.dream_AN_v1 import dream_zs, gelman_rubin

MU = np.array([1.0, -2.0, 0.5, 3.0])
A = np.array([[1, .6, 0, 0], [.6, 1, .3, 0], [0, .3, 1, -.5], [0, 0, -.5, 1]]) * 0.25
PREC = np.linalg.inv(A)


def lp(x):
    r = x - MU
    return -0.5 * r @ PREC @ r


def test_dream_recovers_gaussian():
    out = dream_zs(lp, [-10] * 4, [10] * 4, n_chains=6, n_gen=4000, n_proc=1, seed=3, verbose=False)
    s = out["chains"][2000:].reshape(-1, 4)
    assert np.allclose(s.mean(0), MU, atol=0.05)
    assert np.allclose(np.cov(s.T), A, atol=0.07)
    assert gelman_rubin(out["chains"][2000:]).max() < 1.1
