#!/usr/bin/env python3
# Part of biowenp_py_AN_v1 (LGPL-3.0-or-later). Figure: original MATLAB code (run in Octave) vs Python port.
import os, sys
import numpy as np
from scipy.io import loadmat
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
O = loadmat(os.path.join(ROOT, "results_AN_v1/octave_out_90d_AN_v1.mat"), squeeze_me=True)
P = np.load(os.path.join(ROOT, "results_AN_v1/python_faithful_90d_AN_v1.npz"))
t = np.arange(len(O["s"])) / 144
pairs = [("Qd", "Qd", "outflow incl. overflow (L/10 min)"), ("Sr", "Sr", "IWS level Sr (m)"),
         ("TPflux", "P_TPflux", "TP outflow load (mg/10 min)"), ("qe_ad", "P_qe_ad", "adsorbed P at 1 cm (mg/kg)")]
fig, ax = plt.subplots(len(pairs), 2, figsize=(12, 10), gridspec_kw=dict(width_ratios=[3, 1.3]))
for i, (ko, kp, lab) in enumerate(pairs):
    a, b = np.asarray(O[ko], float), P[kp]
    if a.ndim == 2:
        a, b = a[0], b[0]
    ax[i, 0].plot(t, a, lw=2.5, color="tab:orange", label="original MATLAB code in Octave")
    ax[i, 0].plot(t, b, lw=0.8, color="k", ls="--", label="Python port")
    ax[i, 0].set_ylabel(lab, fontsize=8)
    d = np.abs(a - b)
    ax[i, 1].semilogy(t, np.maximum(d, 1e-18), lw=0.6)
    ax[i, 1].set_title("|difference|, max %.1e" % d.max(), fontsize=8)
ax[0, 0].legend(fontsize=8); ax[-1, 0].set_xlabel("day"); ax[-1, 1].set_xlabel("day")
fig.suptitle("BioWEN-P equivalence: original code vs Python port, 90-day SYNTHETIC 10-min forcing, Table 6 parameters\n"
             "(code-consistency check, not validation)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures_AN_v1/equivalence_octave_vs_python_SYNTHETIC_AN_v1.png"), dpi=150)
print("ok")
