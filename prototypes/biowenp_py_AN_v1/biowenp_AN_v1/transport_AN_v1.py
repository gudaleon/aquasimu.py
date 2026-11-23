# Copyright 2020-2023 Bowen Zhou (original MATLAB BioWEN-P); Python port 2026 Alex Neumann / EnviroManagement.
# Part of biowenp_py_AN_v1, a Python port of BioWEN-P (Zhou et al. 2023, doi:10.1016/j.scitotenv.2023.162749;
# code: FRDR doi:10.20383/103.0640). GNU LGPL v3.0 or later (see COPYING.LESSER, COPYING). NO WARRANTY.
"""1-D advection-dispersion(-reaction) Crank-Nicolson step, ported from the matrix assembly in ``BioC_P_RTM.m``.

The original builds dense N x N matrices and solves with ``A\\B``; the matrix is tridiagonal, so here it is
solved with ``scipy.linalg.solve_banded`` (same algebra, O(N)).  Notation as in the MATLAB code:

  q  = v_x/dx            advection number (v_x already in cm per time step)
  s  = D_x/dx^2          dispersion number (D_x = alpha_L v_x, Table 4 eq. T4.10)
  s1 = s, or 0 when there is no infiltration (upper boundary switched off)
  alpha = retardation-like diagonal factor (1 for kinetic sorption; R for equilibrium sorption)
  p  = linear reaction coefficient (already multiplied by dt/2, negative for a sink)
  Rn = explicit source term per step (mg/L)
  Upper boundary: Dirichlet (inflow concentration, mean of t and t+1 via the C0(t)+C0(t+1) term);
  lower boundary: zero-gradient (Neumann) via a mirrored ghost node.
"""
import numpy as np
from scipy.linalg import solve_banded


def cn_matrices(C, C0_t, C0_t1, q, s, s1, alpha, p, Rn, eps=0.5):
    """Return (ab, B): banded LHS (3 x N, for solve_banded with l=u=1) and RHS, exactly as BioC_P_RTM.m."""
    N = C.size
    alpha = np.broadcast_to(np.asarray(alpha, float), (N,))
    p = np.broadcast_to(np.asarray(p, float), (N,))
    Rn = np.broadcast_to(np.asarray(Rn, float), (N,))
    diag = alpha + 2 * eps * s - p
    lower = np.full(N - 1, -eps * s - q / 4)
    upper = np.full(N - 1, -eps * s + q / 4)
    # row 1 (upper boundary)
    diag = diag.copy()
    diag[0] = alpha[0] + 2 * eps * s1 - p[0]
    upper[0] = -eps * s1 + q / 4
    # row N (Neumann)
    lower[-1] = -2 * eps * s
    B = np.empty(N)
    B[1:-1] = (((1 - eps) * s + q / 4) * C[:-2] + (alpha[1:-1] - 2 * (1 - eps) * s + p[1:-1]) * C[1:-1]
               + ((1 - eps) * s - q / 4) * C[2:] + Rn[1:-1])
    B[0] = (((1 - eps) * s1 + q / 4) * (C0_t + C0_t1) + (alpha[0] - 2 * (1 - eps) * s1 + p[0]) * C[0]
            + ((1 - eps) * s1 - q / 4) * C[1] + Rn[0])
    B[-1] = 2 * (1 - eps) * s * C[-2] + (alpha[-1] - 2 * (1 - eps) * s + p[-1]) * C[-1] + Rn[-1]
    ab = np.zeros((3, N))
    ab[0, 1:] = upper
    ab[1, :] = diag
    ab[2, :-1] = lower
    return ab, B


def cn_step(C, C0_t, C0_t1, q, s, s1, alpha=1.0, p=0.0, Rn=0.0, eps=0.5):
    ab, B = cn_matrices(C, C0_t, C0_t1, q, s, s1, alpha, p, Rn, eps)
    return solve_banded((1, 1), ab, B)


def dense_matrix(ab):
    """Expand the banded form to a dense matrix (used by tests to compare with the MATLAB layout)."""
    N = ab.shape[1]
    A = np.diag(ab[1]) + np.diag(ab[0, 1:], 1) + np.diag(ab[2, :-1], -1)
    return A
