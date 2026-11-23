"""Numba kernel for the mixed and advective-diffusive reactor.

Method of lines: the reactor is split into n completely mixed cells (n = 1 gives a
CSTR; n > 1 gives an advective-dispersive plug-flow compartment, i.e. a finite-volume
discretisation / tanks-in-series).  State variables per cell: NO3-N concentration C_N,
PO4-P concentration C_P and sorbed P on the media q.  The time integration is a
fixed-step classical RK4 with `nsub` sub-steps per forcing interval.

Mass-balance bookkeeping: the cumulative outflow mass and cumulative process sink of
every constituent are carried as EXTRA STATES of the same ODE system and are advanced
by the same RK4 stages as the concentrations.  Runge-Kutta methods conserve linear
invariants exactly, so  (inflow) = (outflow) + (sink) + (change in stored mass)  holds to
round-off.  The tests additionally compare the discretisation against independent
analytical solutions and a scipy LSODA solution, so closure is not the only check.

Units are user-defined but must be consistent (the main model uses d, m3, g/m3).
"""
import numpy as np
from numba import njit

# ---- indices into the float parameter vector `p` --------------------------------
P_V, P_DISP, P_QDES = 0, 1, 2
P_KN20, P_VMAXN, P_KN, P_THN = 3, 4, 5, 6
P_KP20, P_KADS, P_KDES, P_QMAX, P_KP, P_KREL, P_THP = 7, 8, 9, 10, 11, 12, 13
P_RHOB_TH = 14                      # kg media per m3 pore water  (= rho_bulk/theta)
P_BIO_A, P_BIO_KV, P_BIO_DF, P_BIO_LF, P_BIO_KL, P_THBIO = 15, 16, 17, 18, 19, 20
P_Q0, P_CN0, P_CP0, P_TREF = 21, 22, 23, 24
P_OUTMODE = 26                      # 0 = zero-gradient outlet (default), 1 = AQUASIM transmission BC (d2C/dx2 = 0)
P_KD = 25                           # linear-isotherm distribution coefficient (m3/kg media... see model)
NP = 27

# ---- output columns -----------------------------------------------------------------
(R_VIN, R_VTR, R_VBY,
 R_NIN, R_NOUT, R_NBY, R_NREM, R_NSTO,
 R_PIN, R_POUT, R_PBY, R_PREM, R_PSORB, R_PSTO, R_PQMASS,
 R_CNL, R_CPL, R_CNF, R_CPF, R_CNM, R_CPM, R_CNL2, R_CPL2) = range(23)
NR = 23


@njit(cache=True, fastmath=False)
def _biofilm_flux_coeff(kv, Df, Lf, kL):
    """Apparent first-order flux coefficient J/C (m/time) of a first-order biofilm with
    external mass-transfer resistance:  1/(1/kL + 1/(sqrt(kv*Df)*tanh(phi))), phi = Lf*sqrt(kv/Df)."""
    if kv <= 0.0 or Df <= 0.0 or Lf <= 0.0:
        return 0.0
    s = np.sqrt(kv * Df)
    phi = Lf * np.sqrt(kv / Df)
    kin = s * np.tanh(phi)
    if kL <= 0.0:
        return kin
    return 1.0 / (1.0 / kL + 1.0 / kin)


@njit(cache=True, fastmath=False)
def _face(y, off, i, n, cin, scheme):
    """Advective concentration at the face between cell i and i+1 (flow towards increasing i)."""
    c0 = y[off + i]
    c1 = y[off + i + 1]
    if scheme == 0:
        return c0
    if scheme == 1:
        return 0.5 * (c0 + c1)
    cm = y[off + i - 1] if i >= 1 else cin           # ghost value upstream of the first cell = inflow conc.
    d1 = c1 - c0
    d0 = c0 - cm
    if d1 * d0 <= 0.0:
        return c0
    r = d0 / d1
    phi = (r + abs(r)) / (1.0 + abs(r))              # van Leer limiter
    return c0 + 0.5 * phi * d1


@njit(cache=True, fastmath=False)
def _rhs(y, dy, n, Qt, T, cinN, cinP, p, scheme, nmode, pmode):
    Vc = p[P_V] / n
    fN = p[P_THN] ** (T - p[P_TREF])
    fP = p[P_THP] ** (T - p[P_TREF])
    E = p[P_DISP] * n * Qt                       # dispersive exchange between cells (m3/time)
    kbio = 0.0
    if p[P_BIO_A] > 0.0:
        kv = p[P_BIO_KV] * p[P_THBIO] ** (T - p[P_TREF])
        kbio = p[P_BIO_A] * _biofilm_flux_coeff(kv, p[P_BIO_DF], p[P_BIO_LF], p[P_BIO_KL])  # 1/time
    remN = 0.0
    remP = 0.0
    remPs = 0.0
    for i in range(n):
        cN = y[i]
        cP = y[n + i]
        q = y[2 * n + i]
        cNp = cN if cN > 0.0 else 0.0
        cPp = cP if cP > 0.0 else 0.0
        # ---- transport (flux in minus flux out, per cell) ---------------------------
        # advective face concentration: scheme 0 upwind, 1 central, 2 TVD (van Leer limiter);
        # explicit dispersive exchange E*(C_i - C_{i+1}) on top of it.
        if i == 0:
            fin_N = Qt * cinN
            fin_P = Qt * cinP
        else:
            fin_N = Qt * _face(y, 0, i - 1, n, cinN, scheme) - E * (cN - y[i - 1])
            fin_P = Qt * _face(y, n, i - 1, n, cinP, scheme) - E * (cP - y[n + i - 1])
        if i == n - 1:
            if p[P_OUTMODE] > 0.5 and n > 1 and E > 0.0:      # transmission BC: linear extrapolation to the outlet face
                fout_N = Qt * (1.5 * cN - 0.5 * y[i - 1]) - E * (cN - y[i - 1])
                fout_P = Qt * (1.5 * cP - 0.5 * y[n + i - 1]) - E * (cP - y[n + i - 1])
            else:
                fout_N = Qt * cN
                fout_P = Qt * cP
        else:
            fout_N = Qt * _face(y, 0, i, n, cinN, scheme) - E * (y[i + 1] - cN)
            fout_P = Qt * _face(y, n, i, n, cinP, scheme) - E * (y[n + i + 1] - cP)
        dN = (fin_N - fout_N) / Vc
        dP = (fin_P - fout_P) / Vc
        # ---- NO3-N removal (g/m3/time) ------------------------------------------------
        rN = 0.0
        if nmode == 1:
            rN = p[P_KN20] * fN * cNp
        elif nmode == 2:
            rN = p[P_VMAXN] * fN * cNp / (p[P_KN] + cNp)
        rN += kbio * cNp
        dN -= rN
        remN += rN * Vc
        # ---- PO4-P removal ---------------------------------------------------------------
        dq = 0.0
        if pmode == 1:                       # first-order, permanent removal
            rP = p[P_KP20] * fP * cPp
            dP -= rP
            remP += rP * Vc
        elif pmode == 2 or pmode == 3 or pmode == 4:   # sorption onto media
            qmax = p[P_QMAX]
            if pmode == 2:                   # kinetic Langmuir
                head = qmax - q
                if head < 0.0:
                    head = 0.0
                dq = fP * (p[P_KADS] * cPp * head - p[P_KDES] * q)
            elif pmode == 3:                 # relaxation to Langmuir equilibrium
                dq = fP * p[P_KREL] * (qmax * cPp / (p[P_KP] + cPp) - q)
            else:                            # relaxation to linear isotherm q_eq = Kd * C
                dq = fP * p[P_KREL] * (p[P_KD] * cPp - q)
            sink = p[P_RHOB_TH] * dq         # g/m3 pore water/time
            dP -= sink
            remP += sink * Vc
            remPs += sink * Vc
        dy[i] = dN
        dy[n + i] = dP
        dy[2 * n + i] = dq
    # cumulative bookkeeping states
    if p[P_OUTMODE] > 0.5 and n > 1 and E > 0.0:
        dy[3 * n + 0] = Qt * (1.5 * y[n - 1] - 0.5 * y[n - 2]) - E * (y[n - 1] - y[n - 2])
        dy[3 * n + 2] = Qt * (1.5 * y[2 * n - 1] - 0.5 * y[2 * n - 2]) - E * (y[2 * n - 1] - y[2 * n - 2])
    else:
        dy[3 * n + 0] = Qt * y[n - 1]        # N outflow of treated stream
        dy[3 * n + 2] = Qt * y[2 * n - 1]    # P outflow of treated stream
    dy[3 * n + 1] = remN
    dy[3 * n + 3] = remP
    dy[3 * n + 4] = remPs


@njit(cache=True, fastmath=False)
def simulate_kernel(Q, T, cinN, cinP, dt, nsub, p, n, scheme, nmode, pmode):
    """Integrate the unit over len(Q) forcing intervals of length dt.  Returns res[k, NR]."""
    nint = Q.shape[0]
    ns = 3 * n + 5
    y = np.zeros(ns)
    for i in range(n):
        y[i] = p[P_CN0]
        y[n + i] = p[P_CP0]
        y[2 * n + i] = p[P_Q0]
    k1 = np.zeros(ns); k2 = np.zeros(ns); k3 = np.zeros(ns); k4 = np.zeros(ns)
    yt = np.zeros(ns)
    res = np.zeros((nint, NR))
    h = dt / nsub
    Vc = p[P_V] / n
    for k in range(nint):
        Qk = Q[k]
        Qt = Qk if Qk < p[P_QDES] else p[P_QDES]
        Qb = Qk - Qt
        # zero the cumulative states at the start of the interval
        for j in range(3 * n, ns):
            y[j] = 0.0
        stoN0 = 0.0; stoP0 = 0.0; mq0 = 0.0
        for i in range(n):
            stoN0 += y[i] * Vc
            stoP0 += y[n + i] * Vc
        for s in range(nsub):
            _rhs(y, k1, n, Qt, T[k], cinN[k], cinP[k], p, scheme, nmode, pmode)
            for j in range(ns):
                yt[j] = y[j] + 0.5 * h * k1[j]
            _rhs(yt, k2, n, Qt, T[k], cinN[k], cinP[k], p, scheme, nmode, pmode)
            for j in range(ns):
                yt[j] = y[j] + 0.5 * h * k2[j]
            _rhs(yt, k3, n, Qt, T[k], cinN[k], cinP[k], p, scheme, nmode, pmode)
            for j in range(ns):
                yt[j] = y[j] + h * k3[j]
            _rhs(yt, k4, n, Qt, T[k], cinN[k], cinP[k], p, scheme, nmode, pmode)
            for j in range(ns):
                y[j] += h * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0
        stoN = 0.0; stoP = 0.0; mq = 0.0
        for i in range(n):
            stoN += y[i] * Vc
            stoP += y[n + i] * Vc
            mq += y[2 * n + i] * p[P_RHOB_TH] * Vc
        res[k, R_VIN] = Qk * dt
        res[k, R_VTR] = Qt * dt
        res[k, R_VBY] = Qb * dt
        res[k, R_NIN] = Qk * cinN[k] * dt
        res[k, R_NOUT] = y[3 * n + 0]
        res[k, R_NBY] = Qb * cinN[k] * dt
        res[k, R_NREM] = y[3 * n + 1]
        res[k, R_NSTO] = stoN - stoN0           # change of dissolved storage in this interval
        res[k, R_PIN] = Qk * cinP[k] * dt
        res[k, R_POUT] = y[3 * n + 2]
        res[k, R_PBY] = Qb * cinP[k] * dt
        res[k, R_PREM] = y[3 * n + 3]
        res[k, R_PSORB] = y[3 * n + 4]
        res[k, R_PSTO] = stoP - stoP0
        res[k, R_PQMASS] = mq                   # sorbed P on media at end of interval (g)
        res[k, R_CNL] = y[n - 1]                # state at END of the interval: last / first cell, mean
        res[k, R_CPL] = y[2 * n - 1]
        res[k, R_CNF] = y[0]
        res[k, R_CPF] = y[n]
        res[k, R_CNL2] = y[n - 2] if n > 1 else y[0]
        res[k, R_CPL2] = y[2 * n - 2] if n > 1 else y[n]
        res[k, R_CNM] = stoN / p[P_V]
        res[k, R_CPM] = stoP / p[P_V]
    return res
