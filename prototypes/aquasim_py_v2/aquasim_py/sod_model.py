"""SEDCOMP port (sedcomp.cpp SEDCOMP::InitCond / SEDCOMP::Delta, highres=FALSE, Qinterf=0)
plus AQCOMP::CalcRates, reduced to an ODE in the inner-cell states.

Unknowns: inner cells j=1..n-2, cell-major: u[(j-1)*nv+i]
  dissolved i : theta_j*C_j   (as AQUASIM stores it)
  particulate : X_j
Algebraic rows of the AQUASIM DAE are solved explicitly:
  Q=0; theta_j = 1-sum X/rho; top: C=topconc, X=inpflux/vsed0, theta=SurfPor;
  bottom: C=given or C_{n-2}, X=X_{n-2}.
"""
import numpy as np
from . import aqufile, formula as F


def _pairs(l, k):
    return [l[i:i + k] for i in range(0, len(l), k)]


class SedModel:
    def __init__(self, path, compname=None, calcname=None):
        t = aqufile.load(path)
        v = aqufile.find_section(t, "VARSYS")
        self.vs = F.VarSys([(v[k], v[k + 1]) for k in range(0, len(v), 2)])
        p = aqufile.find_section(t, "PROCSYS")
        procs = {p[k + 1][1]: p[k + 1] for k in range(0, len(p), 2)}
        c = aqufile.find_section(t, "COMPSYS")
        comps = [c[k + 1] for k in range(0, len(c), 2) if c[k] == "SEDCOMP"]
        sc = comps[0] if compname is None else [x for x in comps if x[1] == compname][0]
        self.sc = sc
        # sc[4] = AQCOMP active variable list; only STATEVARs are unknowns (FORMVARs in
        # the list, e.g. S_HS_plus, are derived and reported as outputs)
        self.vars = [x for x in sc[4] if self.vs.defs[x][0] == "STATEVAR"]
        self.derived = [x for x in sc[4] if self.vs.defs[x][0] != "STATEVAR"]
        nv = self.nv = len(self.vars)
        self.procs = [procs[n] for n in sc[5]]
        self.init = _pairs(sc[7], 3)
        self.topbc = dict(_pairs(sc[8], 2)); self.botbc = dict(_pairs(sc[9], 2))
        self.n = int(sc[10]); self.highres = sc[11] == "TRUE"
        self.Zt, self.Zb = float(sc[16]), float(sc[17])
        self.diss = dict(_pairs(sc[18], 2))
        self.part = {a: (r, b) for a, r, b in _pairs(sc[19], 3)}
        self.surfpor, self.ratepor = sc[20], sc[21]
        n = self.n
        dz = (self.Zb - self.Zt) / (n - 2)
        self.grid = np.r_[self.Zt, self.Zt + (np.arange(1, n - 1) - 0.5) * dz, self.Zb]
        self.dz = abs(self.grid[2] - self.grid[1])
        self.ispart = np.array([x in self.part for x in self.vars])
        cs = aqufile.find_section(t, "CALCSYS")
        calcs = [cs[k + 1] for k in range(0, len(cs), 2)]
        cal = [x for x in calcs if x[7] == "TRUE"][0] if calcname is None else [x for x in calcs if x[1] == calcname][0]
        self.calcnum = float(cal[3]); self.t0 = float(cal[4])
        st = np.array(cal[6], float).reshape(-1, 2)
        self.tout = self.t0 + np.cumsum(np.repeat(st[:, 0], st[:, 1].astype(int)))
        self.calcname = cal[1]
        # state independent part: rho (AQUASIM: varsys->Reset(), evaluated without grid prog vars)
        e = self._env({}, self.t0, 0.0, 0.0)
        self.rho = np.array([float(e.ev(self.part[x][0])) if x in self.part else 0.0 for x in self.vars])
        self.nfev = 0

    def _env(self, state, t, z, theta):
        return self.vs.env(state, {"T": t, "Z": z, "EPSL": theta, "CALC": self.calcnum, "COMP": 0.0})

    # ---------------------------------------------------------------- InitCond
    def initial(self):
        n, nv = self.n, self.nv
        Y = np.zeros((n, nv))
        e = self._env({x: 0.0 for x in self.vars}, self.t0, self.grid, 0.0)
        for zone, var, expr in self.init:
            if int(zone) == 0 and var in self.vars:
                Y[:, self.vars.index(var)] = np.broadcast_to(e.ev(expr), n)
        s = (Y[:, self.ispart] / self.rho[self.ispart]).sum(1)
        theta = 1 - s
        Y[1:-1, ~self.ispart] *= theta[1:-1, None]
        return Y[1:-1].ravel()

    # ---------------------------------------------------------------- full grid
    def expand(self, u, t):
        """inner unknowns -> full-grid conc C (n,nv), theta (n), vsed etc."""
        n, nv, P = self.n, self.nv, self.ispart
        Yi = u.reshape(n - 2, nv)
        th_i = 1 - (Yi[:, P] / self.rho[P]).sum(1)
        C = np.empty((n, nv)); C[1:-1] = Yi; C[1:-1, ~P] /= th_i[:, None]
        e0 = self._env({x: 0.0 for x in self.vars}, t, self.Zt, 0.0)
        thetasurf = float(e0.ev(self.surfpor))
        top = np.zeros(nv)
        for i, x in enumerate(self.vars):
            if x in self.topbc:
                top[i] = float(e0.ev(self.topbc[x]))
        vsed0 = (top[P] / self.rho[P] / (1 - thetasurf)).sum()
        C[0] = top; C[0, P] = top[P] / vsed0
        C[-1] = C[-2]
        for i, x in enumerate(self.vars):
            if x in self.botbc and not P[i]:
                C[-1, i] = float(e0.ev(self.botbc[x]))
        theta = np.r_[thetasurf, th_i, 1 - (C[-1, P] / self.rho[P]).sum()]
        return C, theta, top, vsed0

    def rhs_full(self, t, u):
        self.nfev += 1
        n, nv, P, dz = self.n, self.nv, self.ispart, self.dz
        C, theta, top, vsed0 = self.expand(u, t)
        e = self._env({x: C[:, i] for i, x in enumerate(self.vars)}, t, self.grid, theta)
        rates = np.zeros((n, nv))
        for pr in self.procs:
            r = e.ev(pr[3])
            for var, st in _pairs(pr[4], 2):
                if var in self.vars:
                    rates[:, self.vars.index(var)] += r * e.ev(st)
        D = np.zeros((n, nv))
        for i, x in enumerate(self.vars):
            expr = self.part[x][1] if P[i] else self.diss.get(x)
            if expr is not None:
                D[:, i] = e.ev(expr)
        rp = np.broadcast_to(e.ev(self.ratepor), n)
        vsed = vsed0 + np.r_[0.0, np.cumsum(rp[1:] * dz / (1 - theta[1:]))]
        flux = np.empty((n - 1, nv))          # flux at lower boundary of cell j, j=0..n-2
        # top interface
        flux[0, ~P] = -D[0, ~P] * (C[1, ~P] - C[0, ~P]) * 2 / dz * theta[0]
        flux[0, P] = top[P]
        # inner interfaces j=1..n-3
        Dm = 0.5 * (D[1:-2] + D[2:-1]); tm = 0.5 * (theta[1:-2] + theta[2:-1])
        dC = (C[2:-1] - C[1:-2]) / dz
        vj = vsed[1:-2, None]
        adv = np.where(vj >= 0, vj * C[1:-2], np.where(vsed[2:-1, None] < 0, vsed[2:-1, None] * C[2:-1], 0.0))
        flux[1:-1] = np.where(P, adv - Dm * dC, -Dm * dC * tm[:, None])
        # last interface j=n-2
        j = n - 2
        flux[-1] = np.where(P, vsed[j] * C[j] - D[j + 1] * (C[j + 1] - C[j]) * 2 / dz,
                            -D[j + 1] * (C[j + 1] - C[j]) * 2 / dz * theta[j + 1])
        du = rates[1:-1] + (flux[:-1] - flux[1:]) / dz
        der = {x: np.broadcast_to(e[x], n).copy() for x in self.derived}
        return du.ravel(), dict(C=C, theta=theta, vsed=vsed, flux=flux, rates=rates, derived=der)

    def rhs(self, t, u):
        return self.rhs_full(t, u)[0]

    def sparsity(self, halfcells=1):
        from scipy.sparse import lil_matrix, kron, diags, csc_matrix
        m = self.n - 2
        B = diags([np.ones(m - abs(k)) for k in range(-halfcells, halfcells + 1)], list(range(-halfcells, halfcells + 1)))
        return csc_matrix(kron(B, np.ones((self.nv, self.nv))))


def make_jac(m, floor=1e-10):
    """Banded (block-tridiagonal) finite-difference Jacobian, cells coloured mod 3
    (same band structure AQUASIM gives DASSL; the long-range vsed coupling is dropped,
    exactly like a banded DASSL Jacobian)."""
    from scipy.sparse import csc_matrix
    S = m.sparsity().tocoo()
    rows, cols = S.row, S.col
    nv, N = m.nv, (m.n - 2) * m.nv
    colcell = np.arange(N) // nv
    groups = [np.where((colcell % 3 == c) & (np.arange(N) % nv == i))[0] for c in range(3) for i in range(nv)]
    def jac(t, u):
        f0 = m.rhs(t, u)
        vals = np.zeros(len(rows))
        for g in groups:
            h = 1.5e-8 * np.maximum(np.abs(u[g]), floor)
            up = u.copy(); up[g] += h
            df = m.rhs(t, up) - f0
            hcol = np.zeros(N); hcol[g] = h
            sel = hcol[cols] != 0
            vals[sel] = df[rows[sel]] / hcol[cols[sel]]
        return csc_matrix((vals, (rows, cols)), shape=(N, N))
    return jac
