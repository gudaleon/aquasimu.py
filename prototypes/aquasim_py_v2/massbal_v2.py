"""Column mass balance over [40806, 41076] d, restarting from the saved state at 40806.
dI/dt = F_top - F_bottom + sum(rates*dz), I = sum over inner cells of (theta*C or X)*dz."""
import sys, numpy as np
sys.path.insert(0, '/workspace/aquasim_py_v2')
from scipy.integrate import solve_ivp
from aquasim_py.sod_model import SedModel, make_jac
m = SedModel(sys.argv[1]); R = np.load(sys.argv[2]); k = 6
C = R['C_%d' % k][:, :m.nv]; th = R['theta_%d' % k]
Y = C[1:-1].copy(); Y[:, ~m.ispart] *= th[1:-1, None]; u0 = Y.ravel()
ts = np.linspace(40806, 41076, 271)
s = solve_ivp(m.rhs, (ts[0], ts[-1]), u0, method='BDF', t_eval=ts, rtol=1e-6, atol=1e-6, jac=make_jac(m))
dz = m.dz; I = []; G = []
for j, t in enumerate(s.t):
    _, d = m.rhs_full(t, s.y[:, j])
    I.append(s.y[:, j].reshape(-1, m.nv).sum(0) * dz)
    G.append(d['flux'][0] - d['flux'][-1] + d['rates'][1:-1].sum(0) * dz)
I = np.array(I); G = np.array(G)
from scipy.integrate import trapezoid
gain = trapezoid(G, s.t, axis=0); dI = I[-1] - I[0]
print('%-14s %12s %12s %12s' % ('var', 'dInventory', 'net_source', 'err/inventory'))
for i, v in enumerate(m.vars):
    print('%-14s %12.4e %12.4e %12.2e' % (v, dI[i], gain[i], abs(dI[i] - gain[i]) / max(abs(I[-1, i]), 1e-30)))
