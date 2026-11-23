import sys, time, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from scipy.integrate import solve_ivp
from aquasim_py.sod_model import SedModel, make_jac
path = sys.argv[1]; out = sys.argv[2]; rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
tends = [-50984, -32734, -14484, 3766, 22016, 38076, 40806, 41076]
m = SedModel(path)
u0 = m.initial()
t0 = time.time(); r0 = m.rhs(m.t0, u0); print('rhs ms', (time.time()-t0)*1000, 'neq', u0.size, flush=True)
T = time.time()
sol = solve_ivp(m.rhs, (m.t0, tends[-1]), u0, method='BDF', t_eval=tends, rtol=rtol, atol=1e-6,
                jac=make_jac(m), first_step=1e-6)
el = time.time() - T
print(sol.message, 'nfev', m.nfev, 'njev', sol.njev, 'nlu', sol.nlu, 'runtime_s', el, flush=True)
res = {'vars': m.vars + m.derived, 'grid': m.grid, 'times': sol.t, 'runtime_s': el, 'nfev': m.nfev}
for k, t in enumerate(sol.t):
    _, d = m.rhs_full(t, sol.y[:, k])
    res['C_%d' % k] = np.c_[d['C'], np.array([d['derived'][x] for x in m.derived]).T]; res['theta_%d' % k] = d['theta']; res['vsed_%d' % k] = d['vsed']
    res['flux_%d' % k] = d['flux']; res['rates_%d' % k] = d['rates']
np.savez_compressed(out, **{k: np.asarray(v) for k, v in res.items()})
