import sys, zipfile, io, numpy as np
sys.path.insert(0, '.')
R = np.load(sys.argv[1]); Z = zipfile.ZipFile(sys.argv[2])
vars_ = list(R['vars']); grid = R['grid']; times = R['times']
tags = ['tm50984', 'tm32734', 'tm14484', 't3766', 't22016', 't38076', 't40806', 't41076']
rep = []; worst = {v: 0.0 for v in vars_}; worstp = {v: 0.0 for v in vars_}
for k, tag in enumerate(tags):
    lines = Z.read('REF_%s.txt' % tag).decode('latin1').splitlines()
    vi = [i for i, l in enumerate(lines) if l.startswith('Variable:')][0]
    names = lines[vi].split('\t')[1:]
    rows = []
    for l in lines[vi + 1:]:
        p = l.split('\t')[1:]
        try: rows.append([float(x) if x.strip() else np.nan for x in p])
        except ValueError: continue
    A = np.array([r + [np.nan] * (len(names) - len(r)) for r in rows if len(r)])
    C = R['C_%d' % k]; th = R['theta_%d' % k]
    for c in range(0, len(names), 2):
        v = names[c + 1]; z = A[:, c]; y = A[:, c + 1]; ok = np.isfinite(z) & np.isfinite(y)
        z, y = z[ok], y[ok]
        if v not in vars_: continue
        py = np.interp(z, grid, C[:, vars_.index(v)])
        scale = max(np.abs(y).max(), 1e-30)
        e = np.abs(py - y).max() / scale
        m = np.abs(y) > 1e-3 * scale
        ep = (np.abs(py - y)[m] / np.abs(y[m])).max() if m.any() else 0
        worst[v] = max(worst[v], e); worstp[v] = max(worstp[v], ep)
        if k == 0 and c == 0: print('ref points', len(z), 'z range', z.min(), z.max())
    rep.append(tag)
print('%-14s %12s %12s' % ('var', 'max|d|/max', 'max pointrel'))
for v in vars_: print('%-14s %12.3e %12.3e' % (v, worst[v], worstp[v]))
# inventories at final time (python): theta*C or X integrated over inner cells
k = len(times) - 1; C = R['C_%d' % k]; th = R['theta_%d' % k]; dz = abs(grid[2] - grid[1])
print('runtime_s', float(R['runtime_s']), 'nfev', int(R['nfev']))
