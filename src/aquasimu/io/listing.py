"""Parser for AQUASIM batch listing files produced by `aquasimb -l`.

This module reads a user-supplied text listing. It does not contain, and this
package does not ship, Eawag source, binaries, `.aqu` files, or listings.
"""
import numpy as np


def read_listing(path):
    """Return (legends, curves) where curves[i] = (x[], y[]) for the i-th (argument, value) column pair."""
    rows, legend = [], None
    with open(path) as f:
        started = False
        for line in f:
            if line.startswith("Legend:"):
                legend = line.rstrip("\n").split("\t")[1:]
                started = True
                continue
            if started:
                if not line.strip():
                    if rows:
                        break
                    continue
                f_ = line.rstrip("\n").split("\t")[1:]          # first field is empty
                rows.append([float(v) if v.strip() else np.nan for v in f_])
    w = max(len(r) for r in rows)
    rows = [r + [np.nan] * (w - len(r)) for r in rows]
    a = np.array(rows)
    curves = []
    for i in range(0, a.shape[1] - 1, 2):
        ok = ~np.isnan(a[:, i]) & ~np.isnan(a[:, i + 1])
        curves.append((a[ok, i], a[ok, i + 1]))
    legend = legend[1::2] if legend and len(legend) >= 2 * len(curves) else legend
    return legend, curves
