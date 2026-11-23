"""Model command used by PEST / PEST++ for the SYNTHETIC bioreactor case.

PEST writes ``params.dat`` from ``params.tpl``, runs this module, and reads
``output.dat`` through ``output.ins``. The mass balance is checked inside the
model: a balance that does not close aborts the process, which PEST records as
a failed run.

Run from the case directory::

    python -m aquasimu.examples.bioreactor.forward_run
"""
from __future__ import annotations

import numpy as np

from aquasimu.examples.bioreactor import synthetic as case
from aquasimu.examples.bioreactor.likelihood import D, FORCING


def main() -> None:
    values = {}
    with open("params.dat") as handle:
        for line in handle:
            key, value = line.split()
            values[key] = float(value)
    parameters = np.array([values[name.lower()] for name in case.NAMES])
    out = case.run_model(parameters, FORCING, tol=1e-6)
    with open("output.dat", "w") as handle:
        for day in D:
            handle.write(f"on_{day:03d} {out['cN'][day]:.10e}\n")
        for day in D:
            handle.write(f"op_{day:03d} {out['cP'][day]:.10e}\n")
        handle.write(f"fkgn {out['kgN']:.10e}\n")
        handle.write(f"fkgp {out['kgP']:.10e}\n")


if __name__ == "__main__":
    main()
