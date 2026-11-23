"""aquasimu: a Python port of the AQUASIM aquatic-system simulator.

This release is the bioreactor-first phase: mixed and advective-diffusive
compartments, a first-order biofilm term, mass-balance closure on every run,
and DREAM(ZS) / PEST++ calibration. Soil, river and lake compartments are not
included.

Comparisons against AQUASIM, and any calibration exercise shipped with the
package, are consistency checks. They are not validation against a field site.
"""

from aquasimu.compartments.advdiff import AdvDiffReactor
from aquasimu.compartments.biofilm import Biofilm
from aquasimu.compartments.mixed import MixedReactor
from aquasimu.compartments.reactor import (
    Result,
    Unit,
    choose_nsub,
    simulate_train,
    simulate_unit,
    train_balance,
)
from aquasimu.core.errors import MassBalanceError
from aquasimu.variables.forcing import Forcing

__all__ = [
    "AdvDiffReactor",
    "Biofilm",
    "Forcing",
    "MassBalanceError",
    "MixedReactor",
    "Result",
    "Unit",
    "choose_nsub",
    "simulate_train",
    "simulate_unit",
    "train_balance",
]

__version__ = "0.1.0"
