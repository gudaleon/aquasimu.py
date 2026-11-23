"""Reactor compartments implemented in this phase.

Mixed and advective-diffusive reactors, plus a first-order biofilm term.
Soil, river and lake compartments are not part of this phase.
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
)

__all__ = [
    "AdvDiffReactor",
    "Biofilm",
    "MixedReactor",
    "Result",
    "Unit",
    "choose_nsub",
    "simulate_train",
    "simulate_unit",
]
