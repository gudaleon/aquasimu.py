"""External forcing series supplied to a compartment."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Forcing:
    """Piecewise-constant forcing, one value per interval of length `dt` (days).

    Concentrations are mass per volume in the same units as the reactor state
    (the bioreactor case uses g/m3, which is mg/L). Flows are volume per time.
    """

    Q: np.ndarray        # inflow m3/d
    T: np.ndarray        # temperature degC
    cinN: np.ndarray     # NO3-N g/m3
    cinP: np.ndarray     # PO4-P g/m3
    dt: float = 1.0
    start: str = "2025-01-01"

    def __post_init__(self):
        self.Q, self.T, self.cinN, self.cinP = [
            np.ascontiguousarray(np.asarray(a, float))
            for a in (self.Q, self.T, self.cinN, self.cinP)
        ]

    def __len__(self):
        return len(self.Q)
