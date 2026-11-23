"""Advective-diffusive compartment (finite volume / tanks in series)."""
from __future__ import annotations

from dataclasses import dataclass

from aquasimu.compartments.reactor import Unit


@dataclass
class AdvDiffReactor(Unit):
    """Plug-flow compartment. Requires ``n_cells >= 2``.

    Dispersion is ``disp_ratio = alpha_L / L``. The face scheme is ``upwind``,
    ``central`` or ``tvd`` (van Leer). The outlet is ``zero_gradient`` or, with
    the central scheme only, ``transmission`` (second derivative of concentration
    is zero at the outlet face).
    """

    def __post_init__(self):
        if self.n_cells < 2:
            raise ValueError(
                "AdvDiffReactor needs n_cells >= 2. A single cell is a MixedReactor."
            )
        super().__post_init__()
