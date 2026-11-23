"""Completely mixed reactor (one cell)."""
from __future__ import annotations

from dataclasses import dataclass

from aquasimu.compartments.reactor import Unit


@dataclass
class MixedReactor(Unit):
    """CSTR. ``n_cells`` is always 1, whatever was passed in."""

    def __post_init__(self):
        self.n_cells = 1
        super().__post_init__()
