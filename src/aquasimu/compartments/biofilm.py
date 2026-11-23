"""First-order biofilm term attached to a reactor.

This is not a full biofilm compartment. There is no biomass growth, detachment
or moving boundary. The term adds a first-order NO3-N flux, limited by
diffusion in the film and by an optional external mass-transfer resistance.
"""
from __future__ import annotations

from dataclasses import dataclass

from aquasimu.compartments.reactor import Unit
from aquasimu.processes.kinetics import biofilm_flux_coefficient


@dataclass
class Biofilm:
    """Biofilm flux applied to NO3-N in the bulk cells of a reactor."""

    area: float                 # m2 biofilm per m3 pore water
    kv20: float                 # intrinsic first-order rate at T_ref (1/d)
    Df: float = 1e-4            # diffusivity in the film (m2/d)
    Lf: float = 1e-3            # thickness (m)
    kL: float = 0.0             # external mass-transfer coefficient (m/d); 0 = none
    theta: float = 1.0

    def flux_coefficient(self, temperature: float, t_ref: float = 20.0) -> float:
        """J/C at `temperature`, times nothing yet: multiply by `area` for a volumetric rate."""
        kv = self.kv20 * self.theta ** (temperature - t_ref)
        return biofilm_flux_coefficient(kv, self.Df, self.Lf, self.kL)

    def apply(self, unit: Unit) -> Unit:
        """Return a copy of `unit` with this biofilm term switched on."""
        return unit.with_(
            bio_area=self.area,
            bio_kv20=self.kv20,
            bio_Df=self.Df,
            bio_Lf=self.Lf,
            bio_kL=self.kL,
            bio_theta=self.theta,
        )
