# Where BioWEN-P can reuse aquasimu (not done in this archive)

This folder is the self-contained port, kept here so tests stay independent of
the packaged library in pull request #1. A later step can move the reusable
parts to `src/aquasimu/models/biowenp/` without changing the numerics.

| BioWEN-P piece | aquasimu counterpart (PR #1) | Note |
|---|---|---|
| `transport_AN_v1.py` Crank–Nicolson ADE | `AdvDiffReactor` (advection–dispersion, tanks-in-series) | BioWEN-P uses CN with a Dirichlet top / Neumann bottom, as in the MATLAB; aquasimu uses a method-of-lines kernel with upwind / central / TVD options |
| Kinetic Langmuir `R_ad = k_ad (C − C_eq(q))`, first-order and Michaelis–Menten N, Arrhenius θ | `processes/kinetics.py` (`relax`, Langmuir, MM, first-order, θ^(T−20)) | Same rate laws; eggshell priors (Lee 2022) should stay explicit and not inherit Zhou's sandy-media Qmax |
| `massbalance_AN_v1.py` `gate` / `MassBalanceError` | `checks` / the hard gate inside `simulate_unit` and `simulate_train` | The flow-through augmentation already uses the same rtol 1e-6 convention |
| `simulate_flowthrough` / `simulate_train`, design-flow bypass | `MixedReactor`, `simulate_train`, bypass remixing | The two-tier hierarchy: reduced reactor (PR #1) for design and calibration; BioWEN-P-derived model for multi-year P pools and spent-media checks |

Hydrology unique to BioWEN-P (snow, Thornthwaite PET, IWS, underdrain rating)
has no aquasimu counterpart and should remain in the BioWEN-P module.
