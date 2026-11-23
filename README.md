# aquasimu

Python port of the [AQUASIM](https://www.eawag.ch/en/department/siam/software/) aquatic-system simulator (Eawag, Reichert). This repository is the **bioreactor-first phase**: a mixed reactor, an advective-dispersive compartment, a first-order biofilm term, a hard mass-balance check on every run, and calibration with DREAM(ZS) and PEST++.

Soil, river and lake compartments are not in this phase. There is no `LICENSE` file yet; the owner has not chosen one. See `NOTICE.md`.

Comparisons against AQUASIM, and the calibration exercise below, are **consistency checks**. They are not validation against a field site. The calibration data are **synthetic** (seeded random numbers), not measurements.

Eawag source, AQUASIM binaries, `.aqu` examples and AQUASIM listings are not in this repository. The AQUASIM 2.1g distribution has no licence file, and Eawag has not answered a licence query.

## Install

Python 3.10, 3.11 or 3.12.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"          # runtime + pytest + pyEMU + matplotlib
# or, runtime only:
pip install -e .
```

Runtime dependencies are NumPy, SciPy, pandas and Numba. pyEMU is the `calib` extra (PEST++ control files). The PEST++ executables themselves are **not** installed by pip and must not be committed. Put `pestpp-glm` and `pestpp-ies` on `PATH`, or set `PESTPP_BIN` to the directory that contains them (PEST++ 5.2.x, from the [pestpp releases](https://github.com/pestpp/pestpp/releases)).

```bash
python -m pytest                 # reactor tests, a small DREAM run, PEST file setup
```

Tests that launch PEST++ are skipped when the binaries are absent. Tests that read AQUASIM listings are skipped unless `AQUASIMU_AQUASIM_LISTINGS` is set. The DREAM unit test uses 6 chains and 4000 generations of a cheap Gaussian target so CI stays short. It does not run the year-long bioreactor calibration.

## Quick start

A single mixed reactor with first-order nitrate removal. `simulate_unit` raises `MassBalanceError` if water, nitrogen or phosphorus fails to close (default relative tolerance `1e-6`).

```python
import numpy as np
from aquasimu import Forcing, MixedReactor, simulate_unit

days = 400
forcing = Forcing(
    Q=np.full(days, 20.0),       # m3/d
    T=np.full(days, 20.0),       # degC
    cinN=np.full(days, 10.0),    # g/m3 NO3-N
    cinP=np.full(days, 0.5),     # g/m3 PO4-P
    dt=1.0,
)
unit = MixedReactor(volume=20.0, n_kinetics="first_order", k_n20=1.0)
result = simulate_unit(unit, forcing)
result.check_balance(1e-6)
print(result.cout_N[-1])                          # steady state 5 g/m3
print(result.balance()["kgN_retained_per_yr"])
print(result.balance_table())                     # monthly water, N and P
```

An advective-diffusive train (woodchip, then eggshell) with a biofilm term on the first unit:

```python
from aquasimu import AdvDiffReactor, Biofilm, simulate_train, train_balance

wood = Biofilm(area=50, kv20=20, Lf=2e-4, kL=0.5).apply(
    AdvDiffReactor(name="woodchip", n_cells=4, volume=30, q_design=50,
                   n_kinetics="mm", vmax_n=7, K_n=3, theta_n=1.08)
)
egg = AdvDiffReactor(
    name="eggshell", n_cells=3, volume=1.5, q_design=50,
    p_kinetics="langmuir", k_ads=0.006, k_des=0.001, q_max=0.9,
    rho_bulk_over_theta=1100 / 0.45,
)
results = simulate_train([wood, egg], forcing)    # checks each unit and the train
print(train_balance(results))
```

`q_design` is the treated flow. Inflow above it bypasses the unit and is mixed back at the outlet by flow.

## Equations

State in each cell: dissolved NO3-N, dissolved PO4-P, and sorbed P. One cell is a CSTR. Several cells are a finite-volume (tanks-in-series) advective-dispersive compartment. Face schemes are upwind, central, and TVD with a van Leer limiter. The outlet is zero-gradient, or a transmission condition (second derivative of concentration is zero; central scheme only).

Per cell of pore volume \(V_c\), treated flow \(Q_t\) and dispersive exchange \(E = (\alpha_L/L)\, n\, Q_t\):

\[
V_c \frac{dC_i}{dt} = Q_t (C^\*_{i-1/2} - C^\*_{i+1/2}) + E (C_{i+1} - 2 C_i + C_{i-1}) - V_c r_i
\]

The inlet face is the total-flux (Danckwerts) condition \(Q_t C_{in}\). Bypass flow is \(Q_b = Q - Q_t\) with \(Q_t = \min(Q, Q_{design})\), and the reported outlet concentration mixes the treated outflow with the bypassed inflow.

Nitrate removal, with Arrhenius factor \(\theta^{T-20}\):

- first order: \(r = k_{20}\,\theta^{T-20}\, C\)
- Michaelis–Menten: \(r = v_{max}\,\theta^{T-20}\, C / (K + C)\)

Phosphorus:

- first order, permanent removal, same shape as nitrate
- kinetic Langmuir sorption with a finite capacity (sorbed store \(q\), g/kg media):

\[
\frac{dq}{dt} = \theta^{T-20}\big(k_{ads} C (q_{max} - q) - k_{des} q\big)
\]

The dissolved sink is \((\rho_b/\theta_w)\, dq/dt\). Two relaxation modes (Langmuir isotherm, linear \(K_d\)) are also implemented; they are the forms used when checking sorption breakthroughs. They are not a soil compartment.

Optional biofilm on nitrate only, first order. Flux coefficient \(J/C\) for film thickness \(L_f\), diffusivity \(D_f\), intrinsic rate \(k_v\) and external transfer \(k_L\):

\[
\phi = L_f \sqrt{k_v/D_f}, \qquad
J/C = \left(\frac{1}{k_L} + \frac{1}{\sqrt{k_v D_f}\,\tanh\phi}\right)^{-1}
\]

(with the \(k_L\) term omitted when \(k_L = 0\)). The volumetric rate is biofilm area per pore volume times \(J/C\) times concentration. Growth, detachment and a moving film boundary are not in this phase.

Time integration is fixed-step RK4. The number of sub-steps comes from a stiffness estimate. Cumulative outflow and cumulative removal are extra states of the same ODE, so a linear balance is conserved to round-off when the step succeeds.

## Mass balance

Every call to `simulate_unit` and `simulate_train` checks, for water, N and P:

```text
inflow = treated outflow + bypass + removed + change in storage
```

Water has no removal and no storage change. The relative residual is the absolute residual divided by the largest term. Above the tolerance (default `1e-6`) the run raises `MassBalanceError` and does not return. `result.check_balance(tol)` and `result.balance()` expose the same numbers; `result.balance_table()` is the monthly table. DREAM and PEST++ forward runs use this gate, so a balance failure is a failed calibration run.

The check catches dropped fluxes and unstable steps. It does not show that the chosen kinetics are the right ones. The tests also compare CSTR and tanks-in-series steady states with the analytical solutions, and an advective-dispersive case with an independent SciPy integration of the same equations.

`check=False` exists only so a test can corrupt a flux and show that the error is raised. Leave it at the default.

## Synthetic calibration

`aquasimu.examples.bioreactor` builds one year of **synthetic** daily flow, temperature and inflow concentrations (14 storms above a design flow of 50 m3/d), runs woodchip (4 cells, Michaelis–Menten, \(\theta\)) into eggshell (3 cells, kinetic Langmuir), and adds heteroscedastic noise \(\sigma = a + b\, y\) on every third day. Calibrated parameters: `vmax_N`, `K_N`, `theta_N`, `k_ads`, `q_max` (log10 for all but \(\theta\)). Geometry and `k_des` are treated as known.

Nothing in that module is a site measurement.

### DREAM(ZS)

`aquasimu.calib.dream` is a NumPy DREAM(ZS) sampler after Vrugt (2016) and ter Braak & Vrugt (2008): an archive of past states, difference-pair jumps, randomised-subspace crossover with a few CR values adapted during burn-in, reflection at uniform prior bounds, and Gelman–Rubin \(\hat R\). The snooker update and outlier-chain replacement are **not** implemented.

```bash
python scripts/run_dream_synthetic.py 8000 14000
```

That is 8 chains. Run 1 uses the known observation \(\sigma\) (the same weights as PEST). Run 2 also estimates \(\sigma = s_0 + s_1 y_{sim}\) for N and for P. Chains land in `outputs/dream/`. A full pair of runs is on the order of 10^5 reactor simulations.

The Gaussian unit test (`tests/test_dream.py`) is the fast check of the sampler. It is not this script.

### PEST++

`aquasimu.calib.pestpp.build_case` writes a pyEMU control file, template and instruction file for the same synthetic observations (weights \(1/\sigma_i\), log transforms, forecasts of kg N and kg P with weight 0). The model command stored in the control file is `python -m aquasimu.examples.bioreactor.forward_run` (the interpreter that built the case). PEST++ removes a leading `/` on that line, so the command is the interpreter's file name and the run scripts put its directory on `PATH`. The forward run reads `params.dat`, checks the mass balance at tolerance `1e-6`, and writes `output.dat`. A balance failure exits non-zero, which PEST records as a failed model run.

```bash
export PESTPP_BIN=/path/to/pestpp/bin     # or put pestpp-glm and pestpp-ies on PATH
python scripts/run_pest_synthetic.py      # pestpp-glm, then pestpp-ies (200 realisations, 4 iterations)
python scripts/compare_synthetic.py       # tables and figures under outputs/
```

`compare_synthetic.py` needs both the DREAM archives and the PEST++ run directories. It reports parameter recovery, \(\hat R\), kg retained per year, and the monthly mass-balance table for the synthetic truth. Those outputs are consistency checks between two calibrators on one noise realisation.

On that one synthetic realisation the source experiment (same equations, same seeds, same sampler settings) recovered every true parameter inside the 95% interval. Fixed-\(\sigma\) DREAM, PEST-GLM and PEST-IES agreed with each other and sat about 3% high on annual kg N and kg P, so the truth fell just outside those intervals. The heteroscedastic DREAM run covered the truth. One realisation is not a coverage test.

| parameter | true | DREAM fixed-\(\sigma\) median [95%] | PEST-GLM [linear 95%] |
|---|---|---|---|
| vmax_N | 7.0 | 7.52 [6.41, 9.11] | 7.49 [6.35, 8.84] |
| K_N | 3.0 | 3.46 [2.48, 4.88] | 3.42 [2.48, 4.71] |
| theta_N | 1.08 | 1.078 [1.072, 1.084] | 1.078 [1.072, 1.083] |
| k_ads | 0.006 | 0.00615 [0.00533, 0.00703] | 0.00619 [0.00539, 0.00712] |
| q_max | 0.9 | 0.922 [0.836, 1.022] | 0.916 [0.828, 1.014] |

Synthetic truth for that forcing: woodchip retains about 25.9 kg N/yr; eggshell retains about 1.57 kg P/yr. `vmax` and `K` are strongly correlated, and so are `k_ads` and `q_max`.

## AQUASIM listings (optional)

```bash
AQUASIMU_AQUASIM_LISTINGS=/path/to/your/listings python scripts/aquasim_consistency.py
```

The same hook is `pytest` marked `aquasim`: with the variable unset the test is skipped. The directory must be one you supply. Expected names are listed in `aquasimu.benchmarks.aquasim_consistency`. The script compares batch first-order and Michaelis–Menten decay, an advective-diffusive tracer (including the transmission outlet), sorption breakthroughs integrated on the advective-diffusive reactor, and a quasi-steady biofilm diffusion–reaction profile. Differences are consistency checks against another implementation, not validation. Mass balance is still required on every one of those runs.

## Layout

```text
src/aquasimu/
  core/            errors (MassBalanceError)
  variables/       forcing series
  processes/       rate equations (nitrate, Langmuir, biofilm coefficient)
  compartments/    mixed reactor, advective-diffusive reactor, biofilm term, numba kernel
  checks/          mass-balance gate
  calib/           DREAM(ZS), PEST++ case builder
  io/              reader for a user-supplied AQUASIM listing
  benchmarks/      optional listing comparison
  examples/bioreactor/   synthetic train, likelihood, PEST forward run
tests/             pytest
scripts/           synthetic DREAM, synthetic PEST++, comparison, listing check
```

## References

Reichert, P. (1994) AQUASIM — a tool for simulation and data analysis of aquatic systems. Water Sci. Technol. 30(2):21–30. Reichert, P. (1998) AQUASIM 2.0 user manual (Eawag). Vrugt, J.A. (2016) Environ. Model. Softw. 75:273–316. ter Braak, C.J.F. and Vrugt, J.A. (2008) Stat. Comput. 18:435–446. Doherty, J. (2015) PEST. White, J.T. et al. (2020) PEST++ version 5. White, J.T. (2018) pyEMU, Environ. Model. Softw. 108:289–304.
