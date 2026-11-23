# bioreactor_model_AN_v1 — woodchip denitrifying bioreactor + eggshell P-sorption reactor (AQUASIM-style, Python)

**Status: research prototype. All calibration runs use SYNTHETIC test data (seeded random numbers), not field data.
Nothing here is validated against a real site.** Comparisons against AQUASIM are *consistency checks of the equations*, not validation.
Every file carries the `_AN_v1` suffix; nothing from earlier work was overwritten.

## 1. Layout

| Path | What |
|---|---|
| `bioreactor_AN_v1/` | Python package (written so it can later be moved to `src/aquasimu/`) |
| `  kernel_AN_v1.py` | numba RK4 method-of-lines kernel (cells, advection/dispersion, kinetics, biofilm, bypass, balance accumulators) |
| `  model_AN_v1.py` | `Unit`, `Forcing`, `Result`, `simulate_unit`, `simulate_train`, mass-balance checks (`MassBalanceError`) |
| `  synthetic_AN_v1.py` | SYNTHETIC forcing, true parameters, noisy observations, treatment-train builder |
| `  calib_common_AN_v1.py` | likelihoods shared by DREAM and PEST++ |
| `  dream_AN_v1.py` | **own numpy implementation of DREAM(ZS)** (not pydream) |
| `  aquasim_io_AN_v1.py` | reader for AQUASIM batch-program listing files |
| `scripts/` | `run_all_AN_v1.sh` (everything), `run_dream_`, `setup_pest_`, `run_pest_`, `compare_`, `run_aquasim_examples_`, `run_benchmarks_` |
| `pest_AN_v1/` | PEST++ control file, template (`.tpl`), instruction (`.ins`) files, `forward_run_AN_v1.py`, run dirs `run_glm/`, `run_ies/` |
| `tests/` | 14 pytest tests |
| `tables_AN_v1/`, `figures_AN_v1/` | CSV tables and PNG figures (see §7) |
| `results_AN_v1/` | DREAM chains (`.npz`) and run logs |
| `benchmarks_AN_v1/aquasim_out/` | numeric *listing outputs* produced by running AQUASIM (see §6 on licence) |

## 2. How to run

```bash
python -m venv --system-site-packages venv && source venv/bin/activate
pip install -r requirements_AN_v1.txt            # numpy scipy pandas matplotlib numba pyemu pytest
export PESTPP_BIN=/path/to/pestpp/bin            # folder with pestpp-glm, pestpp-ies (v5.2.16 used)
python -m pytest -q tests                        # 14 tests
bash scripts/run_all_AN_v1.sh                    # everything below, in order (~20 min on 8 cores)
```
or individually: `scripts/run_dream_AN_v1.py 8000 14000`, `scripts/setup_pest_AN_v1.py`, `scripts/run_pest_AN_v1.py [--ies-only]`,
`scripts/compare_AN_v1.py`, `scripts/run_aquasim_examples_AN_v1.py`, `scripts/run_benchmarks_AN_v1.py`.
PEST++ notes: the model command in the `.pst` is `python forward_run_AN_v1.py`, so the venv's `bin/` must be first on `PATH`
when PEST++ starts workers (done in the scripts). The forward run finds the package by searching parent folders or `BIOREACTOR_AN_V1_ROOT`.

Minimal use of the model:
```python
from bioreactor_AN_v1 import Unit, Forcing, simulate_unit
r = simulate_unit(Unit(n_cells=4, volume=30, q_design=50, n_kinetics="mm", vmax_n=7, K_n=3, theta_n=1.08), forcing)
r.check_balance(1e-6)          # raises MassBalanceError if water/N/P do not close
r.balance_table()              # monthly + annual table;  r.balance()["kgN_retained_per_yr"]
```

## 3. Model

State per cell *i* (n cells, pore volume V, cell volume V_c = V/n), daily forcing Q(t), T(t), C_in,N(t), C_in,P(t):

* Flow routing / **bypass**: Q_t = min(Q, Q_design) enters the unit; Q_b = Q − Q_t bypasses and is re-mixed with the outlet (flow-weighted).
* Mass balance per cell (NO3-N and PO4-P, same form):
  `V_c dC_i/dt = Q_t (C_{i-1} − C_i) [advection] + E (C_{i+1} − 2C_i + C_{i−1}) [dispersion, E = (α_L/L)·n·Q_t] − V_c r_i − (biofilm flux)`
  Inlet is a total-flux (Danckwerts) condition Q_t·C_in. n = 1 is a CSTR; n > 1 is an advective-dispersive plug-flow reactor
  (finite volume). Schemes: `upwind`, `central`, `tvd` (van Leer limiter). Outlet: `zero_gradient`, or AQUASIM's `transmission`
  (d²C/dx² = 0, AQUASIM manual eq. 3.43b; requires `central`).
* **NO3-N kinetics** (user-selectable): first order `r = k20 θ^(T−20) C`; Michaelis–Menten `r = vmax θ^(T−20) C/(K+C)`
  (Arrhenius/θ temperature dependence, T_ref = 20 °C).
* **PO4-P kinetics**: first order; **kinetic Langmuir sorption with finite capacity** `dq/dt = θ^(T−20)(k_ads C (q_max − q) − k_des q)`
  with a sorbed-P store q (g P/kg eggshell) that fills (capacity depletes) — the solid mass is `ρ_b/θ_w · q` per m³ pore water;
  plus `relax` (Langmuir-equilibrium relaxation) and `relax_linear` (linear Kd) modes (used for the AQUASIM soil benchmarks).
* **Optional biofilm** (NO3-N): first-order intrinsic rate k_v in a film of thickness L_f, diffusivity D_f; flux uses the
  effectiveness factor tanh(φ)/φ (φ = L_f √(k_v/D_f)) in series with an external mass-transfer resistance k_L. First-order only.
* **Time integration**: fixed-step RK4 with an automatically chosen number of sub-steps from a stiffness estimate.
  Cumulative outflow, removed and sorbed masses are integrated as extra ODE states, so linear conservation laws hold to round-off
  (RK methods preserve linear invariants).

### Mass balance first
Every simulation computes water, N and P budgets: `in = treated out + bypass + removed + Δstorage` with residuals, relative to the
largest term. `check_balance(tol)` raises `MassBalanceError` (default tol = 1e-6). Because the accumulators are part of the ODE system
this check is partly tautological, so it is cross-checked by (i) analytic CSTR / tanks-in-series steady states, (ii) an
independent scipy `solve_ivp` solution of the same equations, (iii) AQUASIM (§6), and (iv) a test that deliberately corrupts a flux and
requires the error to be raised. Every DREAM/PEST forward run also executes the check (tol 1e-6).
**Results (synthetic year):** max relative residual 7e-17 on the true-parameter run; worst of 200 posterior runs 2.3e-16.
`tables_AN_v1/mass_balance_table_AN_v1.csv` has monthly + annual rows per unit (m³, kg N, kg P).
Synthetic truth: **woodchip removes 25.91 kg N/yr** (of 85.6 kg N in; 10.7 kg N bypassed); **eggshell retains 1.570 kg P/yr** (of 5.16 kg P in).

## 4. Synthetic calibration case (NOT field data)
Seeded one-year daily Q (with 14 storms above the 50 m³/d design flow), T, C_in; train = woodchip (4 tanks, MM, θ) → eggshell (3 tanks,
kinetic Langmuir). Observations: train-outlet NO3-N and PO4-P every 3rd day (121 + 121 = 242), noise σ = a + b·y
(N: 0.15 + 0.05 y; P: 0.01 + 0.05 y). Calibrated: vmax_N 7.0, K_N 3.0, θ_N 1.08, k_ads 0.006, q_max 0.9 (log10 for all but θ). Known: geometry, k_des.

### DREAM(ZS) — own numpy implementation (`dream_AN_v1.py`)
After Vrugt (2016, Environ. Model. Softw. 75:273–316) and ter Braak & Vrugt (2008): 8 chains, archive Z (initialised from the prior,
thinned every K = 10 generations), δ = 1–3 difference pairs, γ = 2.38/√(2δd′) with γ = 1 for 20 % of jumps, λ ~ U(−0.1, 0.1),
ε ~ N(0, 1e-6²), randomised-subspace crossover with 3 CR values adapted by normalised jump distance during the first 10 % of generations,
Metropolis acceptance, reflection at the (uniform) prior bounds, R̂ (Gelman–Rubin) on the second half of the chains; posterior = second half.
**Not implemented:** snooker update, outlier-chain replacement (a simplification of the published algorithm — say so if cited).
Two runs: (1) **fixed-σ Gaussian** (same σᵢ as the PEST weights; 8 000 gen × 8 chains = 64 000 model runs, 154 s, acceptance 7 %);
(2) **heteroscedastic Gaussian** with σ = s0 + s1·y_sim estimated for N and P (9 parameters; 112 000 runs, 485 s, acceptance 5 %).

### PEST++ (v5.2.16, driven with pyEMU 1.7.0)
Template/instruction files + `forward_run_AN_v1.py` in `pest_AN_v1/`; log10 transforms; weights 1/σᵢ (same σ as DREAM run 1);
deliberately wrong start values (vmax 3, K 8, θ 1.04, k_ads 0.002, q_max 2.0). **pestpp-glm** (Gauss–Levenberg–Marquardt, 5 iterations, Φ: 9117 → 197.5
for 242 weighted obs) + linear/FOSM (pyEMU `Schur`, `ErrVar`) uncertainty, identifiability and forecast (kg N, kg P) uncertainty.
**pestpp-ies** (200 realisations, 4 iterations, median Φ 202). The binaries installed fine, so no fallback was needed.

## 5. Results on the synthetic data (all from `tables_AN_v1/`)

**Parameter recovery** (truth → DREAM fixed-σ median [95 % CI], R̂ | PEST-GLM best [linear 95 % CI]):

| par | true | DREAM fixed σ | R̂ | PEST-GLM | DREAM hetero σ |
|---|---|---|---|---|---|
| vmax_N | 7.0 | 7.52 [6.41, 9.11] | 1.011 | 7.49 [6.35, 8.84] | 7.36 [6.47, 8.44] |
| K_N | 3.0 | 3.46 [2.48, 4.88] | 1.012 | 3.42 [2.48, 4.71] | 3.35 [2.56, 4.33] |
| θ_N | 1.08 | 1.0780 [1.072, 1.084] | 1.009 | 1.0776 [1.072, 1.083] | 1.0778 [1.073, 1.083] |
| k_ads | 0.006 | 0.00615 [0.00533, 0.00703] | 1.010 | 0.00619 [0.00539, 0.00712] | 0.00620 [0.00534, 0.00707] |
| q_max | 0.9 | 0.922 [0.836, 1.022] | 1.008 | 0.916 [0.828, 1.014] | 0.905 [0.825, 1.009] |

All true values lie inside every 95 % interval. Max R̂ = 1.012 (run 1), 1.016 (run 2, incl. error-model parameters).
Error-model recovery (run 2): s_N0 0.137 [0.086, 0.207] (true 0.15), s_N1 0.043 (0.05), s_P0 0.0086 (0.01), s_P1 0.051 (0.05) — all within CI.
Strong correlations (identifiability limits): vmax–K 0.985 (DREAM) / 0.99 (PEST-GLM), θ–vmax 0.76, k_ads–q_max −0.96. Only the products/ratios are tightly constrained.

**kg retained per year** (truth 25.913 kg N, 1.570 kg P):

| method | kg N (95 %) | kg P (95 %) | truth inside? |
|---|---|---|---|
| DREAM fixed σ | 26.67 [26.06, 27.31] | 1.618 [1.584, 1.656] | **no** (N and P) |
| DREAM hetero σ | 26.42 [25.91, 26.98] | 1.603 [1.569, 1.636] | yes |
| PEST-GLM linear | 26.69 [26.06, 27.32] | 1.620 [1.584, 1.656] | **no** |
| PEST-IES | 26.69 [26.07, 27.33] | 1.622 [1.584, 1.666] | **no** |

Honest reading: all four methods agree with each other (the PEST-vs-DREAM like-for-like comparison, fixed σ, is essentially identical), but
for this single noise realisation they sit ~3 % (N) / ~3 % (P) above the truth and the fixed-σ intervals just miss it. The intervals
only represent parameter uncertainty given the data realisation; one realisation cannot say whether coverage is nominal (that needs repeated
synthetic experiments, not done). The hetero-σ run, which lets the error model be estimated, covers the truth.
The PEST-GLM linear and IES intervals match DREAM closely, i.e. the problem is close to linear-Gaussian in log-parameters.

## 6. AQUASIM consistency check (not validation)
AQUASIM 2.1g source and `tutor20.zip` example files were downloaded **from Eawag** (https://www.eawag.ch/en/department/siam/software/)
into `/workspace/aquasim_examples_AN` (outside this folder) and the original batch program `aquasimb` was compiled on Linux. The
examples were *run*, so benchmark numbers are real AQUASIM output, not manuals/fabricated. **No Eawag source or example input files are in this
project**; they are referenced by path/URL and re-downloaded by `scripts/run_aquasim_examples_AN_v1.py` (local read-only copy of the source: `/workspace/aquasim_source_original`). Only numeric
*listing outputs* (`benchmarks_AN_v1/aquasim_out/*.txt`) are kept here; because the Eawag licence is unclear these are **excluded from the
zip** — regenerate them with the script. Edits to Eawag inputs (done in the work folder, not here): listing precision raised; grid points 52 → 502 for a refinement run; `proc_b` initial conditions changed
to isolate the Monod step.

Errors as % of scale (`tables_AN_v1/aquasim_benchmark_errors_AN_v1.csv`, `figures_AN_v1/fig_aquasim_consistency_AN_v1.png`):

| case | Python vs AQUASIM |
|---|---|
| B1 batch first-order (proc_a) | RMSE 1e-4 % |
| B2 batch Monod (proc_b mod.) | RMSE 6e-5 % |
| B3 advective-dispersive tracer (advdif_a), n = 100, transmission outlet | RMSE 0.008 %, max 0.018 % (n=50: 0.054 %, n=200: 0.037 %) |
| B3 same with zero-gradient outlet | RMSE 5.7 %, max 9.6 % — **the outlet BC matters**; AQUASIM uses d²C/dx²=0 |
| B4 soil_a/soil_b columns (no sorption, linear, Langmuir ×2 pulses), 502 pts | RMSE 0.004–0.28 %, max 0.01–4 % (largest at the sharp Langmuir shock) |
| B4 same, native 52-pt grid | RMSE 0.09–1.0 %, max up to 11 % at the Langmuir shock |
| B5 biofilm film_a: quasi-steady Monod diffusion–reaction BVP vs AQUASIM's pore-water profile S(z) | S2: RMSE 0.014 / 0.07 / 0.28 %; S1: 1.2 / 2.9 / 1.2 % (days 10 / 25 / 50) |

Caveats: B5 is an equation-level check of the diffusion–reaction profile using AQUASIM's biomass profiles; the package kernel's biofilm term is first-order only
and **full film growth/detachment dynamics (film_a time course, film_b) and the Freundlich-isotherm case were not reproduced**. TVD ("van Leer") was
added because upwind's numerical diffusion gave ~28 % max error versus AQUASIM's high-resolution scheme. I did not find a usable published
numerical dataset for chapter-17 biofilm cases (its DOI page was behind a bot check, not bypassed), so nothing is compared to literature numbers.
Mass balances closed in every benchmark run (relative residuals ≤ 2e-13).

## 7. Figures and tables
Figures: `fig_dream_posterior` (corner plot, truth marked), `fig_dream_convergence` (R̂ trace, chains), `fig_fit_and_forecast` (fit + kg retained histograms),
`fig_pest_vs_dream` (estimates, correlations, identifiability, forecast intervals), `fig_mass_balance`, `fig_aquasim_consistency`.
Tables: `parameter_recovery`, `error_model_recovery`, `forecast_kg_retained`, `pest_vs_dream_parameters` (SDs, identifiability, info gain),
`correlation_{DREAM,PESTGLM,PESTIES}`, `mass_balance_table`, `aquasim_benchmark_errors`, `run_summary_AN_v1.json`.

## 8. Assumptions and limitations
* Daily forcing, constant-parameter, 1-D, isothermal within a day, no dead zones/preferential flow, no gas (N2O/N2) or DOC limitation,
  no clogging/ageing of media, bypass assumed instantaneous and perfectly mixed with outlet. Woodchip N removal is lumped (MM or first order),
  not an explicit carbon/biomass model. Eggshell P-sorption is kinetic Langmuir; no pH/Ca-precipitation pathway.
* Mass balance is exact by construction of the accumulators; it detects coding/flux mistakes and time-step blow-ups but not wrong physics.
* DREAM omits snooker and outlier handling; priors are wide uniform boxes (log-uniform for scale parameters).
* PEST/DREAM comparison is like-for-like only for the fixed-σ run; σᵢ for the fixed run are computed from the *observed* values (σ = a + b·y_obs).
* Coverage of the intervals was checked on one noise realisation only (see §5).
* Biofilm term first-order only; no full film growth/detachment port yet.

## 9. What is still needed from field data (to go beyond synthetic tests)
Daily (or finer) inflow and outflow flow, **plus bypass/overflow records and the design flow**; water temperature in the bed; influent and effluent NO3-N, NH4-N,
total/dissolved PO4-P (flow-proportional or event-sampled, not just grab samples), DO and DOC/BOD of the effluent; reactor geometry, media porosity and bulk density, tracer test(s)
for residence-time distribution (to set n cells and dispersion); eggshell Langmuir batch isotherm/kinetic tests (q_max, k_ads, k_des, pH, Ca) and sorbed-P at start; woodchip age and media
depletion; measurement uncertainties and detection limits (to replace the synthetic σ); and at least one dry and one wet season to test the transferability of θ and the bypass logic.

## 10. Porting into `src/aquasimu/` (gudaleon/aquasimu.py — not pushed)
`kernel_AN_v1.py` (numerics), `model_AN_v1.py` (units/mass balance) and `dream_AN_v1.py` (sampler) are free of project-specific paths and can be moved as modules
(drop the `_AN_v1` suffix on import). `synthetic_AN_v1.py`, `calib_common_AN_v1.py` and `scripts/` are the case-specific parts. Nothing was pushed to any repository.

## References
Vrugt, J.A. (2016) Environ. Model. Softw. 75:273–316. ter Braak & Vrugt (2008) Stat. Comput. 18:435–446. Reichert (1994) AQUASIM, Water Sci. Technol. 30(2):21–30 and
AQUASIM 2.0 user manual/tutorial (Eawag). Doherty (2015) PEST; White et al. (2020) PEST++ v5 (USGS TM 7-C26); White (2018) pyEMU, Environ. Model. Softw. 108:289–304.
