# biowenp_py_AN_v1 — Python port of BioWEN-P (Zhou et al. 2023) + flow-through eggshell/woodchip augmentation

**Status:** research prototype, 23 Nov 2026. LGPL-3.0-or-later (see `NOTICE_AN_v1.md`, `COPYING.LESSER`).
Not pushed to any repository. **All runs use SYNTHETIC forcing.** Every comparison below is a *consistency
check* (code vs code, code vs paper); nothing here is validated against field data.

Original: Zhou, Shafii, Parsons, Passeport, Rezanezhad, Lisogorsky, Van Cappellen (2023), *Sci. Total Environ.* 876,
162749, doi:10.1016/j.scitotenv.2023.162749; MATLAB code FRDR doi:10.20383/103.0640 (read-only snapshot in
`/workspace/parsons_bioretention_AN_v1/`, never modified).

## 1. Layout

| Path | Ported from / purpose |
|---|---|
| `biowenp_AN_v1/params_AN_v1.py` | `Setparameters.m`; defaults = calibrated values of paper Table 6 (the calibration `.out` files are not in FRDR) |
| `biowenp_AN_v1/hydrology_AN_v1.py` | hydrologic loop of `ModelMainloop.m`, `InflowEstimator.m`, `Discharge.m` (snow, runoff, catchment storage, infiltration, ET, leakage, ponding/overflow, IWS, underdrain, exfiltration) |
| `biowenp_AN_v1/transport_AN_v1.py` | Crank–Nicolson ADE matrix assembly of `BioC_P_RTM.m` (banded solver) |
| `biowenp_AN_v1/p_pools_AN_v1.py` | `BioC_P_RTM.m`: rate laws of Table 4 as functions + statement-by-statement time loop (`run_p_faithful`) |
| `biowenp_AN_v1/massbalance_AN_v1.py` | `Waterbalancechecking.m`, `Pmassbalancechecking.m`, extended water ledger, `MassBalanceError`/`gate` |
| `biowenp_AN_v1/flowthrough_AN_v1.py` | **augmentation**: saturated flow-through units, bypass, N, kinetic Langmuir with capacity depletion, gate |
| `biowenp_AN_v1/priors_AN_v1.py` | **augmentation**: eggshell priors (Lee et al. 2022), woodchip priors (Addy et al. 2016, Schipper et al. 2010) |
| `biowenp_AN_v1/synthetic_AN_v1.py` | SYNTHETIC forcing generators (seeded) |
| `scripts/octave_AN_v1/run_biowenp_octave_AN_v1.m` | runs the ORIGINAL code in GNU Octave 9.4 on the same forcing |
| `scripts/run_equivalence_octave_AN_v1.py` | Octave-vs-Python array comparison |
| `scripts/run_crosscheck_soilp_AN_v1.py` | port vs the Soil P thread's Table 4/6 reimplementation |
| `scripts/run_demo_AN_v1.py`, `plot_equivalence_AN_v1.py`, `run_all_AN_v1.sh` | demo and figures |
| `tests/test_biowenp_AN_v1.py` | 16 pytest tests |

Style follows `/workspace/bioreactor_model_AN_v1` (aquasimu.py prototype): dataclass units (`FlowThroughUnit`),
forcing objects, `simulate_flowthrough` / `simulate_train`, `Result.check_balance(rtol)`, `Result.balance_table()`,
same units (d, m³, g/m³) and θ^(T−20) convention, so the flow-through module can later move to `src/aquasimu/`.

## 2. Run
```bash
pip install -r requirements_AN_v1.txt        # numpy scipy pandas matplotlib pytest (+ octave for equivalence)
python3 -m pytest -q tests                   # 16 passed
bash scripts/run_all_AN_v1.sh                # tests, Octave equivalence, cross-check, demo
```

## 3. Equivalence with the original code (consistency check)
90-day SYNTHETIC 10-min forcing (12,960 steps; storms, a snow event, an overflow storm), Table 6 parameters,
published switches except the rating curve (see gaps). The **unmodified** MATLAB code was run in Octave and compared
with the port for 32 output arrays (hydrology and all P pools, depth × time):
worst relative difference **6.7 × 10⁻¹⁰** (the tiny precipitation-flux array); all others ≤ 5 × 10⁻¹¹
(`tables_AN_v1/equivalence_octave_vs_python_90d_AN_v1.csv`, figure `figures_AN_v1/equivalence_octave_vs_python_SYNTHETIC_AN_v1.png`).
Octave 10.9 s, Python 3.0 s.

Cross-check against the independent reimplementation from the paper's Table 4/6 (Soil P thread,
`03_parsons_zhou2023`, no-uptake branch, saturated 0.1 m/d, 1 yr): media accumulation 81.1 % vs 80.0 % of inflow,
bedding 20.0 % vs 19.5 %, outflow 0.27 % vs 0.26 %, max adsorbed/Qmax 0.2559 vs 0.2559, media pool shares within
0.3 points (`tables_AN_v1/crosscheck_vs_soilp_table4_reimpl_AN_v1.csv`).

## 4. Findings about the original code (reproduced by the port on purpose, flagged `QUIRK`)
1. **PET daylength factor bug**: `b_PET.*(...)` broadcasts a row against a column; `PET_Month(m)` then always uses
   January's factor 0.82 (`monthly_pet(..., matlab_quirk=False)` gives the intended value). Found by the Octave run.
2. **P closure is not exact** in the original algorithm: solid pools gain the explicit `k_ad·dt·Cⁿ` (and explicit
   filtration) while the porewater loses the Crank–Nicolson average; inflow enters through a Dirichlet boundary.
   The authors' own `Pmassbalancechecking` residual: −4.8 % of inflow over the 90-day synthetic run; −1.7 % in the
   saturated 1-yr cross-check configuration.
3. **Plant uptake booking**: porewater uptake is applied at all depths, but only cells above `Pu_Depth` (14 cm) are
   credited to the plant pool (39,694 mg = 1.7 % of P inflow in the 90-day run, booked nowhere).
4. Water: `Waterbalancechecking` uses catchment runoff `Q_in` (before depression storage) and omits ponding storage,
   so its residual (20 % here) is catchment retention; the extended ledger closes to −1.2 × 10⁻⁵ (clamps; s is never
   capped at 1).

## 5. Augmentation for eggshell + woodchip reactors (`flowthrough_AN_v1.py`)
Saturated flow-through units (n mixed cells; Danckwerts inlet; dispersion E = (α_L/L) n Q_t) with
* design-flow **bypass** and flow-weighted remix (v4 eqs 4–5), passed exactly to the next unit of a train;
* **N module**: first-order or Michaelis–Menten denitrification, θ^(T−20) (eqs 9–11), Addy 2016 priors
  (R20 = 8.3 g N/m³ bed/d, θ = 1.080; pore basis by ÷ porosity);
* **SRP sorption with capacity depletion**: BioWEN-P kinetic-Langmuir form `R_ad = k_ad f_P (C − C_eq(q))` or the v4
  form `dq/dt = f_P (k_ads C (q_max − q) − k_des q)`; optional BioWEN-P zero-order sorbed→mineral transfer (frees sites);
* **particulate P** colloid filtration `λ v C_PP` (BioWEN-P T4.12);
* **eggshell priors from Lee et al. 2022** (q_max 11.5 g/kg column … 108.2 g/kg batch, CES-800), **not** Zhou's
  Qmax (683.5 mg/kg sandy media). K_L and k_ad are not in the Lee abstract → wide priors, must come from Maria's batch tests;
* **hard mass-balance gate** on every unit and on the train (rtol 1e-6, raises `MassBalanceError`); residuals and
  full input tables written with each run.
Every reaction is a booked transfer, limited so it cannot overshoot the closed-cell equilibrium or empty a pool.

Demo (3 SYNTHETIC years, woodchip bed 50 m³ → eggshell-woodchip mix 20 m³ with 4 wt % eggshell, prior medians):
N retained 47.9 kg/yr; P retained 0.64 / 0.73 / 0.88 kg/yr for q_max = 11.5 / 35 / 108.2 g/kg; capacity used
after 3 yr only 4.3 / 4.2 / 2.8 % because at ~0.1 g/m³ inflow the Langmuir equilibrium q = q_max K_L C/(1+K_L C)
is ~4 % of q_max for K_L = 0.5 m³/g — **K_L at field concentrations, not q_max, controls effective capacity**
(the same mechanism behind Zhou's "no saturation after 7 years"). Max relative residual of all units 4.9 × 10⁻¹⁴.
Figures: `figures_AN_v1/demo_timeseries_SYNTHETIC_AN_v2.png`, `demo_capacity_depletion_SYNTHETIC_AN_v1.png`,
`demo_mass_balance_SYNTHETIC_AN_v1.png` (`demo_timeseries_..._AN_v1.png` is superseded: legend colours).

## 6. Gaps and assumptions
* Elm Drive input data and the calibration `BioC_sbest.out` files are not public → no rerun of the published
  2012–2019 simulation; Table 6 values used. `Modelled_Data.mat` (listed in FRDR checksums) is absent from the snapshot.
* The published run used a MATLAB neural-network rating curve (`Qd_Sr_Net.mat`), not readable outside MATLAB; the
  port (and the Octave check) use the code's alternative nonlinear rating curve (FlowEstimateMethod 1).
* Not ported: Green–Ampt (unused in paper), Priestley–Taylor PET, SWMM-input branch, error/plot scripts, DDS calibration
  (to be replaced by PEST++/DREAM from bioreactor_model_AN_v1).
* Augmentation: no carbon/DO/redox, N₂O, ammonium, woodchip P leaching, pH/Ca chemistry of eggshell, clogging, or
  woodchip ageing; constant pore volume; daily forcing. Bulk density of the mixed bed (350 kg/m³) is an assumption.
