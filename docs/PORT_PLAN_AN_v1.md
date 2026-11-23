# AQUASIM → Python full port: plan (AN_v1)

Prepared 2026-11-22 for Alex Neumann. **Planning only; no port code has been written.**
Target repo: https://github.com/gudaleon/aquasimu.py (private; contains only README.md, tree sha 6defcb3…). Nothing has been pushed there by this work.

---
## 0. What this plan is based on (provenance)

| Item | Location / value |
|---|---|
| **Primary source (version 2.1g, Eawag/Reichert)** | Eawag SIAM page https://www.eawag.ch/en/department/siam/software/ → `aquasim-source_code.zip` (1,895,378 B, SHA-256 3a1a0222ff0789a252cab5183d95eaa9a8f30e71a0c53142d4bc7ff27669a621, 306 files). Copy: `/workspace/aquasim_source_original/folder1/` (extracted) and `folder2/01_Source_Audit_Eawag/aquasim-source_code.zip`. |
| **GitHub copy** | https://github.com/jsta/AQUASIM (Joseph Stachelek, archived, last commit `83dbdfe625a428c18d1a6ec8672d8146ea2ff0e7`, 2018-05-09 "typo"; 6 commits). A *mirror* of the Eawag zip, not an Eawag repo. Cloned read-only to `/workspace/aquasim_source_original/github_AQUASIM`. No GitHub licence is declared. No other GitHub AQUASIM repo (Eawag, Eawag-SIAM, Reichert, or Alex's `gudaleon`) hosts the C++ source; `gudaleon/software-warehouse` (04_AQUASIM) holds Alex's own audit notes only. |
| Drive folder 1 (13FJAr_…) | aqcore / batch / char / win32 = the Eawag zip tree. Drive listing: 300 files, all names + byte sizes equal the zip's (only `win32/res/Thumbs.db` absent from Drive). 23 files downloaded straight from Drive are md5-identical. `/workspace/aquasim_source_original/folder1` is the zip extraction (see README_PROVENANCE.txt there). |
| Drive folder 2 (1OLzZ…) | Alex's earlier audit folder: the original Eawag zips (SHA-256 verified against its own SHA256SUMS.txt), manuals, build notes, tutorial check, licence audit, a *draft, unsent* Eawag clarification email. Copied in full to `folder2/`. `03_Process_Library_Sediment_Diagenesis` is empty. |
| GitHub vs Eawag/Drive diff | `aqcore/`, `batch/`, `char/`, `win32/`: **byte-identical** (diff -rq empty). Examples/tutorials in GitHub differ **only by CRLF line endings**. GitHub lacks `tutor10/` and 10 non-.aqu tutor20 files (data/.job/.par/.def/.plt). So GitHub adds nothing beyond the Eawag zip; Eawag zip is the authority. |
| Docs | `reichert_1998_aquasim_manual.pdf` (+ text `/workspace/aquasim_examples_AN/man.txt`), tutorial, Reichert 1994 technical report (all in folder2). |
| Reference executable | Linux batch build `aquasimb` (g++ 13, `-std=gnu++98 -fpermissive`, no source change): `/workspace/aquasim_examples_AN/build_aqb/`; tutor20 outputs already generated in `/workspace/aquasim_examples_AN/run/out/`; tutorial check PASS (`folder2/04_Benchmarks/TUTORIAL_CHECK.md`). |
| Earlier partial port | `/workspace/aquasim_py_v2/` (SEDCOMP + formula evaluator + SOD_All case, 369 lines, within 1.8e-3 of reference, 198 s vs 515 s) with refs in `/workspace/aquasim_ref_v1/`. Treated as a prototype to fold into the sediment module. |
| Bioreactor model | `/workspace/bioreactor_model_AN_v1` (being built by another worker; empty when checked; not touched). It folds into Phase 2–4 (mixed, advective-dispersive, biofilm) and Phase 9 (PEST++/DREAM) below. |

Version-currency caveat (Alex's statement, recorded in folder2 SOURCE_AUDIT.md): the posted 2.1g is not Reichert's latest state; the sediment reactive-transport model Maria Dittrich used (Lake Zug) is not in it. Our benchmarks cover the posted 2.1g only.

---
## 1. Licence implications (read before choosing the approach)

* The zip has **no licence/COPYING file**. Eawag's page says AQUASIM is "freely available" since 2013, source of all versions posted, no technical support. Redistribution, modification-and-publication, or commercial/SaaS terms are **not stated**; the GitHub mirror declares none either (`license: null`). Silence is not a licence (folder2 LICENCE_AUDIT.md). Only copyright-like string in the code: `%%Copyright: P. Reichert & J. Ruchti` (drawpts.cpp output line).
* Bundled numerics: `dassl.for` (DDASSL, Petzold), `linpack.f`, `eispack.f`, `gamma.f` come from netlib (public-domain lineage). Use the **netlib originals** when we need them, not the AQUASIM copies, so the numeric kernels are licence-clean.
* Consequence for the repo:
  1. **Do not commit the Eawag C++ source, binaries, or the Eawag example/tutorial files into `aquasimu.py`.** Keep them out of git (`reference/` is gitignored; a `scripts/fetch_reference.py` downloads the official zips from Eawag and checks SHA-256).
  2. A line-by-line translation is a derivative work. Two tracks: **Track A (internal, private repo):** faithful port studying the source, for numerical identity. **Track B (anything published / open-sourced):** re-derive each module from the *published* equations (Reichert 1994/1998 manuals, papers) with a written spec per module (`docs/spec/`), written by someone who implements from the spec. The plan builds the spec docs as phase deliverables so Track B costs little extra.
  3. Public release, a public repo, a PyPI package, or a paper's "code available" statement stay **blocked until Eawag answers in writing** (a draft email exists in `folder2/05_Legal_Provenance/EAWAG_CLARIFICATION_DRAFT.md`; its recipient is unverified; it has not been sent and I did not send it). Alex decides whether to send.
  4. Repo `LICENSE`: leave unset (all rights reserved) while private. Add `NOTICE.md` crediting Reichert/Eawag and stating provenance.

---
## 2. What the original is (inventory)

**Language:** C++ (pre-C++11 / 2004 "ANSI C++"; builds with `-std=gnu++98 -fpermissive`), plus Fortran 77 numeric kernels. Four front-ends over one engine.

| Dir | Files | Lines | Role |
|---|---|---|---|
| `aqcore/` | 46 .cpp, 45 .h, 4 .f/.for, aquasim.ini | ≈94.5k C++ + ≈6.9k Fortran | the engine (all science) — **port target** |
| `batch/` | `main_b.cpp` + 11 platform makefiles | 783 | batch CLI (job files) — port as CLI |
| `char/` | 19 .cpp/18 .h | 13.7k | text-menu UI — **not ported** (replaced by Python API/CLI) |
| `win32/` | 48 .cpp/50 .h + MFC resources | 48.9k | Windows GUI — **not ported** (optional later: plotting UI) |

No unit tests exist in the source. Only test-like assets are the example/tutorial `.aqu` systems and, in Eawag's zips, `.job` / `.par` / `.def` files (the printed numbers in the 1998 Tutorial are the only published expected results).

**Dependencies:** none external (own vectors, linked lists, symbol tables, formula parser, random numbers, distributions). `.ini` file = message catalogue.

### 2.1 Class hierarchy (from the headers)
* `FILEIO` (save/load of nested-brace text blocks) → everything saveable.
* `SYMBOL` (named object) → `AQVAR`, `AQPROC`, `AQCOMP`, `AQLINK`, `CALC`, `FIT`, `PLOT`, `BIF`, `PINF`, `PARZONE`; `SYMBOLLIST` containers.
* **Variables `AQVAR`:** `STATEVAR` (dissolved/particulate, EQN type), `PROGVAR` (t, z, r, Δt, …), `CONSTVAR`, `REALLISTVAR` (interpolated data lists, 1-D curves), `VARLISTVAR`, `PROBEVAR`, `FORMVAR` (formula strings → `Fnode` parse tree: number, const, var, function, + − × ÷ idiv mod ^, comparisons, and/or/not, `if-then-else`).
* **Processes `AQPROC`:** `DYNPROC` (rate × stoichiometry), `EQUPROC` (equilibrium: algebraic equations).
* **Compartments `AQCOMP`** (virtual interface: `NumEq`, `EqType`, `InitCond`, `Delta`, `GridValue`, `SpaceValue`, `AdvInQ/OutQ`, `DiffC`): `MIXCOMP`, `FILMCOMP` (biofilm + bulk), `ADVCOMP` (advective-diffusive), `SOILCOMP` (saturated soil; `SERZONE`, `PARZONE`), `RIVCOMP` (+`RIVSEDLAYER`), `LAKECOMP` (+`LAKESEDLAYER`, `PINF`), `SEDCOMP` (sediment; ≈3.2k lines).
* **Links `AQLINK`:** `ADVLINK` (with bifurcation `BIF`), `DIFFLINK`.
* **Systems:** `VARSYS`, `PROCSYS`, `COMPSYS`, `LINKSYS`, `CALCSYS`/`CALC`, `FITSYS`/`FIT`/`FITDATA`, `NUMPAR`/`NUMINFO`, `PLOTSYS`, `STATESYS`/`STATE`, composed in `AQUASYS`.
* **Numerics:** `dassl.for/.cpp` (**DASSL**, BDF/Gear variable-step variable-order DAE solver; dense `dgefa` or banded `dgbfa` LINPACK), `zero.cpp` (`ZeroNewton`, initial-condition solve), `fit.cpp` (`MinSimplex`, `FitSimplex`, `FitSecant` = Ralston-Jennrich Gauss-Newton), `stat.cpp` (distributions, random numbers), `gamma`, `eispack`, `intpol`, `combin`.
* **Analyses in `AQUASYS`:** `Simulate`, `SensAnal` + `SensRank` (finite-difference, step = 1 % of parameter σ, per manual eq. 4.12), `Fit` / `SetFit` / `CalcChi2` / `IdAnal`, `CalcRes` / `CalcSingleRes` (uncertainty, "UNCSIM"), **`MaxPost` and `MCMC`** (added 2000: maximum posterior and Markov chain — so AQUASIM already has a Bayesian path that the DREAM layer should be checked against).
* **Batch CLI flags:** `-s` simulate, `-a` sensitivity, `-e` estimate, `-f` set of estimations, `-c` χ², `-r`/`-u` results/uncertainty, `-b` max posterior, `-m` Markov chain, `-p` plot, `-l` list, or a job file.
* **File format:** `.aqu` = `{AQUASYS}{{OPTIONS}…{VARSYS}…{PROCSYS}…{COMPSYS}…{LINKSYS}…{CALCSYS}…{FITSYS}…{NUMPAR}…{PLOTSYS}…}` nested-brace text with per-block version numbers (old versions must load); results embedded in saved systems; `.par/.def/.fit/.sen/.ran/.plt/.sim/.job` side files.

### 2.2 Examples / tests shipped by Eawag (counts from the actual `.aqu` files)
| Compartment / feature | Systems |
|---|---|
| Mixed reactor | examples: dyphora, sludgeb; tutor10: ex1–ex4 (12 files); tutor20: flow_a–c, numtim_a–c, parest_a–c, proc_a–c, modsel_a, sens_a–b |
| Biofilm reactor (**biofilter-like**) | examples: **hetaut** (heterotroph/autotroph competition), **xylene** (biofilm + mixed + diffusive link); tutor10: ex5a–b; tutor20: **film_a, film_b** |
| Advective-diffusive reactor | tutor20: advdif_a–b, numspa_a–d |
| Saturated soil column | tutor20: soil_a–d |
| River section | examples: **streeter** (Streeter-Phelps, 2 river sections); tutor10: ex6a–b; tutor20: river_a–c |
| Lake | tutor20: lake1_a–c, lake2_a–c, lake3_a–c |
| Links | advective: iawprc1, flow_*, ex1; diffusive: ex2c, proc_c, xylene |
| Sediment (`SEDCOMP`) | none in the zip; use Alex's `SOD_All`, `zug`, `HighSO4` systems already in `/workspace/aq` (+ reference profiles in `aquasim_ref_v1`) |
| Estimation / sensitivity / uncertainty | parest_a–c (+parest1/2.dat), modsel (model selection), sens_a–b, numspa/numtim (discretisation accuracy) |

---
## 3. Target architecture (repo `aquasimu.py`, import name `aquasimu`)

`src/` layout, Python ≥ 3.11, NumPy/SciPy core; strict typing (mypy), ruff, pytest.

```
aquasimu.py/
├─ pyproject.toml            # hatchling/setuptools; extras: [dae] [calib] [dev] [plot]
├─ README.md  NOTICE.md  CITATION.cff   (LICENSE deliberately absent until Eawag clarifies)
├─ .github/workflows/ci.yml  # see §6
├─ docs/  spec/ (one equation spec per module, Track B anchor)  adr/ (decisions)  parity/ (feature matrix)
├─ scripts/ fetch_reference.py  gen_reference_outputs.py (runs Eawag aquasimb → tests/ref)  
├─ reference/                # gitignored: Eawag zips, aquasimb build, original examples
├─ src/aquasimu/
│  ├─ core/        symbol.py fileio.py numpar.py types.py errors.py units.py
│  ├─ variables/   base.py (AQVAR) state.py prog.py const.py reallist.py varlist.py probe.py formula/ (lexer, parser, ast, eval; vectorised over grid)
│  ├─ processes/   base.py dyn.py equ.py stoich.py
│  ├─ compartments/ base.py mixed.py biofilm.py advdiff.py soil.py river.py lake.py sediment.py grid.py
│  ├─ links/       base.py advective.py diffusive.py bifurcation.py network.py
│  ├─ system/      aquasys.py (AQUASYS) varsys.py procsys.py compsys.py linksys.py calcsys.py statesys.py
│  ├─ solvers/     dae.py (interface) dassl_f2py.py ida.py scipy_bdf.py newton.py jacobian.py (banded FD + sparse)
│  ├─ analysis/    simulate.py sensitivity.py (+rank) estimation.py (simplex, secant) chi2.py idanal.py uncertainty.py maxpost.py mcmc.py
│  ├─ calib/       pestpp.py (pyEMU/PEST++ glm, ies, sen) dream.py (DREAM(ZS)) priors.py likelihood.py
│  ├─ checks/      massbalance.py  (per-run closure; see §5)
│  ├─ io/          aqu_reader.py aqu_writer.py (all historical block versions) job.py (job-file + CLI flags) plt_sim_par_def.py results.py (netCDF/parquet)
│  ├─ plot/        matplotlib replacements for PLOTSYS (optional extra)
│  └─ cli.py       `aquasimu -s|-a|-e|-f|-c|-r|-u|-b|-m|-p|-l|jobfile` (same flags as original batch)
├─ tests/
│  ├─ unit/        formula parser, interpolation, distributions, grid, stoichiometry
│  ├─ benchmark/   one file per compartment, comparing vs ref outputs (see §5)
│  ├─ massbalance/ closure tests run on every simulation
│  ├─ io/          round-trip every example .aqu
│  └─ ref/         small reference result tables (derived numbers only; no Eawag files)
└─ examples/       our own systems (bioreactor, etc.)
```

**Bioreactor model fold-in:** `/workspace/bioreactor_model_AN_v1` becomes `examples/bioreactor/` plus thin wrappers that build `MixedReactor`/`AdvDiffReactor`/`BiofilmReactor` objects; its DREAM/PEST++ scripts move into `calib/` once Phase 9 exists. Until then it stays independent and its equations become the first cross-check of Phase 2–4 (mass balance first).

### 3.1 Solver choice (decision with fallbacks)
AQUASIM is a method-of-lines index-1 DAE system (algebraic rows: flows, porosity, equilibrium processes, boundary conditions; banded Jacobian via finite differences). Plan a `DAESolver` interface with three back-ends, chosen by benchmark:
1. **`dassl_f2py` (reference):** wrap netlib DDASSL via f2py/gfortran (gfortran, f2py present on the box). Same algorithm as original → closest step-by-step parity, banded linear algebra, same error-test semantics. Used to make parity claims.
2. **`ida`:** SUNDIALS IDA through `scikits.odes` or `assimulo` (needs a compile step; test on CI) — fast, sparse/banded, supports sensitivities natively.
3. **`scipy_bdf` / Radau:** reduce to ODEs by explicitly eliminating the algebraic rows (the approach that already worked for SEDCOMP in `aquasim_py_v2`: 198 s vs 515 s original) — zero extra dependencies; fine for small/medium systems. SciPy has no native DAE solver, so this only works where the algebraic part is explicitly solvable; otherwise use 1 or 2.
Residuals evaluated vectorised over the grid (NumPy now; Numba/JAX optional later). Consistent-initial-condition solve = port of `ZeroNewton`/`InitCond` semantics.

### 3.2 Sensitivity / estimation / Bayesian
* `analysis/sensitivity.py` reproduces the original method (finite-difference with 1 % σ step, four sensitivity functions 4.9a–d, ranking, uncertainty contributions) and adds an optional forward-sensitivity mode (IDA) for speed.
* `analysis/estimation.py` ports Nelder-Mead simplex and Ralston-Jennrich secant with the original transformations/bounds, χ², standard errors, correlations (checked against Tutorial §3.2 numbers).
* `calib/pestpp.py`: write PEST++ control/instruction/template files from an `AquaSystem` + `FitSystem`, run `pestpp-glm`/`ies`/`sen` (binaries already at `/workspace/pestpp_bin_AN/bin`), read back with pyEMU. `calib/dream.py`: DREAM(ZS) sampler against the same log-posterior as AQUASIM's `-m` MCMC, with the same priors, so DREAM vs original MCMC is an independent consistency check (not validation).

---
## 4. Mapping: original file/class → Python module

| Original (lines) | Python module | Notes |
|---|---|---|
| `all`, `types`, `vectors`, `sdlist`, `symbol`, `combin`, `mdi`, `sort`, `readwords`, `fortran`, `init` (≈4k) | `core/*` | mostly replaced by Python builtins/NumPy; keep `SYMBOL` naming/scoping rules |
| `fileio`, `aquasys.cpp` Load/Save (≈0.5k + part of 8.6k) | `io/aqu_reader.py`, `aqu_writer.py` | must read all `aquavers` block versions (code has case 2,3,… branches) |
| `aqvar`, `othervar` (3.5k), `formvar` (1.9k), `intpol` | `variables/*`, `variables/formula/` | formula grammar from `Fnode` classes; real-list variables with LINEAR/SPLINE interpolation |
| `aqproc`, `procsys` (1.2k) | `processes/*` | stoichiometric matrix assembly |
| `aqcomp`, `compsys` (1.5k) | `compartments/base.py` | the virtual interface above |
| `mixcomp` (1.7k) | `compartments/mixed.py` | Phase 2 |
| `advcomp` (3.1k) | `compartments/advdiff.py` | FD/FV scheme, dispersion, numerical-diffusion options |
| `filmcomp` (6.5k) | `compartments/biofilm.py` | biofilm grid + bulk phase, growth/detachment, moving boundary |
| `soilcomp` (5.3k) | `compartments/soil.py` | porewater/sorption/exchange zones |
| `rivcomp` (6.4k) | `compartments/river.py` | open-channel flow (kinematic/diffusive wave), sediment layers |
| `lakecomp` (9.1k) | `compartments/lake.py` | stratified water column, mixing model (k-ε style options), sediment layers, inflow/outflow |
| `sedcomp` (3.2k) | `compartments/sediment.py` | fold in `aquasim_py_v2` |
| `aqlink`, `linksys` (2.6k) | `links/*` | advective, bifurcations, diffusive |
| `calcsys`, `statesys`, `numpar` (1.2k) | `system/calcsys.py`, `statesys.py`, `core/numpar.py` | |
| `dassl.for/.cpp` (3.4k+2.4k), `linpack`, `zero`, `eispack` | `solvers/*` | use netlib originals (§1) |
| `aquasys.cpp` analyses (rest of 8.6k) | `analysis/*`, `system/aquasys.py` | |
| `fit`, `fitsys` (2.3k) | `analysis/estimation.py` | |
| `stat`, `gamma` (5.9k) | `core/distributions.py` (+ scipy.stats cross-check) | beta/gamma/normal/lognormal/uniform etc. + RNG; seeds will differ from the original (compare distributions, not draws) |
| `plotsys`, `plotpts`, `drawpts` (4.5k) | `plot/` (matplotlib) | optional; `.plt` command files still parsed |
| `batch/main_b.cpp` (0.8k) | `cli.py`, `io/job.py` | |
| `char/*`, `win32/*` (62k) | — | not ported |

---
## 5. Testing strategy (Alex's first principle: mass balance first; comparisons are consistency checks, not validation)

1. **Mass-balance closure on every run (all modules, hard gate).** `checks/massbalance.py` hooks into `Simulate()` and every benchmark: for each state variable and compartment, ∑(storage change) = ∑(advective + diffusive + boundary in/out) + ∑(process rate × stoichiometry × volume) within a stated tolerance; reports the closure residual (absolute and relative) in the result object and the log; **a run whose residual exceeds tolerance fails the test and is flagged in output**. Conservative-tracer cases check total mass exactly; processes check element conservation via stoichiometry (N, P, C, O columns). Water balance (Q closure) checked first, then solute, then particulates.
2. **Reference runs from the original:** `scripts/gen_reference_outputs.py` runs the Eawag-built `aquasimb` on every example/tutorial job and stores compact reference tables (state profiles/time series, χ², parameters, sensitivity rankings) in `tests/ref/` as derived numbers. Tolerances per case: ODE-level cases rtol ≈ 1e-6 vs same-solver, ≈ 1e-3 vs different step history (the SEDCOMP port reached ≤1.8e-3 of profile max with 4-sig-digit references).
3. **Per-compartment benchmarks** (minimum set; all must close mass balance *and* match):

| Phase | Benchmark systems (Eawag originals unless noted) | Compared quantities |
|---|---|---|
| Mixed + processes + links | tutor20 `proc_a–c`, `flow_a–c`, `numtim_a–c`; `dyphora`, `sludgeb`, `iawprc1` | concentrations(t); Monod/batch analytic solutions; tracer mass |
| Advective-diffusive | `advdif_a–b`, `numspa_a–d` | profiles vs analytic Gaussian/advection solutions + reference; numerical diffusion vs mesh |
| **Biofilm / biofilter** | **`film_a–b`, `hetaut`, `xylene`, tutor10 `ex5a–b`** | substrate/biomass profiles in z, bulk concentrations, biofilm thickness; flux balance bulk↔film |
| Saturated soil | `soil_a–d` | breakthrough curves, sorption/pore-exchange mass |
| River | `river_a–c`, `streeter`, tutor10 `ex6a–b` | Q, depth, C(x,t); Streeter-Phelps DO sag vs analytic |
| Lake | `lake1_a–c` (mixing), `lake2_a–c`, `lake3_a–c` | C(z,t), Kz profile, tracer mass vs analytic Gaussian diffusion |
| Sediment | Alex's `SOD_All`, `zug`, `HighSO4` + `aquasim_ref_v1` profiles | porewater profiles (existing prototype passes 1.8e-3) |
| Sensitivity | `sens_a–b`, parest-based | sensitivity functions + ranking vs original files |
| Estimation | `parest_a–c` (Tutorial §3.2: K=1.38, rmax=1.15, χ²=14.6) , `modsel_a` | parameters, χ², std errors, correlation |
| Uncertainty / Bayes | `-u`, `-b`, `-m` on parest cases | posterior mode/variance; DREAM vs AQUASIM MCMC marginals |
| File I/O | every `.aqu` (61 shipped) | read → write → read round-trip identical; original binary can reload our output (`aquasimb -s` on written file) |

4. **Property/analytic tests:** conservative tracers, pure decay, steady-state advection-reaction, symmetric diffusion, grid-refinement convergence orders.
5. **Solver parity:** same case on all back-ends; differences reported, not hidden.
6. Every benchmark records which back-end/tolerance was used; a failing parity test blocks the phase gate.

---
## 6. Repo mechanics (tests, CI)

* **CI (`.github/workflows/ci.yml`, GitHub Actions, private repo minutes):** matrix Python 3.11/3.12/3.13 on ubuntu-latest; jobs: `lint` (ruff, mypy), `unit` (pytest -m "not slow"), `benchmark-fast` (cases < 30 s, uses committed `tests/ref/` numbers; needs no Eawag files), `benchmark-full` (nightly/manual; downloads Eawag zips via `scripts/fetch_reference.py`, builds `aquasimb`, regenerates references, runs everything incl. SEDCOMP SOD ≈ minutes), `f2py` build job (gfortran), coverage report. Branch protection: mass-balance tests are required checks.
* Pre-commit: ruff, mypy, a check that blocks Eawag files/binaries (`aqcore/`, `*.aqu` from Eawag, `aquasimb`) from being committed.
* Conventions: versioned filenames (`_AN_vX`) for non-code deliverables; ADRs for each solver/scheme decision; each phase tagged (`phase-N-complete`) after its gate.
* **Commit workflow:** a separate coding agent commits to feature branches + PRs into `main`; this plan itself does not push. Suggested first PR: scaffold + CI + `docs/spec` stubs + this plan copied as `docs/PORT_PLAN.md`.

---
## 7. Phases, effort, gates

Effort = focused coding-agent work (implementation + tests) plus a human-review allowance; estimates are ranges, not commitments. Calendar assumes Alex reviews at each gate (≈1–2 days each). Lines-of-original are the main sizing driver (≈1 working day per ~1–1.5k original lines of tricky numerical code including tests).

| # | Phase | Scope | Effort | Gate (all include mass balance) |
|---|---|---|---|---|
| 0 | Foundations & repo | scaffold, CI, fetch/reference-gen scripts, NOTICE, spec stubs, `aquasimb` reference runs for all examples; licence decision recorded | 2–3 d | CI green; reference tables for all 61 shipped systems |
| 1 | Core model + I/O | symbols, variables, formula parser, interpolation, processes/stoichiometry, `.aqu` reader/writer (all versions), job files, results container, `massbalance.py` skeleton | 8–12 d | all `.aqu` round-trip; formula tests; reload by original binary |
| 2 | Solver layer | DAE interface; `dassl_f2py`, scipy-BDF reduction, IDA option; `ZeroNewton`/consistent initial conditions; banded Jacobian | 5–8 d | solver parity on mixed-reactor ODE cases; step-count comparable to DASSL |
| 3 | Mixed reactor + links + network | MIXCOMP, ADVLINK (junction/bifurcation), DIFFLINK, calc series | 4–6 d | tutor20 proc/flow/numtim, dyphora, sludgeb, iawprc1 pass + closure |
| 4 | Advective-diffusive + **bioreactor fold-in** | ADVCOMP; fold in the bioreactor model as example; verify its equations in the framework (mass balance first) | 6–9 d | advdif/numspa pass; bioreactor consistency check vs framework |
| 5 | **Biofilm reactor** | FILMCOMP (biofilm grid, moving boundary, bulk coupling) | 10–15 d | film_a–b, hetaut, xylene, ex5 pass; flux balance bulk↔film closes |
| 6 | Saturated soil column | SOILCOMP zones, sorption, pore exchange | 8–12 d | soil_a–d pass |
| 7 | River section | RIVCOMP hydraulics + transport + sediment layers | 10–15 d | river_a–c, streeter, ex6 pass (Streeter-Phelps analytic check) |
| 8 | Lake + sediment | LAKECOMP (mixing, ≈9k original lines), SEDCOMP (fold in v2 prototype), sediment layers | 15–22 d | lake1–3, SOD_All/zug match prior reference |
| 9 | Sensitivity, estimation, uncertainty | SensAnal/SensRank, simplex+secant, χ²/IdAnal, CalcRes (UNCSIM), MaxPost, MCMC | 10–15 d | Tutorial §3 numbers reproduced; ranking files match |
| 10 | PEST++ and DREAM layers | pyEMU/PEST++ glm/ies/sen adapter, DREAM(ZS), shared priors/likelihood; compare against ported MaxPost/MCMC; apply to bioreactor | 6–9 d | synthetic-recovery tests; DREAM vs MCMC marginals consistent; PEST++ vs secant consistent |
| 11 | CLI, plotting, docs, hardening | CLI flags, plot replacements, parity matrix, performance (Numba/sparse), packaging (private), user docs | 6–10 d | `aquasimu` CLI runs every original job file; parity matrix complete |
| | **Total** | | **≈ 90–135 agent-days ≈ 4–6 calendar months** with reviews; **core usable subset (phases 0–6 + 9–10: mixed, adv-diff, biofilm, soil, estimation, PEST++/DREAM) ≈ 55–80 d ≈ 2.5–3.5 months.** Lake and river are the largest/riskiest remaining blocks. |

Suggested order if Alex wants the bioreactor application early: 0 → 1 → 2 → 3 → 4 → 5 → 9 → 10 (this delivers biofilter/bioreactor with Bayesian calibration), then 6, 7, 8, 11.

---
## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| **Licence unclear** (no terms in zip or GitHub mirror) | Private repo; no Eawag files in git; Track A/B split with spec docs; written Eawag clarification before any publication; use netlib DDASSL/LINPACK originals |
| Source is 2.1g, not Reichert's latest; sediment RTM absent | Benchmarks scoped to 2.1g; second question in the Eawag email; sediment benchmark from Alex's SOD/Zug systems |
| Numerical parity (DASSL step history, error control, different Jacobians) | DASSL via f2py as reference back-end; define tolerances up-front; compare profiles not step counts |
| Performance (the SOD case: 7,440 unknowns, 515 s in C++) | banded/sparse Jacobians, vectorised residuals, IDA/Numba options; benchmark budgets in CI |
| Hidden behaviour in 9k-line lake/6k-line river/film code (no tests, sparse comments) | spec-first per module; golden outputs from the original binary for many parameter variations, not just the shipped examples; generate extra stress cases |
| Old `.aqu` block-version zoo | systematic reader for each version branch; round-trip tests against all shipped files; reload test with original |
| Random numbers (MCMC, uncertainty) cannot match bitwise | compare distributions/summary statistics with stated tolerances |
| Formula language edge cases (units, `if` semantics, program variables) | exhaustive parser tests from manual; evaluation vectorised and cross-checked against the C++ evaluator on random expressions |
| Scope creep (GUI, plotting) | GUI explicitly out of scope; plotting last |
| Parallel work (bioreactor model, older `aquasim_py_*` prototypes) | fold-in plan above; no edits to others' directories |

---
## 9. Open decisions for Alex

1. Send (or not) the Eawag licence-clarification email? (Draft exists, recipient unverified, unsent.) Until answered: private repo, Track A internal + Track B specs.
2. Priority order: bioreactor-first (0–5, 9, 10) or compartment-complete order?
3. Solver: accept the three-back-end plan (DASSL-f2py reference + IDA + SciPy)? IDA requires a SUNDIALS build on CI.
4. Is the GUI/plotting UI wanted at all (currently matplotlib-only, last phase)?
5. Whether to ask Maria Dittrich for her compiled sediment-RTM model as a benchmark-only artefact (Alex's call; no contact made).
