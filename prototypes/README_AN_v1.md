# Prototypes archived from the box (November 2026)

These folders are the working trees used to build the bioreactor-first phase
of aquasimu.py and the BioWEN-P port. They are archived here so the code is
not only on a workspace disk. The installable package itself is in pull
request #1 (`src/aquasimu/`).

| Path | What it is |
|---|---|
| `bioreactor_model_AN_v1/` | Original woodchip/eggshell reactor + DREAM(ZS) + PEST++ prototype (`bioreactor_AN_v1`), tests, synthetic figures and tables |
| `biowenp_py_AN_v1/` | Python port of BioWEN-P (Zhou et al. 2023) plus flow-through N/P augmentation. **LGPL-3.0-or-later**; see that folder's `LICENSE` and `NOTICE_AN_v1.md` |
| `aquasim_py_v1/`, `aquasim_py_v2/` | Early Python sketches of an AQUASIM engine (our code: `.aqu` reader, formula evaluator, SEDCOMP SOD experiment). No Eawag source is vendored |

`docs/PORT_PLAN_AN_v1.md` is the port plan. `docs/maria_modelling_section/` holds
Markdown drafts of the proposal modelling section (v1–v4). Word files and
base64 email attachments are not committed.

## Exclusions

- Eawag AQUASIM source, binaries, `.aqu` examples, and output listings
  (the listings that lived in `bioreactor_model_AN_v1/benchmarks_AN_v1/aquasim_out/`
  are omitted; see the README in that folder)
- The FRDR BioWEN-P MATLAB snapshot
- BioWEN-P equivalence binaries `results_AN_v1/*.mat` and `*.npz` (~28 MB and
  ~50 MB); regenerate with `scripts/run_equivalence_octave_AN_v1.py` and a
  local FRDR checkout (documented in that folder)
- Secrets, tokens, and PEST++ binaries
