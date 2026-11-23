# Equivalence-run artefacts

The 90-day Octave vs Python consistency check writes two large binaries that are **not** in this repository
(see `.gitignore`):

| File | Typical size | How it is made |
|---|---|---|
| `octave_out_<N>d_<tag>.mat` | ~28 MB for 90 d | GNU Octave running the original BioWEN-P MATLAB from a local FRDR snapshot |
| `python_faithful_<N>d_<tag>.npz` | ~50 MB for 90 d | the Python port on the same synthetic forcing |

The small inputs and the comparison tables **are** committed:

- `synthetic_forcing_10min_90d_AN_v1.csv`
- `params_matlab_names_AN_v1.json`
- `../tables_AN_v1/equivalence_octave_vs_python_90d_AN_v1.csv`
- `../tables_AN_v1/equivalence_summary_90d_AN_v1.json`

## Regenerate (does not overwrite an existing tag)

You need GNU Octave (>= 7; 9.4.0 was used) and the original model code from
[FRDR doi:10.20383/103.0640](https://doi.org/10.20383/103.0640), unpacked so that
`ModelMainloop.m` and `BioC_P_RTM.m` sit in a snapshot directory. The driver defaults to
`/workspace/parsons_bioretention_AN_v1/frdr_8c090e43_model_code` (edit `SNAP` in
`scripts/run_equivalence_octave_AN_v1.py` if your path differs). The original `.m` files
are not redistributed here.

```bash
python3 scripts/run_equivalence_octave_AN_v1.py 90 new_tag
python3 scripts/plot_equivalence_AN_v1.py   # if you also want a new figure; use a new filename
```

The script refuses to overwrite an existing `octave_out_90d_<tag>.mat`. This is a
code-equivalence (consistency) check on **synthetic** forcing, not a validation against field data.
