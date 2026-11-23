# Exclusions

`aquasim_out/` (AQUASIM 2.1g list files) is not in this repository. The licence
status of AQUASIM and its output listings is unclear, so they are not committed.

To regenerate listings locally, `scripts/run_aquasim_examples_AN_v1.py` downloads
the public AQUASIM distribution into a local working directory outside this tree
(default `/workspace/aquasim_examples_AN`) and writes listings under `aquasim_out/`.
That script is ours; it does not vendor Eawag source, binaries, or `.aqu` examples.
