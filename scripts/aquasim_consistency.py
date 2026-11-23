"""Compare the reactor with user-supplied AQUASIM listings.

    AQUASIMU_AQUASIM_LISTINGS=/path/to/listings python scripts/aquasim_consistency.py

The path is not searched inside this repository. Listings, Eawag source and
.aqu files are not shipped. The table written to outputs/tables is a
consistency check, not a validation. This does not run a soil, river or lake
compartment: the sorption-column rows use the advective-diffusive reactor.
"""
import sys
from pathlib import Path

from aquasimu.benchmarks.aquasim_consistency import LISTING_ENV, listing_dir_from_env, run_consistency

listing_dir = listing_dir_from_env()
if listing_dir is None:
    sys.exit(f"Set {LISTING_ENV} to a directory of AQUASIM listing files.")
table = run_consistency(listing_dir)
out = Path("outputs/tables")
out.mkdir(parents=True, exist_ok=True)
destination = out / "aquasim_consistency.csv"
table.to_csv(destination, index=False, float_format="%.6g")
columns = [c for c in ("case", "quantity", "RMSE_pct_of_scale", "max_err_pct_of_scale") if c in table.columns]
print(table[columns].to_string())
print(f"wrote {destination}")
