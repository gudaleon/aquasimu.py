"""Optional consistency check against user-supplied AQUASIM listings.

Skipped unless AQUASIMU_AQUASIM_LISTINGS points at a directory the user filled.
This repository does not contain those listings, Eawag source, or .aqu files.
The comparison is not a validation.
"""
import os

import pytest

from aquasimu.benchmarks.aquasim_consistency import LISTING_ENV, run_consistency


@pytest.mark.aquasim
def test_aquasim_listing_consistency():
    listing_dir = os.environ.get(LISTING_ENV)
    if not listing_dir:
        pytest.skip(
            f"Set {LISTING_ENV} to a directory of AQUASIM listing files. "
            "Listings are not shipped. The comparison is a consistency check, not a validation."
        )
    table = run_consistency(listing_dir)
    residuals = table[["N_rel_resid", "P_rel_resid", "water_rel_resid"]].dropna()
    assert (residuals < 1e-6).all().all()
    b1 = table.loc[table.quantity.str.startswith("C_A(t) python")]
    assert float(b1.RMSE_pct_of_scale.iloc[0]) < 0.01
