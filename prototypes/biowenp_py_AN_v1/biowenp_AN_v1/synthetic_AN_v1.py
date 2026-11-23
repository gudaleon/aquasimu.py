# Copyright 2020-2023 Bowen Zhou (original MATLAB BioWEN-P); Python port 2026 Alex Neumann / EnviroManagement.
# Part of biowenp_py_AN_v1, a Python port of BioWEN-P (Zhou et al. 2023, doi:10.1016/j.scitotenv.2023.162749;
# code: FRDR doi:10.20383/103.0640). GNU LGPL v3.0 or later (see COPYING.LESSER, COPYING). NO WARRANTY.
"""SYNTHETIC forcing generators. Every array produced here is invented (seeded random numbers) and is labelled
SYNTHETIC in tables and figures. Nothing here is Elm Drive monitoring data (that input data set is not public)."""
import numpy as np
from datetime import datetime, timedelta


def synthetic_stormwater_10min(n_days=60, start="2013-04-01", seed=7, storm_rate_per_day=0.25,
                               cold_spell=True, extra_storms=((6, 20.0, 36), (20, 90.0, 36), (45, 60.0, 24))):
    """10-minute precipitation (mm/10 min) and air temperature (degC) for the BioWEN-P hydrology.
    Storms: Poisson arrivals, exponential depth (mean 12 mm), 2-8 h duration. Temperature: seasonal + diurnal;
    optional 3-day cold spell below the snow threshold to exercise snow accumulation/melt."""
    rng = np.random.default_rng(seed)
    n = n_days * 144
    t0 = datetime.fromisoformat(start)
    times = [t0 + timedelta(minutes=10 * k) for k in range(n)]
    doy = np.array([tt.timetuple().tm_yday + tt.hour / 24 + tt.minute / 1440 for tt in times])
    temp = 8 + 10 * np.sin(2 * np.pi * (doy - 110) / 365) + 4 * np.sin(2 * np.pi * (doy % 1 - 0.375))
    if cold_spell and n_days > 10:
        temp[5 * 144:8 * 144] -= 12.0
    precip = np.zeros(n)
    for d in range(n_days):
        for _ in range(rng.poisson(storm_rate_per_day)):
            start_k = d * 144 + rng.integers(0, 144)
            dur = rng.integers(12, 48)
            depth = rng.exponential(12.0)
            shape = np.sin(np.linspace(0, np.pi, dur)) + 1e-3
            k1 = min(n, start_k + dur)
            precip[start_k:k1] += depth * shape[:k1 - start_k] / shape.sum()
    for day, depth, dur in extra_storms:          # deterministic events: snow during the cold spell, overflow storm
        if day < n_days:
            k0 = day * 144 + 60
            k1 = min(n, k0 + dur)
            precip[k0:k1] += depth / dur
    month = np.array([tt.month for tt in times])
    datenum = np.array([tt.toordinal() + 366 + (tt.hour * 60 + tt.minute) / 1440 for tt in times])
    return dict(precip=precip, temp=temp, month=month, datenum=datenum, times=times)


def synthetic_reactor_daily(n_days=365, seed=11, q_design=50.0):
    """Daily flow (m3/d), water temperature (degC), nitrate-N, SRP and particulate P (g/m3) at a wetland outflow.
    Includes storm peaks above the design flow so the bypass is exercised. SYNTHETIC."""
    rng = np.random.default_rng(seed)
    d = np.arange(n_days)
    T = 11 + 9 * np.sin(2 * np.pi * (d - 110) / 365) + rng.normal(0, 0.8, n_days)
    T = np.clip(T, 0.5, None)
    base = 25 + 10 * np.sin(2 * np.pi * (d - 60) / 365)
    storms = np.zeros(n_days)
    for k in rng.choice(n_days, 16, replace=False):
        peak = rng.uniform(40, 90)
        for j in range(4):
            if k + j < n_days:
                storms[k + j] += peak * np.exp(-j)
    Q = np.clip(base + storms + rng.normal(0, 2, n_days), 2, None)
    CN = np.clip(6 + 3 * np.sin(2 * np.pi * (d - 20) / 365) + 0.04 * storms + rng.normal(0, 0.4, n_days), 0.2, None)
    CP = np.clip(0.08 + 0.002 * storms + rng.normal(0, 0.01, n_days), 0.005, None)
    CPP = np.clip(0.05 + 0.003 * storms + rng.normal(0, 0.01, n_days), 0.0, None)
    return dict(Q=Q, T=T, C_N=CN, C_SRP=CP, C_PP=CPP)
