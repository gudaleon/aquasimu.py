"""Build the PEST++ control file with pyEMU (template / instruction files, weights, options)."""
import sys, pathlib, shutil, os
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import pyemu
from bioreactor_AN_v1 import synthetic_AN_v1 as S
from bioreactor_AN_v1 import calib_common_AN_v1 as C

P = HERE / "pest_AN_v1"
START = dict(vmax_N=3.0, K_N=8.0, theta_N=1.04, k_ads=0.002, q_max=2.0)   # deliberately wrong start values


def build():
    # template file
    with open(P / "params_AN_v1.tpl", "w") as f:
        f.write("ptf ~\n")
        for n in S.NAMES:
            f.write(f"{n.lower()} ~ {n.lower():<10s} ~\n")
    with open(P / "params_AN_v1.dat", "w") as f:           # example input (start values)
        for n in S.NAMES:
            f.write(f"{n.lower()} {START[n]:.10e}\n")
    # instruction file (reads  name value  lines)
    names = [f"on_{d:03d}" for d in C.D] + [f"op_{d:03d}" for d in C.D] + ["fkgn", "fkgp"]
    with open(P / "output_AN_v1.ins", "w") as f:
        f.write("pif ~\n")
        for n in names:
            f.write(f"l1 w !{n}!\n")
    # dummy output so pyEMU can verify the instruction file
    ref = S.run_model(START_VEC, C.FORCING)
    with open(P / "output_AN_v1.dat", "w") as f:
        for d in C.D: f.write(f"on_{d:03d} {ref['cN'][d]:.10e}\n")
        for d in C.D: f.write(f"op_{d:03d} {ref['cP'][d]:.10e}\n")
        f.write(f"fkgn {ref['kgN']:.10e}\n"); f.write(f"fkgp {ref['kgP']:.10e}\n")

    pst = pyemu.Pst.from_io_files([str(P / "params_AN_v1.tpl")], [str(P / "params_AN_v1.dat")],
                                  [str(P / "output_AN_v1.ins")], [str(P / "output_AN_v1.dat")], pst_path=".")
    par = pst.parameter_data
    for i, n0 in enumerate(S.NAMES):
        n = n0.lower()
        par.loc[n, "parval1"] = START[n0]
        par.loc[n, "parlbnd"] = S.LO[i]; par.loc[n, "parubnd"] = S.HI[i]
        par.loc[n, "partrans"] = "log" if S.LOG[i] else "none"
        par.loc[n, "pargp"] = "model"
        par.loc[n, "parchglim"] = "factor"
    obs = pst.observation_data
    for d, o, sg in zip(C.D, C.OBS["obsN"], C.SIG_N):
        obs.loc[f"on_{d:03d}", ["obsval", "weight", "obgnme"]] = [o, 1.0 / sg, "nitrate"]
    for d, o, sg in zip(C.D, C.OBS["obsP"], C.SIG_P):
        obs.loc[f"op_{d:03d}", ["obsval", "weight", "obgnme"]] = [o, 1.0 / sg, "phosphate"]
    obs.loc["fkgn", ["obsval", "weight", "obgnme"]] = [C.OBS["true_kgN"], 0.0, "forecast"]
    obs.loc["fkgp", ["obsval", "weight", "obgnme"]] = [C.OBS["true_kgP"], 0.0, "forecast"]
    pst.model_command = ["python forward_run_AN_v1.py"]   # needs numpy/numba on PATH python (see README)
    pst.control_data.noptmax = 12
    pst.pestpp_options["forecasts"] = "fkgn,fkgp"
    pst.pestpp_options["uncertainty"] = True
    pst.pestpp_options["max_run_fail"] = 3
    pst.pestpp_options["ies_num_reals"] = 200
    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options["ies_no_noise"] = False
    pst.rectify_pgroups()
    pst.write(str(P / "pest_AN_v1.pst"), version=1)
    return pst


START_VEC = np.array([START[n] for n in S.NAMES])
if __name__ == "__main__":
    pst = build(); print(pst.npar, pst.nobs, "obs;", "written")
