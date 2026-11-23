"""PEST++ control files for the SYNTHETIC bioreactor case, via pyEMU.

The PEST++ executables (``pestpp-glm``, ``pestpp-ies``) are not part of this
package and must not be committed. Put them on ``PATH``, or set ``PESTPP_BIN``
to the directory that contains them. Functions that launch a binary raise
``FileNotFoundError`` when it is missing; the test suite skips those cases.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

START = {
    "vmax_N": 3.0,
    "K_N": 8.0,
    "theta_N": 1.04,
    "k_ads": 0.002,
    "q_max": 2.0,
}  # deliberately away from the synthetic truth

PST_NAME = "bioreactor.pst"
PARAM_TPL = "params.tpl"
PARAM_DAT = "params.dat"
OUTPUT_INS = "output.ins"
OUTPUT_DAT = "output.dat"


def model_command() -> str:
    """Command PEST++ runs in the case directory.

    PEST++ strips a leading ``/`` on the model-command line, so an absolute
    interpreter path is not usable. The command is the interpreter's file name
    (``python`` or ``python3``). ``run_glm`` and ``run_ies`` put that
    interpreter's directory on ``PATH`` before launching PEST++.
    """
    return f"{Path(sys.executable).name} -m aquasimu.examples.bioreactor.forward_run"


def _ensure_interpreter_on_path() -> None:
    bindir = str(Path(sys.executable).parent)
    parts = os.environ.get("PATH", "").split(os.pathsep)
    if parts[:1] != [bindir]:
        os.environ["PATH"] = bindir + os.pathsep + os.environ.get("PATH", "")


def pest_executable(name: str) -> str | None:
    """Absolute path of a PEST++ binary, or None when it is not available."""
    directory = os.environ.get("PESTPP_BIN")
    if directory:
        candidate = Path(directory) / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return shutil.which(name)


def require_pest(*names: str) -> dict[str, str]:
    found = {name: pest_executable(name) for name in names}
    missing = [name for name, path in found.items() if path is None]
    if missing:
        raise FileNotFoundError(
            "PEST++ binaries not found: "
            + ", ".join(missing)
            + ". Install them and put the directory on PATH or in PESTPP_BIN. "
            "Do not commit the binaries."
        )
    return found  # type: ignore[return-value]


def build_case(workdir: Path, noptmax: int = 12, ies_num_reals: int = 200):
    """Write template, instruction, input and control files for the synthetic case.

    Requires pyEMU. Does not launch PEST++. ``noptmax=0`` asks PEST++ to run the
    model once at the starting parameters.
    """
    import pyemu

    from aquasimu.examples.bioreactor import synthetic as case
    from aquasimu.examples.bioreactor.likelihood import D, FORCING, OBS, SIG_N, SIG_P

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    start_vector = np.array([START[name] for name in case.NAMES])

    with open(workdir / PARAM_TPL, "w") as handle:
        handle.write("ptf ~\n")
        for name in case.NAMES:
            handle.write(f"{name.lower()} ~ {name.lower():<10s} ~\n")
    with open(workdir / PARAM_DAT, "w") as handle:
        for name in case.NAMES:
            handle.write(f"{name.lower()} {START[name]:.10e}\n")

    names = [f"on_{day:03d}" for day in D] + [f"op_{day:03d}" for day in D] + ["fkgn", "fkgp"]
    with open(workdir / OUTPUT_INS, "w") as handle:
        handle.write("pif ~\n")
        for name in names:
            handle.write(f"l1 w !{name}!\n")

    reference = case.run_model(start_vector, FORCING)
    with open(workdir / OUTPUT_DAT, "w") as handle:
        for day in D:
            handle.write(f"on_{day:03d} {reference['cN'][day]:.10e}\n")
        for day in D:
            handle.write(f"op_{day:03d} {reference['cP'][day]:.10e}\n")
        handle.write(f"fkgn {reference['kgN']:.10e}\n")
        handle.write(f"fkgp {reference['kgP']:.10e}\n")

    pst = pyemu.Pst.from_io_files(
        [str(workdir / PARAM_TPL)],
        [str(workdir / PARAM_DAT)],
        [str(workdir / OUTPUT_INS)],
        [str(workdir / OUTPUT_DAT)],
        pst_path=".",
    )
    parameters = pst.parameter_data
    for index, name in enumerate(case.NAMES):
        key = name.lower()
        parameters.loc[key, "parval1"] = START[name]
        parameters.loc[key, "parlbnd"] = case.LO[index]
        parameters.loc[key, "parubnd"] = case.HI[index]
        parameters.loc[key, "partrans"] = "log" if case.LOG[index] else "none"
        parameters.loc[key, "pargp"] = "model"
        parameters.loc[key, "parchglim"] = "factor"
    observations = pst.observation_data
    for day, value, sigma in zip(D, OBS["obsN"], SIG_N):
        observations.loc[f"on_{day:03d}", ["obsval", "weight", "obgnme"]] = [value, 1.0 / sigma, "nitrate"]
    for day, value, sigma in zip(D, OBS["obsP"], SIG_P):
        observations.loc[f"op_{day:03d}", ["obsval", "weight", "obgnme"]] = [value, 1.0 / sigma, "phosphate"]
    observations.loc["fkgn", ["obsval", "weight", "obgnme"]] = [OBS["true_kgN"], 0.0, "forecast"]
    observations.loc["fkgp", ["obsval", "weight", "obgnme"]] = [OBS["true_kgP"], 0.0, "forecast"]
    pst.model_command = [model_command()]
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["forecasts"] = "fkgn,fkgp"
    pst.pestpp_options["uncertainty"] = True
    pst.pestpp_options["max_run_fail"] = 3
    pst.pestpp_options["ies_num_reals"] = ies_num_reals
    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options["ies_no_noise"] = False
    pst.rectify_pgroups()
    pst.write(str(workdir / PST_NAME), version=1)
    return pst


def run_glm(workdir: Path) -> subprocess.CompletedProcess:
    """Run pestpp-glm on the control file already written in ``workdir``.

    The iteration limit is the ``noptmax`` stored in that control file
    (``noptmax=0`` is a single forward run).
    """
    _ensure_interpreter_on_path()
    binary = require_pest("pestpp-glm")["pestpp-glm"]
    return subprocess.run(
        [binary, PST_NAME],
        cwd=workdir,
        check=True,
        env=os.environ.copy(),
    )


def run_ies(
    master_dir: Path,
    template_dir: Path,
    *,
    noptmax: int = 4,
    num_reals: int = 200,
    num_workers: int = 8,
    port: int = 4004,
):
    """Run pestpp-ies with pyEMU workers. The template directory must contain a built case.

    This is the synthetic-case ensemble smoother. It is not a validation.
    """
    import pyemu

    _ensure_interpreter_on_path()
    binary = require_pest("pestpp-ies")["pestpp-ies"]
    template_dir = Path(template_dir)
    master_dir = Path(master_dir)
    pst = pyemu.Pst(str(template_dir / PST_NAME))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_subset_size"] = min(8, num_reals)
    pst.pestpp_options["ies_lambda_mults"] = "0.5,1.0,2.0"
    pst.pestpp_options["lambda_scale_fac"] = "0.7,1.0"
    pst.pestpp_options.pop("uncertainty", None)
    pst.write(str(template_dir / PST_NAME), version=1)
    worker_root = master_dir.parent / "workers_tmp"
    worker_root.mkdir(exist_ok=True)
    try:
        pyemu.os_utils.start_workers(
            str(template_dir),
            binary,
            PST_NAME,
            num_workers=num_workers,
            worker_root=str(worker_root),
            master_dir=str(master_dir),
            port=port,
            verbose=False,
        )
    finally:
        shutil.rmtree(worker_root, ignore_errors=True)
