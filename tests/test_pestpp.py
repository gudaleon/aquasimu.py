"""PEST++ setup for the synthetic case.

Building the control file needs pyEMU and does not need the PEST++ binaries.
Launching pestpp-glm or pestpp-ies is skipped unless both are on PATH (or in
PESTPP_BIN). The binaries are not part of the repository.
"""
import shutil
import subprocess
import sys

import pytest

from aquasimu.calib.pestpp import PST_NAME, build_case, pest_executable, run_glm


def test_build_synthetic_control_file_and_forward_run(tmp_path):
    pytest.importorskip("pyemu")
    work = tmp_path / "case"
    pst = build_case(work, noptmax=0)
    assert pst.npar == 5
    assert pst.nobs == 244
    assert "forward_run" in pst.model_command[0]
    assert (work / PST_NAME).is_file()
    subprocess.run(
        [sys.executable, "-m", "aquasimu.examples.bioreactor.forward_run"],
        cwd=work,
        check=True,
    )
    lines = (work / "output.dat").read_text().splitlines()
    assert len(lines) == 244
    assert lines[-2].startswith("fkgn ")
    assert lines[-1].startswith("fkgp ")


@pytest.mark.pest
def test_pestpp_glm_single_forward_run(tmp_path):
    if pest_executable("pestpp-glm") is None or pest_executable("pestpp-ies") is None:
        pytest.skip(
            "pestpp-glm and pestpp-ies are not on PATH. "
            "Set PESTPP_BIN to their directory. Binaries are not committed."
        )
    pytest.importorskip("pyemu")
    work = tmp_path / "glm"
    build_case(work, noptmax=0)
    run_glm(work)
    # A successful single forward run leaves a residuals file.
    written = list(work.glob("*.rei")) + list(work.glob("*.rec"))
    assert written, "pestpp-glm exited without writing a run record"
    assert shutil.which("pestpp-glm") or pest_executable("pestpp-glm")
