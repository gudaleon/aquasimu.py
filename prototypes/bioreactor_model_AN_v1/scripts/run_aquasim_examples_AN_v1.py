"""Fetch the public AQUASIM 2.1g source + tutorial examples from Eawag, build the BATCH
version on Linux, run the examples used as consistency benchmarks and write the listings to
benchmarks_AN_v1/aquasim_out/.

Source of everything: https://www.eawag.ch/en/department/siam/software/
  aquasim-source_code.zip, tutor20.zip (AQUASIM 2.0 tutorial example files)
Only the *listing outputs* produced by running the real program are stored in this project;
the Eawag input/source files are downloaded on demand into /workspace/aquasim_examples_AN.

Modifications made to the Eawag example files (clearly: NOT to the physics):
  * listing precision (PLOTLISOPT) raised from 4 to 10 significant digits
  * soil_a/soil_b "_g502": number of grid points of the soil column raised from 52 to 502 (grid-refinement run)
  * proc_b: initial conditions changed to C_A(0)=0, C_B(0)=1 so that the Monod step B -> C is
    isolated (file procB_mm_mod.aqu)
"""
import os, re, subprocess, shutil, sys, urllib.request, zipfile, pathlib

BASE = "https://www.eawag.ch/fileadmin/Domain1/Abteilungen/siam/software/aquasim/"
WORK = pathlib.Path("/workspace/aquasim_examples_AN")
HERE = pathlib.Path(__file__).resolve().parent.parent
OUT = HERE / "benchmarks_AN_v1" / "aquasim_out"


def fetch(name):
    WORK.mkdir(exist_ok=True)
    f = WORK / name
    if not f.exists():
        urllib.request.urlretrieve(BASE + name, f)
    return f


def build():
    exe = WORK / "build_aqb" / "aquasimb"
    if exe.exists():
        return exe
    z = fetch("aquasim-source_code.zip"); ini = fetch("aquasim_program.zip")
    src = WORK / "x_aquasim-source_code"
    zipfile.ZipFile(z).extractall(src)
    b = WORK / "build_aqb"; b.mkdir(exist_ok=True)
    for d in ("aqcore", "batch"):
        for p in (src / d).iterdir():
            shutil.copy(p, b)
    zipfile.ZipFile(ini).extract("aquasim.ini", b)
    mk = (b / "mkaquab.lnx").read_text().replace("CFLAGS =\n", "CFLAGS = -fpermissive -w -std=gnu++98 -fcommon\n", 1)
    (b / "mkaquab.lnx").write_text(mk)
    subprocess.run(["make", "-f", "mkaquab.lnx", "-j8"], cwd=b, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return exe


def prep(run, name, modify=None):
    raw = (run / "tutor20" / f"{name}.aqu").read_text(errors="ignore")
    if modify:
        raw = modify(raw)
    head, body = raw.split("{AQUASYS}", 1)          # (line breaks inside tokens are not significant)
    body = "{AQUASYS}" + body.replace("\n", "")
    body = body.replace("{PLOTLISOPT}{{1}{4}{TAB}}", "{PLOTLISOPT}{{1}{10}{TAB}}")
    assert "{PLOTLISOPT}{{1}{10}{TAB}}" in body
    return head + body + "\n"


def mod_procB(raw):
    head, body = raw.split("{AQUASYS}", 1)
    body = "{AQUASYS}" + body.replace("\n", "")
    body = body.replace("{{0}{C_A}{C_Aini}}", "{{0}{C_A}{C_Aini}{0}{C_B}{C_Bini}}")
    body = body.replace("{FORMVAR}{{1}{k_A}", "{FORMVAR}{{1}{C_Bini}{Initial B}{mg/l}{1}}{FORMVAR}{{1}{k_A}")
    body = body.replace("{{1}{C_Aini}{Initial concentration of substance A}{mg/l}{1}}",
                        "{{1}{C_Aini}{Initial concentration of substance A}{mg/l}{0}}")
    return head + body + "\n"


def mod_soil_precision(raw):
    return raw


def main():
    exe = build()
    z = fetch("tutor20.zip")
    run = WORK / "run_AN_v1"; run.mkdir(exist_ok=True)
    zipfile.ZipFile(z).extractall(run)
    shutil.copy(exe.parent / "aquasim.ini", run / "aquasim.ini")
    OUT.mkdir(parents=True, exist_ok=True)
    cases = {  # file -> (aqu name to write, modify, plots)
        "proc_a": ("proc_a", None, ["Conc"]),
        "procB_mm_mod": ("proc_b", mod_procB, ["Conc"]),
        "advdif_a": ("advdif_a", None, ["C_t"]),
        "soil_a": ("soil_a", None, ["plot"]),
        "soil_b": ("soil_b", None, ["plot"]),
        "soil_a_g502": ("soil_a", lambda r: r.replace("{TRUE}{52}{TRUE}", "{TRUE}{502}{TRUE}"), ["plot"]),
        "soil_b_g502": ("soil_b", lambda r: r.replace("{TRUE}{52}{TRUE}", "{TRUE}{502}{TRUE}"), ["plot"]),
        "film_a": ("film_a", None, ["LF_t", "S_z", "X_z"]),
        "film_b": ("film_b", None, ["LF_t", "S_z", "X_z"]),
    }
    for tag, (src, mod, plots) in cases.items():
        (run / f"{tag}.aqu").write_text(prep(run, src, mod))
        r = subprocess.run([str(exe), "-s", f"{tag}.s.log", f"{tag}.aqu", f"{tag}.res.aqu"], cwd=run, capture_output=True, text=True)
        if r.returncode != 0 or "Number of steps" not in (run / f"{tag}.s.log").read_text():
            print("AQUASIM simulation failed for", tag); continue
        for p in plots:
            lis = OUT / f"{tag}__{p}.txt"
            if lis.exists(): lis.unlink()
            (run / "l.cmd").write_text(f"{p} {lis}\n")
            subprocess.run([str(exe), "-l", f"{tag}.l.log", f"{tag}.res.aqu", "l.cmd"], cwd=run, capture_output=True)
        print("ok", tag)
    shutil.copy(run / "aquasim.ini", OUT / "aquasim_AN_v1_used.ini") if False else None


if __name__ == "__main__":
    main()
