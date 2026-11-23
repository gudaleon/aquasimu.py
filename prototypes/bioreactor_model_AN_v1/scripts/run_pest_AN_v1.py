"""Run PEST++ (pestpp-glm, then pestpp-ies) on the SYNTHETIC case, driven from Python with pyEMU.
PEST++ binaries: https://github.com/pestpp/pestpp/releases (v5.2.16 linux build)."""
import sys, os, shutil, subprocess, pathlib, time
HERE = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "scripts"))
import pyemu
import setup_pest_AN_v1 as setup

BIN = pathlib.Path(os.environ.get("PESTPP_BIN", "/workspace/pestpp_bin_AN/bin"))
P = HERE / "pest_AN_v1"
FILES = ["forward_run_AN_v1.py", "params_AN_v1.tpl", "params_AN_v1.dat", "output_AN_v1.ins", "pest_AN_v1.pst"]
os.environ["PATH"] = str(pathlib.Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]   # 'python' with numba

setup.build()
# ---- 1. pestpp-glm: Gauss-Levenberg-Marquardt + linear (FOSM) uncertainty ---------------------------
if "--ies-only" not in sys.argv:
    g = P / "run_glm"
    shutil.rmtree(g, ignore_errors=True); g.mkdir()
    for f in FILES: shutil.copy(P / f, g / f)
    t = time.time()
    subprocess.run([str(BIN / "pestpp-glm"), "pest_AN_v1.pst"], cwd=g, check=True, stdout=open(g / "stdout.txt", "w"))
    print(f"pestpp-glm finished in {time.time()-t:.0f} s")

# ---- 2. pestpp-ies: iterative ensemble smoother (200 realisations, 4 iterations) ----------------------
w = P / "worker_template"
shutil.rmtree(w, ignore_errors=True); w.mkdir()
for f in FILES: shutil.copy(P / f, w / f)
pst = pyemu.Pst(str(w / "pest_AN_v1.pst"))
pst.control_data.noptmax = 4
pst.pestpp_options["ies_num_reals"] = 200
pst.pestpp_options["ies_subset_size"] = 8
pst.pestpp_options["ies_lambda_mults"] = "0.5,1.0,2.0"
pst.pestpp_options["lambda_scale_fac"] = "0.7,1.0"
pst.pestpp_options.pop("uncertainty", None)
pst.write(str(w / "pest_AN_v1.pst"), version=1)
(P / "workers_tmp").mkdir(exist_ok=True)
m = P / "run_ies"
shutil.rmtree(m, ignore_errors=True)
t = time.time()
pyemu.os_utils.start_workers(str(w), str(BIN / "pestpp-ies"), "pest_AN_v1.pst", num_workers=8,
                             worker_root=str(P / "workers_tmp"), master_dir=str(m), port=4004, verbose=False)
print(f"pestpp-ies finished in {time.time()-t:.0f} s")
shutil.rmtree(P / "workers_tmp", ignore_errors=True); shutil.rmtree(w, ignore_errors=True)
