"""Colab exec 1 of 4: install the engine and polyphar on the VM. Idempotent.

Run with `colab exec -f bin/colab/setup.py --timeout 3600`.

Uni-Dock comes from conda-forge via micromamba because Colab has no Docker daemon; this
is the install path validated on an L4 on 2026-11-03. polyphar declares no dependencies,
so the worker's own imports (pandas, rdkit) are installed explicitly here.

POLYPHAR_COMMIT pins the worker to the commit that exported the bundle. Do not relax it
to a branch: the unit keys record the code version, and a worker running different code
produces draws that `dock-merge` cannot vouch for.
"""
import os
import subprocess
import sys
import time

POLYPHAR_REPO = "https://github.com/costbio/bioinf_polyphar_workflow"
POLYPHAR_COMMIT = os.environ.get("POLYPHAR_COMMIT", "").strip()

UD_PREFIX = "/content/ud"


def sh(cmd, timeout=3600, check=True):
    print(f"\n$ {cmd}", flush=True)
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    out, err = (p.stdout or "")[-4000:], (p.stderr or "")[-2000:]
    if out.strip():
        print(out, flush=True)
    if err.strip():
        print("STDERR:", err, flush=True)
    if check and p.returncode != 0:
        raise SystemExit(f"FAILED (rc={p.returncode}): {cmd}")
    return p.returncode


t0 = time.time()
print("=" * 72, flush=True)
sh("nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv", timeout=120)
print(f"python: {sys.version}", flush=True)

# --- Uni-Dock -------------------------------------------------------------
if not os.path.exists(f"{UD_PREFIX}/bin/unidock"):
    sh("cd /content && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest "
       "| tar -xj bin/micromamba", timeout=600)
    sh(f"cd /content && MAMBA_ROOT_PREFIX=/content/mamba ./bin/micromamba create -y "
       f"-p {UD_PREFIX} -c conda-forge unidock", timeout=3600)
else:
    print("unidock already installed", flush=True)

os.environ["PATH"] = f"{UD_PREFIX}/bin:" + os.environ.get("PATH", "")
sh("unidock --version 2>&1 | head -5", timeout=120, check=False)

# unidock_tools hardcodes cmd = ["unidock"], with no shutil.which, so the binary has to be
# on PATH of whatever process calls it -- launch.py re-exports this prefix for that reason.
sh(f"{sys.executable} -m pip install -q unidock-tools", timeout=1200)
sh("which unidocktools", timeout=180)

# --- polyphar -------------------------------------------------------------
sh(f"{sys.executable} -m pip install -q pandas rdkit requests", timeout=1800)
if not POLYPHAR_COMMIT:
    raise SystemExit("set POLYPHAR_COMMIT to the commit that exported the bundle")
sh(f"{sys.executable} -m pip install -q 'polyphar @ git+{POLYPHAR_REPO}@{POLYPHAR_COMMIT}'",
   timeout=1800)
sh(f"{sys.executable} -m polyphar.cli --help > /dev/null 2>&1 || polyphar --help | head -3",
   timeout=180, check=False)

print(f"\nsetup ok in {time.time()-t0:.0f}s", flush=True)
