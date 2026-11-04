"""Colab exec 2 of 4: unpack a dock-export bundle and start `polyphar dock-worker`
DETACHED on the VM.

Run with `colab exec -f bin/colab/launch.py --timeout 600` after uploading bundle.tgz.

Detaching is the point. The kernel websocket is fragile -- the validation session on
2026-11-03 lost its connection while the work was still running -- so the job must not be
a child of the kernel. Progress is read afterwards with short `poll.py` execs, which
double as the liveness ping that stops Colab reclaiming the VM.
"""
import os
import shutil
import subprocess

BUNDLE_TGZ = "/content/bundle.tgz"
BUNDLE = "/content/bundle"
LOG = "/content/worker.log"
UD_BIN = "/content/ud/bin"

if not os.path.exists(BUNDLE_TGZ):
    raise SystemExit(f"{BUNDLE_TGZ} not found - upload it first")

shutil.rmtree(BUNDLE, ignore_errors=True)
os.makedirs(BUNDLE, exist_ok=True)
p = subprocess.run(["tar", "xzf", BUNDLE_TGZ, "-C", BUNDLE, "--strip-components=1"],
                   capture_output=True, text=True)
if p.returncode:
    raise SystemExit(f"untar failed: {p.stderr[-2000:]}")
print("extracted:", sorted(os.listdir(BUNDLE)), flush=True)

for need in ("units.tsv", "run_key.json", "library.sdf", "receptors"):
    if not os.path.exists(os.path.join(BUNDLE, need)):
        raise SystemExit(f"bundle is missing {need}")
os.makedirs(os.path.join(BUNDLE, "results"), exist_ok=True)

with open(os.path.join(BUNDLE, "units.tsv")) as fh:
    n_units = sum(1 for _ in fh) - 1
print(f"units to work: {n_units}", flush=True)

# unidock_tools calls the bare name "unidock", so the prefix must be on the detached
# process's own PATH -- inheriting this exec's environment is not enough.
env_path = f"{UD_BIN}:" + os.environ.get("PATH", "")
cmd = (f"cd {BUNDLE} && PATH={env_path} nohup polyphar dock-worker --bundle {BUNDLE} "
       f"> {LOG} 2>&1 & echo $!")
p = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
print("launched pid:", p.stdout.strip(), p.stderr.strip(), flush=True)

print(subprocess.run(["bash", "-c", f"sleep 20; tail -20 {LOG}"],
                     capture_output=True, text=True).stdout, flush=True)
