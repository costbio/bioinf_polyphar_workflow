"""Colab exec 2 of 4: unpack a dock-export bundle and start `polyphar dock-worker`
DETACHED on the VM.

Run with `colab exec -f bin/colab/launch.py --timeout 600` after uploading bundle.tgz.

Detaching is the point. The kernel websocket is fragile -- the validation session on
2026-11-03 lost its connection while the work was still running -- so the job must not be
a child of the kernel. Progress is read afterwards with short `poll.py` execs, which
double as the liveness ping that stops Colab reclaiming the VM.
"""
import os
import shlex
import shutil
import subprocess
import time

BUNDLE_TGZ = "/content/bundle.tgz"
BUNDLE = "/content/bundle"
LOG = "/content/worker.log"
PIDFILE = "/content/worker.pid"
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

# Pushover, if the credentials were uploaded. They travel as a FILE rather than through
# `colab exec --env`: colab-cli records request bodies in ~/.config/colab-cli/colab.log at
# DEBUG level, and a token passed as an argument can end up there. Without this the
# detached worker runs silent -- it still calls notify(), but notify() finds no
# credentials, logs one warning and skips, and that log dies with the VM.
PUSHOVER_ENV = "/content/pushover.env"
exports = ""
if os.path.exists(PUSHOVER_ENV):
    keys = []
    with open(PUSHOVER_ENV) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            k, v = k.strip(), v.strip().strip("'\"")
            if k in ("PUSHOVER_APP_TOKEN", "PUSHOVER_USER_KEY"):
                exports += f"{k}={shlex.quote(v)} "
                keys.append(k)
    print("pushover:", ", ".join(sorted(keys)) or "file present but no usable keys", flush=True)
else:
    print("pushover: no /content/pushover.env - the worker will run silent", flush=True)

# unidock_tools calls the bare name "unidock", so the prefix must be on the detached
# process's own PATH -- inheriting this exec's environment is not enough.
env_path = f"{UD_BIN}:" + os.environ.get("PATH", "")
# setsid + all three streams detached from this kernel. `nohup ... &` alone is NOT
# enough: capture_output= gives bash a pipe, the backgrounded worker inherits it, and
# subprocess.run() then blocks until the WORKER exits -- which hangs the kernel for the
# whole job. That is worse than it sounds: poll.py and pack.py queue behind it, and the
# polls are also the keep-alive that stops Colab reclaiming the VM.
cmd = (f"cd {BUNDLE} && PATH={env_path} {exports}setsid polyphar dock-worker "
       f"--bundle {BUNDLE} < /dev/null > {LOG} 2>&1 & echo $! > {PIDFILE}")
subprocess.Popen(["bash", "-lc", cmd],
                 stdin=subprocess.DEVNULL,
                 stdout=subprocess.DEVNULL,
                 stderr=subprocess.DEVNULL,
                 start_new_session=True)

# the worker is now on its own; this exec must return promptly
for _ in range(30):
    time.sleep(1)
    if os.path.exists(PIDFILE):
        break
pid = open(PIDFILE).read().strip() if os.path.exists(PIDFILE) else "?"
alive = subprocess.run(["bash", "-c", f"kill -0 {pid} 2>/dev/null && echo RUNNING || echo GONE"],
                       capture_output=True, text=True).stdout.strip()
print(f"worker pid={pid} {alive}", flush=True)
print(subprocess.run(["bash", "-c", f"tail -20 {LOG} 2>/dev/null"],
                     capture_output=True, text=True).stdout, flush=True)
