"""Colab exec 2 of 4: unpack a dock-export bundle and start `polyphar dock-worker`
DETACHED on the VM.

Run with `colab exec -f bin/colab/launch.py --timeout 600` after uploading bundle.tgz.

Detaching is the point. The kernel websocket is fragile -- the validation session on
2026-11-03 lost its connection while the work was still running -- so the job must not be
a child of the kernel. Progress is read afterwards with short `poll.py` execs, which
double as the liveness ping that stops Colab reclaiming the VM.
"""
import hashlib
import os
import pathlib
import shlex
import shutil
import subprocess
import time

BUNDLE_TGZ = "/content/bundle.tgz"
BUNDLE = "/content/bundle"
LOG = "/content/worker.log"
PIDFILE = "/content/worker.pid"
UD_BIN = "/content/ud/bin"

PUSHOVER_ENV = "/content/pushover.env"
WANT = ("PUSHOVER_APP_TOKEN", "PUSHOVER_USER_KEY")
exports, keys = "", []
if os.path.exists(PUSHOVER_ENV):
    with open(PUSHOVER_ENV) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            k, v = k.strip(), v.strip().strip("'\"")
            if k in WANT and v:
                exports += f"{k}={shlex.quote(v)} "
                keys.append(k)

# Refuse to start silent. A detached worker that cannot notify runs for up to 20 hours and
# reports nothing: it can finish, stall or die and the first you hear of it is when you
# poll. That is a worse failure than not starting, and the VM is fresh every session, so
# the upload is exactly the step that gets forgotten on session 20 of 40.
if sorted(keys) != sorted(WANT):
    if os.environ.get("ALLOW_SILENT") == "1":
        print(f"pushover: NOT configured ({', '.join(keys) or 'no keys'}) - starting "
              f"silent because ALLOW_SILENT=1", flush=True)
    else:
        raise SystemExit(
            f"refusing to start a silent worker: {PUSHOVER_ENV} "
            + ("is missing" if not os.path.exists(PUSHOVER_ENV)
               else f"has only {', '.join(keys)} - it needs both "
                    f"{' and '.join(WANT)}" if keys
               else f"has neither {' nor '.join(WANT)}")
            + ".\n\nUpload the credentials first:\n"
            "    colab upload ~/.polyphar/pushover.env /content/pushover.env -s <session>\n\n"
            "Then run this again. To start without notifications anyway, pass\n"
            "    colab exec -f bin/colab/launch.py --env ALLOW_SILENT=1 ...")
else:
    print("pushover:", ", ".join(sorted(keys)), flush=True)

if not os.path.exists(BUNDLE_TGZ):
    raise SystemExit(f"{BUNDLE_TGZ} not found - upload it first")

# Refuse to start a second worker on this bundle. The VM is reused when the session is
# reused, and on 2026-11-04 an earlier worker was still running when a NEW bundle was
# uploaded over it: two workers then docked into one results/ tree, the older one working
# a units.tsv that no longer existed, writing draws whose keys match nothing. Never two.
def _alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False

running = subprocess.run(["bash", "-c", "pgrep -f 'polyphar dock-worker' || true"],
                         capture_output=True, text=True).stdout.split()
running = [q for q in running if _alive(q)]
if running:
    raise SystemExit(
        f"a dock-worker is already running on this VM (pid {', '.join(running)}).\n"
        f"Starting a second one would put two workers in {BUNDLE}/results with different\n"
        f"unit keys. Either let it finish, or stop this session and start a fresh one:\n"
        f"    colab stop -s <session>")

# Keep finished units across a RELAUNCH OF THE SAME BUNDLE -- that is the resume path, and
# wiping it cost a finished 1,960-ligand unit on 2026-11-04. But a DIFFERENT bundle's
# results are not ours: their unit keys belong to a plan this bundle does not contain, and
# dock-merge would reject them one opaque rejection at a time.
def _bundle_id():
    try:
        return hashlib.sha256(
            pathlib.Path(BUNDLE, "run_key.json").read_bytes()).hexdigest()
    except OSError:
        return None

was = _bundle_id()
os.makedirs(BUNDLE, exist_ok=True)
p = subprocess.run(["tar", "xzf", BUNDLE_TGZ, "-C", BUNDLE, "--strip-components=1"],
                   capture_output=True, text=True)
if p.returncode:
    raise SystemExit(f"untar failed: {p.stderr[-2000:]}")
print("extracted:", sorted(os.listdir(BUNDLE)), flush=True)

now = _bundle_id()
if was is not None and was != now:
    shutil.rmtree(os.path.join(BUNDLE, "results"), ignore_errors=True)
    print("a different bundle was here: cleared its results/", flush=True)
os.makedirs(os.path.join(BUNDLE, "results"), exist_ok=True)

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
