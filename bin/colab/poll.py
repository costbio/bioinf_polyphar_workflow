"""Colab exec 3 of 4: report the detached worker's progress, and keep the VM alive.

Run with `colab exec -f bin/colab/poll.py --timeout 300` as often as needed. Short and
cheap on purpose: it is both the progress read and the activity that holds the session.
It must stay fast -- a poll that blocks is a poll that stops being a keep-alive.
"""
import glob
import os
import subprocess

BUNDLE = "/content/bundle"
LOG = "/content/worker.log"
PIDFILE = "/content/worker.pid"

pid = ""
if os.path.exists(PIDFILE):
    pid = open(PIDFILE).read().strip()
if pid:
    rc = subprocess.run(["bash", "-c", f"kill -0 {pid} 2>/dev/null"]).returncode
    alive = "RUNNING" if rc == 0 else "GONE"
else:
    alive = "NO-PIDFILE"

done = len(glob.glob(f"{BUNDLE}/results/**/*.done", recursive=True))
# the worker writes r<i>.sdf per replicate, never .csv -- counting .csv made
# result_files=0 on a run that was producing results normally
sdfs = len(glob.glob(f"{BUNDLE}/results/**/r*.sdf", recursive=True))
try:
    with open(f"{BUNDLE}/units.tsv") as fh:
        total = sum(1 for _ in fh) - 1
except OSError:
    total = "?"

print(f"STATE={alive} pid={pid or '-'}  done={done}/{total}  result_sdf={sdfs}")
print(subprocess.run(["bash", "-c", f"tail -c 3000 {LOG} 2>/dev/null"],
                     capture_output=True, text=True).stdout)
