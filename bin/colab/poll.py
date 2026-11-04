"""Colab exec 3 of 4: report the detached worker's progress, and keep the VM alive.

Run with `colab exec -f bin/colab/poll.py --timeout 300` as often as needed. Short and
cheap on purpose: it is both the progress read and the activity that holds the session.
"""
import glob
import subprocess

BUNDLE = "/content/bundle"
LOG = "/content/worker.log"

alive = subprocess.run(
    ["bash", "-c", "pgrep -f 'dock-worker' >/dev/null && echo RUNNING || echo GONE"],
    capture_output=True, text=True).stdout.strip()
done = len(glob.glob(f"{BUNDLE}/results/*.done"))
csvs = len(glob.glob(f"{BUNDLE}/results/*.csv"))
try:
    with open(f"{BUNDLE}/units.tsv") as fh:
        total = sum(1 for _ in fh) - 1
except OSError:
    total = "?"
print(f"STATE={alive}  done={done}/{total}  result_files={csvs}")
print(subprocess.run(["bash", "-c", f"tail -c 3000 {LOG} 2>/dev/null"],
                     capture_output=True, text=True).stdout)
