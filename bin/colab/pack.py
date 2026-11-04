"""Colab exec 4 of 4: pack only the results for the trip home.

Run with `colab exec -f bin/colab/pack.py --timeout 600`, then download
/content/results.tgz and run `polyphar dock-merge --bundle <the same bundle>` on the host.

Only results/ travels back. The receptors and the library came from the run folder and
merging re-reads them there, so sending them again would just be slower.
"""
import glob
import os
import subprocess

BUNDLE = "/content/bundle"
OUT = "/content/results.tgz"

done = sorted(glob.glob(f"{BUNDLE}/results/*.done"))
csvs = sorted(glob.glob(f"{BUNDLE}/results/*.csv"))
print(f"results: {len(csvs)} csv, {len(done)} done", flush=True)
if not done:
    raise SystemExit("nothing finished - do not pack an empty bundle")

p = subprocess.run(["tar", "czf", OUT, "-C", BUNDLE, "results"],
                   capture_output=True, text=True)
if p.returncode:
    raise SystemExit(f"tar failed: {p.stderr[-2000:]}")
print(f"{OUT}  {os.path.getsize(OUT)/1e6:.2f} MB", flush=True)
print(subprocess.run(["bash", "-c", f"tar tzf {OUT} | head -20"],
                     capture_output=True, text=True).stdout, flush=True)
