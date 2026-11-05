"""Colab exec: pack the units the host does not have yet, for an INCREMENTAL pull.

Replaces the old pack.py, which tarred the whole results/ tree once, at the very end. That
made a reclaimed VM cost every unit drawn since the session started -- up to 200 units and
~20 GPU-hours on net1. The host now calls this every poll, so a lost VM costs at most one
poll interval.

Reads /content/have.txt, one `<label>/<rN>` per line: the units the host already holds.
Writes /content/outbox.tgz with the rest, and prints `FETCH units=<n> bytes=<n>`.

Selection is recomputed from scratch on every call and nothing is marked on the VM, so a
dropped download, a killed wrapper or a re-attach all self-correct on the next poll.

Only units carrying a `.done` marker travel. The worker writes that marker last and
atomically (dock.py `_dock_unit`), so a unit still being written is simply not selected,
and `unit_state` on the host would reject it anyway. Within each unit the marker is also
written into the archive LAST, so an interrupted extraction cannot leave a certificate
without the record it certifies.
"""
import os
import tarfile

# The paths are overridable so the script can be exercised against a temporary directory
# in the tests; on the VM they are the real ones and nothing sets these.
BUNDLE = os.environ.get("FETCH_BUNDLE", "/content/bundle")
RESULTS = os.path.join(BUNDLE, "results")
HAVE = os.environ.get("FETCH_HAVE", "/content/have.txt")
OUT = os.environ.get("FETCH_OUT", "/content/outbox.tgz")

# Per unit, in the order they enter the archive. The marker is last on purpose.
SUFFIXES = (".sdf", ".log", ".prov.json", ".filtered.tsv", ".done")

# Small bundle-root files that pack.py never sent home, which is why a merged manifest
# recorded `engine_worker: "unknown"` and the worker's own unit failures were invisible on
# the host. They are re-sent whenever an archive is built; all are tiny.
EXTRAS = ("worker.json", "failed.tsv", "logs/dock-worker.log", "results/discarded.tsv")


def have():
    try:
        with open(HAVE, encoding="utf-8") as fh:
            return {ln.strip() for ln in fh if ln.strip()}
    except OSError:
        return set()


def units():
    """[(unit id, directory, stem)] for every unit with a .done marker, sorted."""
    out = []
    for label in sorted(os.listdir(RESULTS)) if os.path.isdir(RESULTS) else []:
        d = os.path.join(RESULTS, label, "replicates")
        if not os.path.isdir(d):
            continue
        for name in sorted(os.listdir(d)):
            if name.endswith(".done"):
                stem = name[: -len(".done")]
                out.append((f"{label}/{stem}", d, stem))
    return out


held = have()
new = [(uid, d, stem) for uid, d, stem in units() if uid not in held]
force = os.environ.get("FETCH_FORCE") == "1"

if not new and not force:
    # Nothing to send. Remove a stale archive so a later failed build cannot be mistaken
    # for a fresh one, and tell the host to skip the download.
    try:
        os.unlink(OUT)
    except OSError:
        pass
    print("FETCH units=0 bytes=0", flush=True)
    raise SystemExit(0)

tmp = OUT + ".tmp"
with tarfile.open(tmp, "w:gz") as tar:
    for _uid, d, stem in new:
        for suffix in SUFFIXES:
            p = os.path.join(d, stem + suffix)
            if os.path.exists(p):
                tar.add(p, arcname=os.path.relpath(p, BUNDLE))
    for rel in EXTRAS:
        p = os.path.join(BUNDLE, rel)
        if os.path.exists(p):
            tar.add(p, arcname=rel)
os.replace(tmp, OUT)

print(f"FETCH units={len(new)} bytes={os.path.getsize(OUT)}", flush=True)
for uid, _d, _stem in new[:20]:
    print(f"  {uid}", flush=True)
if len(new) > 20:
    print(f"  ... and {len(new) - 20} more", flush=True)
