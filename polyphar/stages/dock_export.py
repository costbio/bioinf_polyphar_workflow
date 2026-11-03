"""`polyphar dock-export`: bundle the docking units that still have to be drawn.

For the case the single-host `dock` cannot cover: the engine is on a machine that does not
have the run folder (a Colab VM, TRUBA, a rented container). This stage plans here --
reads final_targets.csv, prepares the receptor PDBQTs, computes every (target, replicate)
unit key -- and writes the outstanding units, their receptors, their boxes and the titled
library into a bundle directory. `dock-worker` runs it, `dock-merge` folds it back in.

It never runs the engine and never needs one.

`--max-units N` takes only as much work as one session can plausibly finish, so a dropped
session loses at most the units in flight. Units are written replicate-major (r0 of every
target, then r1, ...), so what a lost session costs is one draw of many targets rather
than every draw of a few. The next export recomputes what is outstanding, which makes
progress monotone with no coordination between sessions.

**A unit already recorded in the run folder is never exported.** Its draw happened; a
second one would be a different sample, not a repeat.
"""

from __future__ import annotations

import time
from pathlib import Path

from polyphar.stages import dock

STAGE = "dock-export"
# defined in dock so RESERVED_DIRS can keep this directory out of the orphan sweep
EXPORT_DIRNAME = dock.EXPORT_DIRNAME


def add_arguments(parser) -> None:
    dock.add_arguments(parser)
    parser.add_argument("--bundle", default=None,
                        help="bundle directory to write (default "
                             "07_docking/export/<timestamp>/)")
    parser.add_argument("--max-units", type=dock._positive_int, default=None,
                        help="export at most this many (target, replicate) units "
                             "(default: all outstanding ones)")
    parser.add_argument("--engine-version", default=None,
                        help="engine version to put in the unit keys when this host has no "
                             "engine binary (e.g. 'Uni-Dock 1.2.0'); without it the local "
                             "binary decides")


def _bundle_dir(layout, args) -> Path:
    if getattr(args, "bundle", None):
        return Path(args.bundle)
    root = layout.docking / EXPORT_DIRNAME
    stamp, k = time.strftime("%Y%m%dT%H%M%S"), 1
    dest = root / stamp
    while dest.exists():
        k += 1
        dest = root / f"{stamp}_{k}"
    return dest


def units(layout, args):
    """Units to export, capped by --max-units. Counted in ligands, like `dock`."""
    from polyphar.common import WorkUnits

    plan, _lib, _engine = dock.plan_for(layout, args)
    n = len(plan.outstanding(layout))
    if getattr(args, "max_units", None):
        n = min(n, int(args.max_units))
    return WorkUnits(n * plan.n_ligands, "ligand")


def run(layout, args, progress) -> dict:
    plan, lib, engine = dock.plan_for(layout, args)
    rows: list[tuple[str, str]] = list(plan.missing)
    prepared = dock.prepare_receptors(layout, plan, rows)
    ligands_sdf = dock.titled_library(layout, lib, plan)
    dest = _bundle_dir(layout, args)
    dest.mkdir(parents=True, exist_ok=True)
    written = dock.write_bundle(layout, args, dest, plan, prepared, ligands_sdf,
                                max_units=getattr(args, "max_units", None))
    for item, reason in rows:
        from polyphar.common import append_failure

        append_failure(dest, item, STAGE, reason)
    print(f"{STAGE}: bundle {dest}", flush=True)
    progress.advance(len(written) * plan.n_ligands)
    return {"units": len(written), "outstanding": len(plan.outstanding(layout)),
            "targets": len(prepared), "ligands": plan.n_ligands,
            "engine": engine, "bundle": str(dest), "failed": len(rows)}
