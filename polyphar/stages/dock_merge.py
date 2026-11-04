"""`polyphar dock-merge --bundle DIR`: fold a worked bundle's results into the run folder.

It copies each finished unit's record (replicates/r<i>.sdf and its marker) into
07_docking/<label>/, then re-merges and rewrites the manifest exactly as `dock` would --
without running the engine, because the merging host need not have one.

Two rules decide what is taken, and both come from the replicate design:

* a result is taken only when the current plan has a unit with the same label, the same
  replicate index AND the same key. A bundle made under other settings describes different
  draws, and mixing them silently is what this check exists to stop. Rejections are
  reported, not swallowed.
* a unit already recorded here is never overwritten by the bundle's copy. The draw on
  record stays the draw on record -- and the bundle's losing draw is written to
  discarded.tsv, because dropping a real sample is a selection effect whether or not
  anything failed.

The engine version for the KEYS comes from the bundle's run_key.json, not from this host:
the keys were decided at export and must be compared as they were computed. The engine
that actually DREW the results is a different thing, read from the bundle's worker.json
and recorded in the merged manifest as `engine_worker` -- for a stage that exists because
the engine is on another machine, the Methods source must not name the planning host's
binary.

Before it imports anything, the merge compares the bundle's run_key.json `options` with
its own arguments and with the run folder's manifest, and WARNS by name on a mismatch. It
does not refuse: a partly foreign bundle must still merge the units that do match, so the
per-unit key check stays the mechanism. The warning exists because a bundle whose options
differ produces N opaque key rejections and nothing saying which flag caused them.
"""

from __future__ import annotations

import sys

from polyphar.stages import dock

STAGE = "dock-merge"

# Options that define what a draw IS. A difference here means the bundle's units are not
# this run's units. `replicates` is deliberately not among them: how MANY draws a run wants
# is not part of what one draw is (it is not in the unit key either), so a 3-replicate
# bundle merging into a 5-replicate run is legitimate and only warned about.
KEY_OPTIONS = ("base_seed", "exhaustiveness", "num_modes", "ligand_addH", "engine_name",
               "search_mode", "batch_size")
FLAG_FOR = {"base_seed": "--seed", "replicates": "--replicates",
            "engine_name": "--engine", "search_mode": "--search-mode",
            "batch_size": "--batch-size"}


def _warn_gpu_span(layout, progress, args) -> list[str]:
    """Warn when one target's replicates were drawn on different GPU models.

    Throughput and plausibly scores vary by device, so replicates split across an L4 and an
    A100 put a potential systematic difference INSIDE what is being averaged. The run is
    not refused -- the draws are real and the raw per-replicate records are kept -- but the
    Methods have to know, and nothing else in the pipeline would ever say so."""
    plan, _lib, _engine_v = dock.plan_for(layout, args)
    spans = []
    for label in sorted({u.label for u in plan.units}):
        models = {}
        for u in plan.units_of(label):
            model = (dock.read_provenance(layout, u).get("gpu_model") or "").strip()
            if model:
                models.setdefault(model, []).append(u.name)
        if len(models) > 1:
            spans.append(f"{label}: " + ", ".join(
                f"{m} ({', '.join(sorted(n))})" for m, n in sorted(models.items())))
    if not spans:
        return []
    message = (f"{len(spans)} target(s) have replicates drawn on DIFFERENT GPU models, so "
               f"a device difference sits inside the mean over replicates: "
               + "; ".join(spans[:5]) + (", ..." if len(spans) > 5 else "")
               + ". The raw per-replicate scores are kept, so this can be checked")
    print(f"{STAGE}: WARNING: {message}", file=sys.stderr, flush=True)
    notifier = getattr(progress, "notifier", None)
    if notifier is not None:
        try:
            notifier(f"{STAGE}: replicates span GPU models", message, priority=1)
        except Exception:
            pass
    return spans


def add_arguments(parser) -> None:
    # the same options as `dock`: merging has to plan the run exactly as `dock` would, or
    # the keys it compares the bundle against would not be the run's keys.
    dock.add_arguments(parser)
    parser.add_argument("--bundle", required=True,
                        help="bundle directory worked by `polyphar dock-worker`")
    parser.add_argument("--control", action="store_true", help=dock.CONTROL_FLAG_HELP)


def _bundle(args):
    # argparse makes --bundle required, so this only fires for a programmatic caller.
    if not getattr(args, "bundle", None):
        raise ValueError("dock-merge needs --bundle <directory>")
    return args.bundle


def _bundle_options(args) -> dict:
    return dict(dock.read_run_key(_bundle(args)).get("options") or {})


def _engine(args) -> str:
    return str(_bundle_options(args).get("engine") or "unavailable")


def _worker_engine(bundle) -> tuple[str | None, dict]:
    """(engine that actually ran the units, worker.json) from the bundle, or (None, {}).

    `dock-worker` has always written this and nothing ever read it, so the merged manifest
    named the engine the PLANNING host happened to have."""
    import json

    try:
        info = json.loads((dock.Path(bundle) / "worker.json").read_text())
    except (OSError, ValueError):
        return None, {}
    engine = info.get("engine")
    return (str(engine) if engine and engine != "unavailable" else None), info


def _check_options(layout, progress, args) -> None:
    """Warn, up front and by name, when the bundle's docking options are not this merge's.

    A difference in a draw-defining option means every unit will be rejected on its key:
    the per-unit rejection is the mechanism (a partly foreign bundle must still merge the
    units that do match), and this says WHICH flag caused it instead of leaving the user
    with N opaque key rejections. The softer differences -- replicate count, the run
    folder's own recorded arguments -- are warned about the same way."""
    bopts = _bundle_options(args)
    mine = dock._options(args, bopts.get("engine") or "unavailable")
    warn = []
    bad = [k for k in KEY_OPTIONS
           if k in bopts and str(bopts[k]) != str(mine.get(k))]
    if bad:
        detail = ", ".join(f"{k}: bundle {bopts[k]!r} vs this merge {mine.get(k)!r}"
                           f"{' (' + FLAG_FOR[k] + ')' if k in FLAG_FOR else ''}"
                           for k in bad)
        warn.append(f"the bundle {_bundle(args)} was made under different docking options, "
                    f"so its draws are NOT this run's draws and every unit will be "
                    f"rejected on its key -- {detail}. Re-run with the bundle's options, "
                    f"or merge it into the run folder it was exported from")
    if "replicates" in bopts and int(bopts["replicates"]) != int(mine["replicates"]):
        warn.append(f"the bundle was exported with --replicates {bopts['replicates']}, "
                    f"this merge plans {mine['replicates']}")
    recorded = (dock._read_manifest(layout).get("args") or {})
    for k, arg in (("base_seed", "seed"), ("replicates", "replicates")):
        if arg in recorded and str(recorded[arg]) != str(mine[k]):
            warn.append(f"this run folder's last manifest recorded {arg}="
                        f"{recorded[arg]!r}, this merge uses {mine[k]!r}")
    for w in warn:
        print(f"{STAGE}: WARNING: {w}", file=sys.stderr, flush=True)
    notifier = getattr(progress, "notifier", None)
    if warn and notifier is not None:
        try:
            notifier(f"{STAGE}: bundle options differ", "; ".join(warn), priority=1)
        except Exception:
            pass


def units(layout, args):
    """Unit records in the bundle that are not yet recorded here (one unit = one unit of
    work: merging is a copy, not a docking run)."""
    from polyphar.common import WorkUnits

    bundle = _bundle(args)
    args.engine_version = _engine(args)
    layout, control_lib = dock.control_view(layout, args)
    plan, _lib, _engine_v = dock.plan_for(layout, args, library=control_lib)
    bl = dock.BundleLayout(bundle)
    want = {(u.label, u.replicate): u for u in plan.units}
    n = 0
    for b in dock.read_bundle_units(bundle):
        mine = want.get((b.unit.label, b.unit.replicate))
        if (mine is not None and mine.key == b.unit.key
                and dock.unit_state(bl, b.unit) == dock.STATE_DONE
                and dock.unit_state(layout, mine) != dock.STATE_DONE):
            n += 1
    return WorkUnits(n, "unit")


def run(layout, args, progress) -> dict:
    bundle = _bundle(args)
    args.engine_version = _engine(args)
    layout, control_lib = dock.control_view(layout, args)
    _check_options(layout, progress, args)
    plan, _lib, _engine_v = dock.plan_for(layout, args, library=control_lib)
    # taken BEFORE import_bundle, which is what discards here: dock.run's own baseline is
    # taken after, so without this counts["discarded"] reported 0 while discarded.tsv grew.
    baseline = dock._count_discards(layout.docking)
    worker_engine, worker = _worker_engine(bundle)
    taken, rejected, dropped = dock.import_bundle(layout, plan, bundle)
    progress.advance(len(taken))
    # the engine is not run here: everything outstanding simply stays outstanding
    spans = _warn_gpu_span(layout, progress, args)
    counts = dock.run(
        layout, args, progress, execute=False, discard_baseline=baseline,
        extra_counts={"merged": len(taken), "rejected": len(rejected),
                      "dropped": len(dropped),
                      "targets_spanning_gpu_models": len(spans)},
        extra_manifest={"merged_from": str(bundle),
                        # the engine the keys were planned under vs the one that drew them
                        "engine_worker": worker_engine or "unknown",
                        "worker_versions": worker.get("versions") or {},
                        # the rest of what the bundle's worker recorded about its host:
                        # worker.json used to be written and never read at all
                        "worker_host": worker.get("host") or "unknown",
                        "worker_engine_name": worker.get("engine_name") or "unknown",
                        "worker_gpu_models": worker.get("gpu_models") or [],
                        "targets_spanning_gpu_models": spans})
    counts = dict(counts)
    for item, reason in rejected:
        from polyphar.common import append_failure

        append_failure(layout.docking, item, STAGE, reason)
    if worker_engine and worker_engine != _engine(args):
        print(f"{STAGE}: the bundle's units were drawn by {worker_engine!r}, while their "
              f"keys were planned under {_engine(args)!r} (manifest: engine_worker)",
              file=sys.stderr, flush=True)
    if dropped:
        print(f"{STAGE}: {len(dropped)} bundle draw(s) were dropped because this run "
              f"folder already has a draw for those units (see 07_docking/discarded.tsv)",
              file=sys.stderr, flush=True)
    if rejected:
        print(f"{STAGE}: {len(rejected)} unit(s) from {bundle} were NOT merged "
              f"(see 07_docking/failed.tsv)", file=sys.stderr, flush=True)
    return counts
