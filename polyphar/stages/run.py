"""`polyphar run`: the stages ids -> uniprot -> structures -> fix -> pockets -> similarity ->
select -> dock in order (or the contiguous range --from .. --to), each through the same ETA
gate as its single-stage subcommand (cli.gate): exact units, ETA printed and notified before
the first unit, milestone progress, done / FAILED.

* `run` accepts the union of the stages' own options (built from each stage's
  add_arguments(); two stages defining one option differently is an error).
* ids: runs when --txt / --xlsx is given; without them an existing 00_ids/ids.txt counts
  as done (skipped), and a missing one is an error before anything starts.
* A stage without calibration on this machine is benchmarked first (`benchmark --stages
  <stage>`), except pockets: DoGSiteScorer is host-only (bin/polyphar-docker calibrates it
  before `run`); uncalibrated, its ETA comes with the 5 % milestone.
* The first failed stage stops the run (priority-1 notification). A stage that finishes with
  counts["failed"] > 0, or dock with counts["incomplete_targets"] > 0, does not stop it but is
  flagged in the final summary (printed, notified, appended to logs/run.log).
* --dry-run prints the `polyphar estimate` table and does nothing else.

Heavy imports stay inside functions (cli.build_parser imports this module at parse time).
"""

from __future__ import annotations

import argparse
from datetime import datetime

from polyphar.benchmark import BENCH_STAGES, ORDER, estimate_rows, format_table, run_benchmark
from polyphar.common import WorkUnits
from polyphar.eta import Calibration, format_duration, polyphar_home

IDS_SOURCE = ("xlsx", "txt", "sheet", "id_column", "reviewed_column")
FLAG_KEYS = ("failed", "incomplete_targets")


# --- options -------------------------------------------------------------------------------

def _signature(action) -> tuple:
    return (type(action).__name__, tuple(sorted(action.option_strings)), action.dest,
            action.nargs, repr(action.const), repr(action.default), action.type,
            tuple(action.choices) if action.choices is not None else None)


def merge_stage_options(parser, stage_hooks) -> None:
    """Add the options of every (stage name, add_arguments) in `stage_hooks` to `parser`,
    one argument group per stage. An option (or dest) defined identically by several stages
    is added once; defined differently, it raises ValueError. Mutually exclusive groups stay
    exclusive but are never required here (each stage validates its own inputs)."""
    owners = {}   # option string or "dest:<dest>" -> (stage, signature)
    for action in parser._actions:
        for key in [*action.option_strings, f"dest:{action.dest}"]:
            owners[key] = ("run", _signature(action))
    for stage, hook in stage_hooks:
        scratch = argparse.ArgumentParser(add_help=False)
        hook(scratch)
        exclusive = {id(a): g for g in scratch._mutually_exclusive_groups
                     for a in g._group_actions}
        group = parser.add_argument_group(f"{stage} options")
        new_exclusive = {}
        for action in scratch._actions:
            sig = _signature(action)
            keys = [*action.option_strings, f"dest:{action.dest}"]
            clash = [(k, owners[k]) for k in keys if k in owners]
            if clash:
                for key, (owner, owner_sig) in clash:
                    if owner_sig != sig:
                        raise ValueError(
                            f"polyphar run: option {key.removeprefix('dest:')} is defined "
                            f"differently by stages {owner!r} and {stage!r}")
                continue
            g = exclusive.get(id(action))
            if g is not None and id(g) not in new_exclusive:
                new_exclusive[id(g)] = group.add_mutually_exclusive_group()
            target = group if g is None else new_exclusive[id(g)]
            target._add_action(action)
            for key in keys:
                owners[key] = (stage, sig)


def add_arguments(parser) -> None:
    from polyphar import cli

    parser.add_argument("--from", dest="from_stage", choices=ORDER, default=ORDER[0],
                        metavar="STAGE", help=f"first stage (default {ORDER[0]}); one of "
                                              f"{', '.join(ORDER)}")
    parser.add_argument("--to", dest="to_stage", choices=ORDER, default=ORDER[-1],
                        metavar="STAGE", help=f"last stage (default {ORDER[-1]})")
    hooks = []
    for name in ORDER:
        module = cli._load_stage(name)
        if module is not None and hasattr(module, "add_arguments"):
            hooks.append((name, module.add_arguments))
    merge_stage_options(parser, hooks)


# --- run -----------------------------------------------------------------------------------

def units(layout, args) -> WorkUnits:
    return WorkUnits(1, "run")


def run(layout, args, progress) -> dict:   # unused: cli.main calls main()
    return {"failed": 0} if main(layout, args) == 0 else {"failed": 1}


def _ensure_calibrated(stage: str) -> None:
    """benchmark --stages <stage> when this machine has no calibration for it. Never
    for pockets (host-calibrated) or trivial stages; a failed benchmark only leaves the
    stage uncalibrated."""
    if stage not in BENCH_STAGES or stage == "pockets":
        return
    calib = Calibration.load(polyphar_home())
    if calib.cost(stage) is not None:
        return
    print(f"{stage}: not calibrated on this machine; running polyphar benchmark "
          f"--stages {stage}", flush=True)
    try:
        run_benchmark([stage], calib, log=lambda m: print(m, flush=True))
    except Exception as exc:
        print(f"{stage}: benchmark failed ({type(exc).__name__}: {exc}); continuing "
              "uncalibrated", flush=True)


def _fmt_s(seconds) -> str:
    if seconds is None:
        return "uncalibrated"
    return "instant" if seconds == 0 else format_duration(seconds)


def _flags(rows) -> list[str]:
    """Failed stages and finished stages with failed units / incomplete targets."""
    out = []
    for r in rows:
        if r["status"] == "FAILED":
            out.append(f"{r['stage']} FAILED")
            continue
        for key in FLAG_KEYS:
            v = (r.get("counts") or {}).get(key)
            if isinstance(v, (int, float)) and v > 0:
                out.append(f"{r['stage']} {key}={v:g}")
    return out


def summary_table(rows, run_name: str) -> str:
    lines = [f"run summary ({run_name})",
             f"{'stage':<12}{'status':<16}{'units':<22}{'estimate':<14}{'actual':<12}"
             f"{'failed':<8}incomplete_targets"]
    for r in rows:
        counts = r.get("counts")
        failed = incomplete = "-"
        if counts is not None:
            failed = str(counts.get("failed", 0))
            if r["stage"] == "dock":
                incomplete = str(counts.get("incomplete_targets", 0))
        est = _fmt_s(r["estimate_s"]) if "estimate_s" in r else "-"
        actual = format_duration(r["elapsed_s"]) if r.get("elapsed_s") is not None else "-"
        lines.append(f"{r['stage']:<12}{r['status']:<16}{r.get('units', '-'):<22}{est:<14}"
                     f"{actual:<12}{failed:<8}{incomplete}")
    flags = _flags(rows)
    lines.append(f"failures: {'; '.join(flags) if flags else 'none'}")
    return "\n".join(lines)


def _finish(layout, rows, notifier, title, priority) -> str:
    from polyphar.notify import run_name

    text = summary_table(rows, run_name(layout.run_dir))
    print(text, flush=True)
    try:
        layout.logs.mkdir(parents=True, exist_ok=True)
        with open(layout.logs / "run.log", "a", encoding="utf-8") as fh:
            fh.write(f"{datetime.now().astimezone().isoformat(timespec='seconds')} "
                     f"{title}\n{text}\n")
    except OSError as exc:
        print(f"could not write logs/run.log: {exc}", flush=True)
    try:
        notifier(title, text, priority=priority)
    except Exception:
        pass   # notifications never fail a run
    return text


def main(layout, args) -> int:
    from polyphar import cli

    if args.dry_run:
        rows, observed = estimate_rows(layout, args.nprocs, Calibration.load(polyphar_home()))
        print(format_table(rows, observed))
        return 0
    first, last = ORDER.index(args.from_stage), ORDER.index(args.to_stage)
    if first > last:
        raise SystemExit(f"polyphar run: --from {args.from_stage} comes after "
                         f"--to {args.to_stage}")
    stages = ORDER[first:last + 1]
    ids_given = any(getattr(args, k, None) for k in IDS_SOURCE)
    if "ids" in stages and not ids_given and not layout.ids_txt.is_file():
        raise SystemExit("polyphar run: 00_ids/ids.txt does not exist; give the ids source "
                         "(--txt or --xlsx with --sheet and --id-column)")

    notifier = cli.make_notifier(layout, args)
    rows = []
    for i, name in enumerate(stages):
        if name == "ids" and not ids_given:
            print("ids: 00_ids/ids.txt exists and no --txt/--xlsx given: skipped (done)",
                  flush=True)
            rows.append({"stage": name, "status": "skipped (done)"})
            continue
        _ensure_calibrated(name)
        row = {"stage": name, "status": "running"}
        rows.append(row)
        try:
            row["counts"] = cli.gate(name, cli._load_stage(name), layout, args, notifier,
                                     record=row)
        except BaseException:
            row["status"] = "FAILED"
            rows.extend({"stage": s, "status": "not run"} for s in stages[i + 1:])
            _finish(layout, rows, notifier, f"run FAILED at {name}", priority=1)
            raise
        row["status"] = "done"
    flags = _flags(rows)
    title = f"run done with failures ({len(flags)} flagged)" if flags else "run done"
    _finish(layout, rows, notifier, title, priority=0)
    return 0
