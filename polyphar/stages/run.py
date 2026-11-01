"""`polyphar run`: the stages ids -> uniprot -> structures -> fix -> pockets -> similarity ->
select -> dock in order (or the contiguous range --from .. --to), each through the same ETA
gate as its single-stage subcommand (cli.gate): exact units, ETA printed and notified before
the first unit, milestone progress, done / FAILED.

* `run` accepts the union of the stages' own options (built from each stage's
  add_arguments(); two stages defining one option differently is an error).
* ids: runs when --txt / --xlsx is given; without them an existing 00_ids/ids.txt counts
  as done (skipped), and a missing one is an error before anything starts.
* Every stage's units() is computed first. A stage with 0 units and a manifest.json is
  "skipped (done)": not benchmarked, not run, no notifications (its manifest counts still
  feed the failure flags). A stage with 0 units but no manifest runs through the gate (it
  only finishes its bookkeeping).
* Options of stages outside --from..--to are warned about (stderr) and ignored.
* A stage with remaining units and no calibration on this machine is benchmarked first (`benchmark --stages
  <stage>`), except pockets: DoGSiteScorer is host-only (bin/polyphar-docker calibrates it
  before `run`); uncalibrated, its ETA comes with the 5 % milestone.
* The first failed stage stops the run (priority-1 notification). SIGTERM (docker stop) is
  turned into a failure for the whole run, so an interrupted benchmark or units() also ends
  in the "run FAILED" summary. A stage that finishes with counts["failed"] > 0, or dock with
  counts["incomplete_targets"] > 0, does not stop the run but is flagged in the final
  summary (printed, notified at priority 1, appended to logs/run.log).
* Exit codes: 0 success; 1 a stage failed (the exception propagates; 143 for SIGTERM);
  2 usage error; 3 the run finished but failures were flagged.
* --dry-run prints the `polyphar estimate` table and does nothing else.

Heavy imports stay inside functions (cli.build_parser imports this module at parse time).
"""

from __future__ import annotations

import argparse
import json
import signal
import sys
from datetime import datetime

from polyphar.benchmark import (BENCH_STAGES, ORDER, STAGE_DIRS, estimate_rows, format_table,
                                run_benchmark)
from polyphar.eta import Calibration, format_duration, polyphar_home

IDS_SOURCE = ("xlsx", "txt", "sheet", "id_column", "reviewed_column")
FLAG_KEYS = ("failed", "incomplete_targets")
EXIT_USAGE = 2
EXIT_FLAGGED = 3
SKIPPED = "skipped (done)"


# --- options -------------------------------------------------------------------------------

def _signature(action) -> tuple:
    return (type(action).__name__, tuple(sorted(action.option_strings)), action.dest,
            action.nargs, repr(action.const), repr(action.default), action.type,
            tuple(action.choices) if action.choices is not None else None)


def merge_stage_options(parser, stage_hooks) -> dict:
    """Add the options of every (stage name, add_arguments) in `stage_hooks` to `parser`,
    one argument group per stage. An option (or dest) defined identically by several stages
    is added once; defined differently, it raises ValueError. Mutually exclusive groups stay
    exclusive but are never required here (each stage validates its own inputs).
    Returns {dest: (option string, default, [stages defining it])}."""
    by_dest = {}
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
                if action.dest in by_dest:
                    by_dest[action.dest][2].append(stage)
                continue
            g = exclusive.get(id(action))
            if g is not None and id(g) not in new_exclusive:
                new_exclusive[id(g)] = group.add_mutually_exclusive_group()
            target = group if g is None else new_exclusive[id(g)]
            target._add_action(action)
            for key in keys:
                owners[key] = (stage, sig)
            if action.option_strings:
                by_dest[action.dest] = (action.option_strings[-1], action.default, [stage])
    return by_dest


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
    # which stage each option belongs to, for the out-of-range warning
    parser.set_defaults(run_option_owners=merge_stage_options(parser, hooks))


# --- run -----------------------------------------------------------------------------------

# `run` has no units()/run(): cli.main calls main() for it (hasattr(stage, "main")).


def _needs_benchmark(stage: str) -> bool:
    """No calibration for `stage` on this machine and it can be benchmarked here (never
    pockets: host-calibrated; never trivial stages)."""
    if stage not in BENCH_STAGES or stage == "pockets":
        return False
    return Calibration.load(polyphar_home()).cost(stage) is None


def _benchmark(stage: str) -> None:
    """benchmark --stages <stage>; a failed benchmark only leaves the stage uncalibrated."""
    print(f"{stage}: not calibrated on this machine; running polyphar benchmark "
          f"--stages {stage}", flush=True)
    try:
        run_benchmark([stage], Calibration.load(polyphar_home()),
                      log=lambda m: print(m, flush=True))
    except Exception as exc:
        print(f"{stage}: benchmark failed ({type(exc).__name__}: {exc}); continuing "
              "uncalibrated", flush=True)


def _fmt_s(seconds) -> str:
    if seconds is None:
        return "uncalibrated"
    return "instant" if seconds == 0 else format_duration(seconds)


def _flags(rows) -> list[str]:
    """Failed stages, and stages (run or skipped as done) with failed units / incomplete
    targets."""
    out = []
    for r in rows:
        if r["status"] == "FAILED":
            out.append(f"{r['stage']} FAILED")
            continue
        for key in FLAG_KEYS:
            v = (r.get("counts") or {}).get(key)
            if isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0:
                out.append(f"{r['stage']} {key}={v:g}")
    return out


def failures_line(rows) -> str:
    flags = _flags(rows)
    return f"failures: {'; '.join(flags) if flags else 'none'}"


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
    lines.append(failures_line(rows))
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
    # notify() truncates at 1024 characters: the failures line goes first
    try:
        notifier(title, failures_line(rows) + "\n" + text, priority=priority)
    except Exception:
        pass   # notifications never fail a run
    return text


def _usage(message: str) -> int:
    print(f"polyphar run: {message}", file=sys.stderr)
    return EXIT_USAGE


def _warn_out_of_range(args, stages) -> None:
    """One stderr line per given option whose stages are all outside the range."""
    for dest, (option, default, owners) in (getattr(args, "run_option_owners", None)
                                             or {}).items():
        if any(o in stages for o in owners):
            continue
        if getattr(args, dest, default) != default:
            print(f"polyphar run: warning: {option} belongs to {', '.join(owners)}, outside "
                  f"--from {args.from_stage} --to {args.to_stage}; ignored", file=sys.stderr)


def _manifest_counts(layout, stage) -> dict | None:
    """counts of the stage's manifest.json, or None when the stage has no manifest."""
    path = layout.run_dir / STAGE_DIRS[stage] / "manifest.json"
    if not path.is_file():
        return None
    try:
        counts = json.loads(path.read_text(encoding="utf-8")).get("counts")
        return counts if isinstance(counts, dict) else {}
    except (OSError, ValueError, AttributeError):
        return {}


def main(layout, args) -> int:
    from polyphar import cli

    if args.dry_run:
        rows, observed = estimate_rows(layout, args.nprocs, Calibration.load(polyphar_home()))
        print(format_table(rows, observed))
        return 0
    first, last = ORDER.index(args.from_stage), ORDER.index(args.to_stage)
    if first > last:
        return _usage(f"--from {args.from_stage} comes after --to {args.to_stage}")
    stages = ORDER[first:last + 1]
    ids_given = any(getattr(args, k, None) for k in IDS_SOURCE)
    if "ids" in stages and not ids_given and not layout.ids_txt.is_file():
        return _usage("00_ids/ids.txt does not exist; give the ids source "
                      "(--txt or --xlsx with --sheet and --id-column)")
    _warn_out_of_range(args, stages)

    notifier = cli.make_notifier(layout, args)
    previous = cli._install_sigterm_handler()   # docker stop -> SystemExit(143) -> summary
    try:
        return _run_stages(cli, layout, args, stages, ids_given, notifier)
    finally:
        if previous is not None:
            signal.signal(signal.SIGTERM, previous)


def _run_stages(cli, layout, args, stages, ids_given, notifier) -> int:
    rows = []
    for i, name in enumerate(stages):
        row = {"stage": name, "status": "running"}
        rows.append(row)
        try:
            if name == "ids" and not ids_given:
                print(f"ids: 00_ids/ids.txt exists and no --txt/--xlsx given: {SKIPPED}",
                      flush=True)
                row.update(status=SKIPPED, counts=_manifest_counts(layout, name))
                continue
            stage = cli._load_stage(name)
            units = stage.units(layout, args)
            done_counts = _manifest_counts(layout, name) if units.count == 0 else None
            if done_counts is not None:
                print(f"{name}: 0 units left and a manifest exists: {SKIPPED}", flush=True)
                row.update(status=SKIPPED, counts=done_counts,
                           units=f"0 ({units.cost_basis})")
                continue
            if units.count > 0 and _needs_benchmark(name):
                _benchmark(name)
            row["counts"] = cli.gate(name, stage, layout, args, notifier, record=row,
                                     units=units)
            row["status"] = "done"
        except BaseException:
            row["status"] = "FAILED"
            rows.extend({"stage": s, "status": "not run"} for s in stages[i + 1:])
            _finish(layout, rows, notifier, f"run FAILED at {name}", priority=1)
            raise
    flags = _flags(rows)
    if flags:
        _finish(layout, rows, notifier, f"run done with failures ({len(flags)} flagged)",
                priority=1)
        return EXIT_FLAGGED
    _finish(layout, rows, notifier, "run done", priority=0)
    return 0
