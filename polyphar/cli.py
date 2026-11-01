"""polyphar command line.

Each subcommand is bound to a stage module in STAGES. A stage module exposes
    units(layout, args) -> WorkUnits
    run(layout, args, progress) -> dict
and may define add_arguments(parser) to add stage-specific CLI options.
To wire a new stage, add or change one line in STAGES; subcommands whose module
does not exist yet print "not implemented" and exit 2.
"""

from __future__ import annotations

import argparse
import importlib
import os
import signal
import sys
import time
from datetime import datetime
from pathlib import Path

from polyphar import __version__
from polyphar.common import RunLayout
from polyphar.eta import (PARALLEL_EFFICIENCY, Calibration, Progress, estimate_stage,
                          format_eta, polyphar_home, slots_for)
from polyphar.notify import notify

# subcommand -> (stage module, help text)
STAGES = {
    "ids": ("polyphar.stages.ids", "build 00_ids/ids.txt from an xlsx or txt source"),
    "uniprot": ("polyphar.stages.uniprot", "fetch UniProt records"),
    "structures": ("polyphar.stages.structures", "download PDB/AlphaFold structures"),
    "fix": ("polyphar.stages.fix", "repair structures with pdbfixer"),
    "pockets": ("polyphar.stages.pockets", "detect pockets (DoGSiteScorer, host only)"),
    "similarity": ("polyphar.stages.similarity", "PocketMatch pocket similarity"),
    "select": ("polyphar.stages.select", "select polypharmacology targets"),
    "dock": ("polyphar.stages.dock", "dock the ligand library with smina"),
    "run": ("polyphar.stages.run", "run all stages in order"),
    "estimate": ("polyphar.stages.estimate", "print ETA for the next and later stages"),
    "benchmark": ("polyphar.stages.benchmark", "per-machine calibration"),
}


def default_run_dir() -> str:
    return "/work" if Path("/work").is_dir() else os.getcwd()


def default_nprocs() -> str:
    """Raw default for --nprocs: $NPROCS if set and non-empty, else "4".

    Returned as a string so argparse runs it through _nprocs_type at parse time:
    a bad $NPROCS then gives a clear usage error, and --help still works."""
    return os.environ.get("NPROCS") or "4"


def _nprocs_type(value: str) -> int:
    try:
        n = int(value)
    except (TypeError, ValueError):
        n = 0
    if n < 1:
        raise argparse.ArgumentTypeError(
            f"must be an integer >= 1 (got {value!r}); check --nprocs and $NPROCS")
    return n


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="polyphar", description="Polypharmacology workflow.")
    parser.add_argument("--version", action="version", version=f"polyphar {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>", required=True)
    for name, (_module, help_text) in STAGES.items():
        p = sub.add_parser(name, help=help_text, description=help_text)
        p.add_argument("--run-dir", default=default_run_dir(),
                       help="run folder (default: /work if it exists, else cwd)")
        p.add_argument("--nprocs", type=_nprocs_type, default=default_nprocs(),
                       help="worker processes (default: $NPROCS or 4)")
        p.add_argument("--dry-run", action="store_true",
                       help="print units and ETA, then exit")
        p.add_argument("--no-notify", action="store_true",
                       help="do not send notifications")
        stage = _load_stage(name)
        if stage is not None and hasattr(stage, "add_arguments"):
            stage.add_arguments(p)   # optional hook: stage-specific options
    return parser


def _load_stage(name: str):
    module_name = STAGES[name][0]
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        # Only "this stage does not exist yet"; a missing dependency inside an
        # existing stage must still raise.
        if exc.name and (module_name == exc.name or module_name.startswith(exc.name + ".")):
            return None
        raise


# Not per-stage gated: estimate and benchmark run directly; run (stages/run.py)
# passes every stage through gate() itself.
UNGATED = frozenset({"run", "estimate", "benchmark"})


def _no_notify(title, message, priority=0):
    return False


STOP_SIGNALS = (signal.SIGTERM, signal.SIGHUP)
FLAG_KEYS = ("failed", "incomplete_targets")
EXIT_FLAGGED = 3


def _install_sigterm_handler():
    """Turn SIGTERM (docker stop, kill) into SystemExit(143) and SIGHUP (ssh drop, closed
    terminal) into SystemExit(129), so an interrupted run unwinds through the failure
    notification, the run summary and the pool cleanup. Returns {signal: previous handler}
    for _restore_handlers, or None when they cannot be installed (not the main thread)."""
    def handler(signum, frame):
        raise SystemExit(128 + signum)
    previous = {}
    try:
        for sig in STOP_SIGNALS:
            previous[sig] = signal.signal(sig, handler)
    except ValueError:
        _restore_handlers(previous)
        return None
    return previous


def _restore_handlers(previous) -> None:
    for sig, handler in (previous or {}).items():
        signal.signal(sig, handler)


def flagged(counts) -> list[str]:
    """`key=value` for every failure flag (failed units, incomplete targets) set in counts."""
    out = []
    for key in FLAG_KEYS:
        v = (counts or {}).get(key)
        if isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0:
            out.append(f"{key}={v:g}")
    return out


def make_notifier(layout, args):
    """notify() bound to the run folder, or a no-op with --no-notify."""
    if args.no_notify:
        return _no_notify

    def notifier(title, message, priority=0):
        return notify(title, message, priority, run_dir=layout.run_dir)
    return notifier


def gate(name, stage, layout, args, notifier, record=None, units=None) -> dict:
    """ETA gate: units -> estimate -> print + notify start -> (dry-run: stop) ->
    run with Progress -> done/fail. Returns the stage's counts ({} for a dry run).
    `record`, when given, receives units / estimate_s / elapsed_s (also on failure).
    `units`, when given, is the stage's units() computed just before (by `run`)."""
    record = {} if record is None else record
    tpj = getattr(args, "threads_per_job", 1) or 1
    if units is None:
        units = stage.units(layout, args)
    calib = Calibration.load(polyphar_home())
    est = estimate_stage(name, units, args.nprocs, calib, PARALLEL_EFFICIENCY,
                         threads_per_job=tpj)
    record.update(units=f"{units.count} ({units.cost_basis})", estimate_s=est.seconds)
    eta = format_eta(est, datetime.now())
    line = f"{name}: {units.count} units ({units.cost_basis}) - ETA {eta}"
    print(line, flush=True)
    if args.dry_run:
        return {}   # print the ETA only; nothing starts, nothing is notified
    notifier(f"{name} starting", line)
    total_size = float(sum(units.sizes)) if units.sizes else None
    progress = Progress(name, units.count, est, notifier, total_size=total_size,
                        calib=calib, slots=slots_for(name, args.nprocs, tpj) or 1.0,
                        log_path=layout.logs / f"{name}.log")
    previous = _install_sigterm_handler()
    t0 = time.monotonic()
    try:
        counts = stage.run(layout, args, progress)
    except BaseException as exc:   # incl. KeyboardInterrupt / SIGTERM (docker stop)
        record["elapsed_s"] = time.monotonic() - t0
        progress.fail(exc)
        raise
    finally:
        _restore_handlers(previous)
    record["elapsed_s"] = time.monotonic() - t0
    progress.done(counts)
    return counts or {}


def _run_gated(name, stage, layout, args) -> int:
    """0, or 3 (like `run`) when the stage finished with failed units / incomplete
    targets (its done notification then goes out at priority 1)."""
    counts = gate(name, stage, layout, args, make_notifier(layout, args))
    if flagged(counts):
        print(f"{name}: finished with failures: {', '.join(flagged(counts))} "
              f"(see {name} failed.tsv); exit {EXIT_FLAGGED}", file=sys.stderr)
        return EXIT_FLAGGED
    return 0


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    stage = _load_stage(args.command)
    if stage is None:
        print(f"polyphar {args.command}: not implemented", file=sys.stderr)
        return 2
    layout = RunLayout(args.run_dir)
    if args.command not in UNGATED:
        return _run_gated(args.command, stage, layout, args)
    if hasattr(stage, "main"):   # run: its own dry-run and multi-stage gate
        return stage.main(layout, args)
    if args.dry_run:
        wu = stage.units(layout, args)
        print(f"{args.command}: {wu.count} units ({wu.cost_basis})")
        return 0
    stage.run(layout, args, print)
    return 0


if __name__ == "__main__":
    sys.exit(main())
