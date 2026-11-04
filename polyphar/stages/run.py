"""`polyphar run`: the stages ids -> uniprot -> structures -> fix -> pockets -> annotate ->
select -> similarity -> ligands -> dock in order, and `stats` after them on request (or the
contiguous range --from .. --to; benchmark.ORDER is the authority and ends at stats, but
--to defaults to **dock**, because `stats` needs a clinical-guideline positive set that
will be absent for most of this project's life -- `--to stats` includes it, and a run that
ends at dock without one says, in one line, what would unblock it). Each stage runs through
the same ETA gate as its single-stage
subcommand (cli.gate): exact units, ETA printed and notified before the first unit,
milestone progress, done / FAILED.

* `run` accepts the union of the stages' own options (built from each stage's
  add_arguments(); two stages defining one option differently is an error).
* ids: runs when --txt / --xlsx is given; without them an existing 00_ids/ids.txt counts
  as done (skipped), and a missing one is an error before anything starts.
* Before the first stage, a preflight checks what later stages need: with ligands or dock
  in range the ligand library exists (and, for ligands with work left, RDKit imports);
  with dock in range `smina --version` works; with similarity in range the PocketMatch
  binaries are there; with stats in range the positive set is there. A problem exits 2
  before anything starts.
* Every stage's units() is computed first. A stage with 0 units and a manifest.json whose
  recorded inputs all still have their recorded sha256 is "skipped (done)": not
  benchmarked, not run, no notifications (its manifest counts still feed the failure
  flags). With a changed or missing input (e.g. a shrunken ids.txt after the stage
  finished), or without a manifest, it runs through the gate and rebuilds its tables from
  the current inputs. Skipping a stage whose manifest comes from another code version
  (git_commit), or an unknown one, prints and notifies a warning (it does not block).
* A code version that is "unknown" (image built without POLYPHAR_GIT_COMMIT) is warned
  about once at the start (Methods need it).
* Options of stages outside --from..--to are warned about (stderr) and ignored.
* A stage with remaining units and no calibration on this machine is benchmarked first (`benchmark --stages
  <stage>`), except pockets: DoGSite3 / DoGSiteScorer is host-only (bin/polyphar-docker
  calibrates it, per pocket tool, before `run`); uncalibrated, its ETA comes with the 5 %
  milestone.
* A finished pockets stage with nothing left to run is rerun when its pocket selection
  (--min-volume / --drugscore-cutoff) changed (pockets.needs_rerun).
* The first failed stage stops the run (priority-1 notification). So does structures, fix
  or pockets when failed / work items exceeds --max-failed-fraction (default 0.2): exit 1,
  "run STOPPED"; --keep-going overrides. Expected exclusions (structures unmapped pairs,
  pockets outside_segment) are not failures and not in the fraction. SIGTERM (docker stop) is
  turned into a failure for the whole run, so an interrupted benchmark or units() also ends
  in the "run FAILED" summary. A stage that finishes with counts["failed"] > 0, or dock with
  counts["incomplete_targets"] > 0, does not stop the run but is flagged in the final
  summary (printed, notified at priority 1, appended to logs/run.log).
* Exit codes: 0 success; 1 a stage failed (the exception propagates; 143 for SIGTERM,
  129 for SIGHUP) or failed too often;
  2 usage error; 3 the run finished but failures were flagged; 4 the run stopped cleanly
  with work still outstanding (`--engine unidock` on a host with no GPU: everything up to
  dock is done, and dock-export / dock-worker / dock-merge finish it). 4 is not 0 because
  the documented workflow is `tmux` + `bin/polyphar-docker run`, where the guidance on
  stdout is what nobody reads and `polyphar run && <next>` is what proceeds.
* --dry-run prints the `polyphar estimate` table and does nothing else.

Heavy imports stay inside functions (cli.build_parser imports this module at parse time).
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime

from polyphar import common
from polyphar.benchmark import (BENCH_STAGES, ORDER, STAGE_DIRS, estimate_rows, format_table,
                                run_benchmark)
from polyphar.eta import Calibration, format_duration, polyphar_home

IDS_SOURCE = ("xlsx", "txt", "sheet", "id_column", "reviewed_column", "gene_column")
FLAG_KEYS = ("failed", "incomplete_targets")
EXIT_USAGE = 2
EXIT_FLAGGED = 3
EXIT_STOPPED = 1
# The run stopped cleanly with work still to do: nothing failed, but the pipeline is NOT
# finished. Distinct from 0, because every machine-readable signal of the GPU-less stop
# otherwise said success -- and `polyphar run && <next>` in the documented tmux workflow is
# exactly where nobody reads the stdout that carried the guidance.
EXIT_INCOMPLETE = 4
SKIPPED = "skipped (done)"
DEFAULT_MAX_FAILED_FRACTION = 0.2
# `stats` is a real pipeline stage and stays in benchmark.ORDER, but it is not in the
# default range: see --to in add_arguments.
DEFAULT_TO_STAGE = "dock"
FRACTION_STAGES = ("structures", "fix", "pockets")


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
    # NOT ORDER[-1]. `stats` needs inputs/positive_sets.csv, which is clinical-guideline
    # reading and will be absent for most of this project's life; making the everyday
    # command refuse to do the eight stages it CAN do inverts the cost asymmetry. Same
    # shape as _gpu_stop: finish what is possible, stop cleanly, say what is missing.
    # `--to stats` and `--from stats` still work, and the preflight check stays fatal when
    # the user puts `stats` in range explicitly.
    parser.add_argument("--to", dest="to_stage", choices=ORDER, default=DEFAULT_TO_STAGE,
                        metavar="STAGE",
                        help=f"last stage (default {DEFAULT_TO_STAGE}; pass --to stats to "
                             f"include the analysis)")
    parser.add_argument("--max-failed-fraction", type=_fraction,
                        default=DEFAULT_MAX_FAILED_FRACTION, metavar="F",
                        help="stop the run when structures, fix or pockets fails more than "
                             f"this fraction of its work items (default "
                             f"{DEFAULT_MAX_FAILED_FRACTION:g})")
    parser.add_argument("--keep-going", action="store_true",
                        help="never stop for --max-failed-fraction")
    hooks = []
    for name in ORDER:
        module = cli._load_stage(name)
        if module is not None and hasattr(module, "add_arguments"):
            hooks.append((name, module.add_arguments))
    # which stage each option belongs to, for the out-of-range warning
    parser.set_defaults(run_option_owners=merge_stage_options(parser, hooks))


def _fraction(value: str) -> float:
    try:
        f = float(value)
    except (TypeError, ValueError):
        f = -1.0
    if not 0.0 <= f <= 1.0:
        raise argparse.ArgumentTypeError(f"must be a number in [0, 1] (got {value!r})")
    return f


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


def summary_table(rows, run_name: str, note=None) -> str:
    """The run's table, closed by `failures: ...` -- or by `note`, when the run ended for a
    reason the failure count does not describe (the GPU-less stop: nothing failed, but the
    pipeline is not finished, and "failures: none" alone reads as success)."""
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
    lines.append(note or failures_line(rows))
    return "\n".join(lines)


def _finish(layout, rows, notifier, title, priority, note=None) -> str:
    from polyphar.notify import run_name

    text = summary_table(rows, run_name(layout.run_dir), note)
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
        notifier(title, (note or failures_line(rows)) + "\n" + text, priority=priority)
    except Exception:
        pass   # notifications never fail a run
    return text


def _notify_quietly(notifier, title, message, priority=0) -> None:
    try:
        notifier(title, message, priority=priority)
    except Exception:
        pass   # notifications never fail a run


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


def _manifest(layout, stage) -> dict | None:
    """The stage's manifest.json, or None when it has none ({} when unreadable)."""
    path = layout.run_dir / STAGE_DIRS[stage] / "manifest.json"
    if not path.is_file():
        return None
    try:
        m = json.loads(path.read_text(encoding="utf-8"))
        return m if isinstance(m, dict) else {}
    except (OSError, ValueError):
        return {}


def _manifest_counts(layout, stage) -> dict | None:
    """counts of the stage's manifest.json, or None when the stage has no manifest."""
    m = _manifest(layout, stage)
    if m is None:
        return None
    counts = m.get("counts")
    return counts if isinstance(counts, dict) else {}


def _input_path(layout, recorded: str):
    """The file a manifest input path names, in this run folder. Paths are recorded as
    the stage saw them (/work/... in the container, the host path outside it): when the
    recorded path does not exist, its longest tail that exists under the run folder is it."""
    from pathlib import Path

    p = Path(recorded)
    if p.is_file():
        return p
    parts = p.parts[1:] if p.is_absolute() else p.parts
    for i in range(len(parts)):
        candidate = layout.run_dir.joinpath(*parts[i:])
        if candidate.is_file():
            return candidate
    return None


def changed_inputs(layout, manifest: dict) -> list[str]:
    """The manifest's recorded inputs whose current file is missing or has another sha256
    ([] when every input is unchanged)."""
    from polyphar.common import sha256_file

    changed = []
    for entry in manifest.get("inputs") or []:
        try:
            recorded, sha = entry["path"], entry["sha256"]
        except (KeyError, TypeError):
            changed.append(str(entry))
            continue
        path = _input_path(layout, recorded)
        if path is None:
            changed.append(f"{recorded} (missing)")
        elif sha256_file(path) != sha:
            changed.append(f"{recorded} (changed)")
    return changed


def failed_fraction(stage: str, counts) -> tuple[int, int] | None:
    """(failed, work items) of structures / fix / pockets from their counts, or None when
    the stage is not checked or has no work items. Expected exclusions are no failures:
    structures' unmapped pairs (not in `failed`, removed from the denominator) and pockets'
    outside_segment pockets (counted apart from `failed`)."""
    counts = counts or {}

    def num(key):
        v = counts.get(key, 0)
        return v if isinstance(v, (int, float)) and not isinstance(v, bool) else 0

    if stage == "structures":
        unmapped = counts.get("unmapped") or {}
        n_unmapped = sum(v for v in unmapped.values() if isinstance(v, (int, float))) \
            if isinstance(unmapped, dict) else 0
        total = num("pdb_entries_listed") - n_unmapped + num("accessions")
    elif stage == "fix":
        total = num("inputs")
    elif stage == "pockets":
        total = num("models")
    else:
        return None
    if total <= 0:
        return None
    return int(num("failed")), int(total)


def preflight_problems(layout, args, stages) -> list[str]:
    """What the stages in range need from outside the run folder's earlier stages: the
    ligand library (ligands and dock), a working smina (dock), RDKit (ligands), the
    PocketMatch binaries (similarity), the positive set (stats)."""
    problems = []
    wants_library = {"dock", "ligands"} & set(stages)
    if wants_library:
        from polyphar.stages import dock, ligands

        who = "/".join(s for s in ("ligands", "dock") if s in wants_library)
        try:
            # the RAW library: whether a prepared one exists is dock's business, not the
            # preflight's, and `ligands` needs the source either way.
            source = dock._raw_library_path(layout, args, copy=False)
        except (OSError, ValueError) as exc:
            source = None
            problems.append(f"{who}: {exc}")
        if "ligands" in stages:
            # a missing RDKit is one missing dependency; caught per record it would make
            # the stage write an unprepared library and report success.
            rdkit = ligands.rdkit_problem()
            if rdkit and (source is None or ligands._finished(layout, source) is None):
                problems.append(f"ligands: {rdkit}")
    if {"dock", "control"} & set(stages):
        from polyphar.stages import dock
        from polyphar.stages import dock_unidock as unidock

        # engine-aware: demanding smina from a Uni-Dock run would stop a perfectly good
        # GPU host that has no CPU engine at all. `control` is in here because it runs the
        # same engine over its own ligand set (polyphar.stages.control).
        if dock.engine_name(args) == dock.ENGINE_UNIDOCK:
            # a GPU-less host is NOT a preflight failure: `run` does everything up to dock
            # and stops cleanly there (see _gpu_stop). What IS a failure is a GPU host
            # carrying a Uni-Dock too old to trust.
            if unidock.gpu_model() is not None:
                try:
                    unidock.check_version()
                except unidock.UniDockError as exc:
                    problems.append(f"dock: {exc}")
        elif dock.smina_version() == "unavailable":
            problems.append(f"dock: smina ({dock.smina_bin()}) not runnable "
                            "(`smina --version` failed)")
    if "stats" in stages:
        from polyphar.stages import stats

        # the positive set is guideline reading, not code: a run that would reach `stats`
        # without one stops here rather than after docking the whole library.
        path = stats._path(layout, getattr(args, "positive_set", stats.POSITIVE_SETS))
        if not path.is_file():
            # `stats` is only in range because the user asked for it (--to stats /
            # --from stats): --to defaults to dock. So this is fatal, as it should be.
            problems.append("stats: " + stats.POSITIVE_SET_MISSING.format(path=path)
                            + " (or drop --to stats: `run` stops after dock by default)")
    if "similarity" in stages:
        from polyphar.stages import similarity

        try:
            similarity._pm_bin_dir(args)
        except (OSError, ValueError) as exc:
            problems.append(f"similarity: {exc}")
    return problems


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

    problems = preflight_problems(layout, args, stages)
    if problems:
        return _usage("preflight failed, nothing was started:\n  " + "\n  ".join(problems))

    notifier = cli.make_notifier(layout, args)
    commit = common._git_commit()
    if commit == "unknown":
        msg = ("code version unknown: this polyphar was built without POLYPHAR_GIT_COMMIT "
               "(rebuild the image with --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD)); "
               "manifests will not say which code produced them")
        print(f"run: warning: {msg}", flush=True)
        _notify_quietly(notifier, "run: code version unknown", msg)
    args.run_git_commit = commit
    # docker stop -> SystemExit(143), ssh drop (SIGHUP) -> SystemExit(129) -> summary
    previous = cli._install_sigterm_handler()
    try:
        return _run_stages(cli, layout, args, stages, ids_given, notifier)
    finally:
        cli._restore_handlers(previous)


def _warn_code_version(name, manifest, current, notifier) -> None:
    recorded = manifest.get("git_commit") or "unknown"
    if recorded != "unknown" and recorded == current:
        return
    msg = (f"skipped stage was produced by code version {recorded}, this is {current}; "
           f"rerun it (delete {STAGE_DIRS[name]}/manifest.json) if the code change matters")
    print(f"{name}: warning: {msg}", flush=True)
    _notify_quietly(notifier, f"{name}: skipped stage from another code version", msg)


def _gpu_stop(args, name: str) -> str | None:
    """Why `run` must stop cleanly before this stage, or None.

    `--engine unidock` needs a CUDA GPU, and POSEIDON -- the host the run folders live on
    -- has none. Everything up to dock is still worth doing there, so the run completes it
    and stops at dock with the export/worker/merge instructions, rather than failing.

    `control` docks too (the co-crystal ligands), so it stops for the same reason."""
    if name not in ("dock", "control"):
        return None
    from polyphar.stages import dock
    from polyphar.stages import dock_unidock as unidock

    if dock.engine_name(args) != dock.ENGINE_UNIDOCK or unidock.gpu_model() is not None:
        return None
    if name == "control":
        # the export/worker/merge path in NO_GPU_MESSAGE is dock's bundle, and there is no
        # control-export: say what actually applies instead of pointing at the wrong tool.
        return (unidock.NO_GPU_MESSAGE.split("Run the dock stage")[0]
                + "`polyphar control` docks too (the co-crystal ligands) and has no "
                  "export/merge split: run it on the GPU host against the same run folder, "
                  "or use `--engine smina` for a CPU run.")
    return unidock.NO_GPU_MESSAGE


def _stats_not_run(layout, args, stages, notifier) -> None:
    """One line, printed and notified, when a run ends at `dock` with no positive set.

    `--to` defaults to dock because the positive set is clinical-guideline reading and will
    be absent for most of this project's life. That is the right default, but it must not
    be silent: a run that stops one stage short of the analysis has to say so, and say what
    would unblock it."""
    if stages[-1] != "dock" or "stats" in stages:
        return
    from polyphar.stages import stats

    path = stats._path(layout, getattr(args, "positive_set", stats.POSITIVE_SETS))
    if path.is_file():
        return
    msg = (f"stats not run: supply {path} (guideline reading) and run "
           f"polyphar run --from stats")
    print(msg, flush=True)
    _notify_quietly(notifier, "run: stats not run", msg)


def _run_stages(cli, layout, args, stages, ids_given, notifier) -> int:
    rows = []
    gpu_stop = None
    current = getattr(args, "run_git_commit", None) or common._git_commit()
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
            manifest = _manifest(layout, name) if units.count == 0 else None
            stale = (stage.needs_rerun(layout, args)
                     if manifest is not None and hasattr(stage, "needs_rerun") else None)
            if manifest is not None:
                changed = changed_inputs(layout, manifest)
                if stale:
                    print(f"{name}: 0 units left but {stale}: running it", flush=True)
                elif changed:
                    print(f"{name}: 0 units left but inputs changed since its manifest "
                          f"({'; '.join(changed[:5])}): running it", flush=True)
                else:
                    print(f"{name}: 0 units left and a manifest exists: {SKIPPED}", flush=True)
                    _warn_code_version(name, manifest, current, notifier)
                    counts = manifest.get("counts")
                    row.update(status=SKIPPED, counts=counts if isinstance(counts, dict) else {},
                               units=f"0 ({units.cost_basis})")
                    if _too_many_failures(layout, args, rows, i, stages, notifier):
                        return EXIT_STOPPED
                    continue
            # AFTER the "0 units left and a manifest exists" branch above, never before
            # it: a run folder whose docking is already complete -- the normal end state
            # of export / worker / merge -- must report `skipped (done)`, not `not run`
            # with instructions to go and export a bundle it has already finished.
            stop = _gpu_stop(args, name)
            if stop is not None:
                gpu_stop = (name, stop)
                rows.pop()
                rows.extend({"stage": s, "status": "not run"} for s in stages[i:])
                break
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
        if _too_many_failures(layout, args, rows, i, stages, notifier):
            return EXIT_STOPPED
    if gpu_stop is not None:
        name, stop = gpu_stop
        print(f"{name}: {stop}", flush=True)
        # The note leads the Pushover notification, so it has to be the advice that
        # applies to THIS stage: `control` has no export/worker/merge split (there is no
        # `polyphar control-export`), and telling its reader to run one contradicts the
        # line printed immediately above it.
        note = (f"stopped before {name}: run dock-export / dock-worker / dock-merge to "
                f"finish" if name == "dock" else
                f"stopped before {name}: rerun `polyphar {name}` on a GPU host against "
                f"this run folder, or use --engine smina")
        _finish(layout, rows, notifier,
                f"run stopped before {name}: this host has no GPU", priority=1, note=note)
        return EXIT_INCOMPLETE
    _stats_not_run(layout, args, stages, notifier)
    flags = _flags(rows)
    if flags:
        _finish(layout, rows, notifier, f"run done with failures ({len(flags)} flagged)",
                priority=1)
        return EXIT_FLAGGED
    _finish(layout, rows, notifier, "run done", priority=0)
    return 0


def _too_many_failures(layout, args, rows, i, stages, notifier) -> bool:
    """Stop the run (summary, priority 1) when the stage just done failed more than
    --max-failed-fraction of its work items, unless --keep-going."""
    row = rows[-1]
    ff = failed_fraction(row["stage"], row.get("counts"))
    limit = getattr(args, "max_failed_fraction", DEFAULT_MAX_FAILED_FRACTION)
    if ff is None or ff[0] <= limit * ff[1]:
        return False
    failed, total = ff
    what = (f"{failed}/{total} work items failed ({failed / total:.0%} > "
            f"--max-failed-fraction {limit:g})")
    if getattr(args, "keep_going", False):
        print(f"{row['stage']}: {what}; continuing (--keep-going)", flush=True)
        return False
    print(f"{row['stage']}: {what}; stopping the run (see {STAGE_DIRS[row['stage']]}/"
          "failed.tsv; --keep-going continues anyway)", flush=True)
    rows.extend({"stage": s, "status": "not run"} for s in stages[i + 1:])
    _finish(layout, rows, notifier, f"run STOPPED at {row['stage']}: {what}", priority=1)
    return True
