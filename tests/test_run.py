"""`polyphar run`: every stage in order through the same ETA gate (fake stages)."""
import argparse
import json
import os
import signal
import sys
import types

import pytest

from polyphar import cli, common, notify as notify_mod
from polyphar.common import WorkUnits
from polyphar.eta import Calibration
from polyphar.stages import run as run_mod

# imported, not duplicated: a local copy silently stops stubbing a stage the moment a
# new one is added to the pipeline
from polyphar.benchmark import ORDER, STAGE_DIRS as DIRS  # noqa: E402
# the single source of truth: zipping ORDER against a hardcoded tuple silently
# shifts every mapping when a stage is added


@pytest.fixture
def fake(monkeypatch, tmp_path):
    """Replace every stage with a fake that keeps the real stage's options. events is
    one ordered log of units / run / notifications / benchmarks / prints."""
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    events, notes = [], []
    # remaining: units left per stage (default 2); a finished stage writes its manifest
    # and has 0 left; a failed one keeps 1 left (part of the work done)
    state = {"fail": None, "counts": {}, "remaining": {}}

    for name in ORDER:
        real = cli._load_stage(name)
        mod = types.ModuleType(f"fake_run_{name}")

        def units(layout, args, name=name):
            events.append(("units", name))
            return WorkUnits(state["remaining"].get(name, 2), "things")

        def run(layout, args, progress, name=name):
            events.append(("run", name))
            for _ in range(state["remaining"].get(name, 2)):
                progress.advance()
            if state["fail"] == name:
                state["remaining"][name] = 1
                raise RuntimeError(f"{name} blew up")
            counts = {"failed": 0, **state["counts"].get(name, {})}
            d = layout.run_dir / DIRS[name]
            d.mkdir(parents=True, exist_ok=True)
            manifest = {"counts": counts, "git_commit": common._git_commit()}
            manifest.update(state.get("manifest_extra", {}).get(name, {}))
            (d / "manifest.json").write_text(json.dumps(manifest))
            if name == "ids":
                layout.ids_txt.write_text("P25815\n")
            state["remaining"][name] = 0
            return counts

        mod.units, mod.run = units, run
        if hasattr(real, "add_arguments"):
            mod.add_arguments = real.add_arguments
        monkeypatch.setitem(sys.modules, mod.__name__, mod)
        monkeypatch.setitem(cli.STAGES, name, (mod.__name__, "fake"))

    def fake_notify(title, message, priority=0, **kw):
        notes.append((title, message, priority))
        events.append(("notify", title))
        return True

    monkeypatch.setattr(cli, "notify", fake_notify)

    def fake_bench(stages, calib, log=print):
        events.append(("benchmark", tuple(stages)))
        if state.get("bench_sigterm"):
            os.kill(os.getpid(), state.get("bench_signal", signal.SIGTERM))
        for s in stages:
            calib.update(s, 1.0)
        calib.save()
        return {s: {"status": "calibrated"} for s in stages}

    monkeypatch.setattr(run_mod, "run_benchmark", fake_bench)
    real_preflight = run_mod.preflight_problems
    monkeypatch.setattr(run_mod, "preflight_problems",
                        lambda layout, args, stages: state.get("preflight", []))
    state["real_preflight"] = real_preflight
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    return types.SimpleNamespace(events=events, notes=notes, state=state, run=run_dir,
                                 home=tmp_path / "home")


def _calibrate_all(home):
    c = Calibration.load(home)
    for s in ("uniprot", "structures", "fix", "pockets", "annotate",
              "similarity", "dock"):
        c.update(s, 1.0)
    c.save()


def _main(fake, *extra):
    return cli.main(["run", "--run-dir", str(fake.run), "--nprocs", "2", *extra])


def _runs(events):
    return [e[1] for e in events if e[0] == "run"]


def test_runs_all_stages_in_order_through_the_gate(fake, capsys):
    _calibrate_all(fake.home)
    (fake.run / "inputs").mkdir()
    (fake.run / "inputs" / "ids.txt").write_text("P25815\n")
    assert _main(fake, "--txt", "inputs/ids.txt") == 0
    # --to defaults to `dock`, not ORDER[-1]: `stats` needs a guideline positive set
    assert _runs(fake.events) == list(ORDER[:-1])
    # per stage: units -> "<stage> starting" notification -> run (the single-stage gate)
    for name in ORDER[:-1]:
        i_units = fake.events.index(("units", name))
        i_start = fake.events.index(("notify", f"{name} starting"))
        i_run = fake.events.index(("run", name))
        assert i_units < i_start < i_run
    start_msgs = {t: m for t, m, _ in fake.notes if t.endswith(" starting")}
    assert "2 units (things) - ETA" in start_msgs["fix starting"]
    out = capsys.readouterr().out
    assert "fix: 2 units (things) - ETA" in out
    # final summary: printed and notified, one row per stage
    title, msg, prio = fake.notes[-1]
    assert title.startswith("run done") and prio == 0
    for name in ORDER[:-1]:
        assert name in msg
    assert "run summary" in out
    assert (fake.run / "logs" / "run.log").read_text().count("run summary") == 1


def test_from_to_select_a_contiguous_range(fake):
    _calibrate_all(fake.home)
    assert _main(fake, "--from", "fix", "--to", "select") == 0
    assert _runs(fake.events) == ["fix", "pockets", "annotate", "select"]


def test_from_after_to_is_a_usage_error(fake, capsys):
    assert _main(fake, "--from", "dock", "--to", "fix") == 2
    assert "--from" in capsys.readouterr().err
    assert _runs(fake.events) == []


def test_ids_skipped_as_done_when_ids_txt_exists_and_no_ids_options(fake, capsys):
    _calibrate_all(fake.home)
    (fake.run / "00_ids").mkdir()
    (fake.run / "00_ids" / "ids.txt").write_text("P25815\n")
    assert _main(fake, "--to", "uniprot") == 0
    assert _runs(fake.events) == ["uniprot"]
    out = capsys.readouterr().out
    assert "skipped (done)" in out
    assert "skipped (done)" in fake.notes[-1][1]


def test_ids_options_required_without_ids_txt(fake, capsys):
    assert _main(fake, "--to", "uniprot") == 2
    assert "--txt or --xlsx" in capsys.readouterr().err
    assert _runs(fake.events) == []


def test_first_failed_stage_stops_the_run(fake, capsys):
    _calibrate_all(fake.home)
    fake.state["fail"] = "fix"
    with pytest.raises(RuntimeError, match="fix blew up"):
        _main(fake, "--from", "uniprot")
    assert _runs(fake.events) == ["uniprot", "structures", "fix"]
    titles = [(t, p) for t, _, p in fake.notes]
    assert ("fix FAILED", 1) in titles                 # from the stage gate
    title, msg, prio = fake.notes[-1]
    assert prio == 1 and "fix" in title and "FAILED" in title
    assert "not run" in msg                             # later stages listed as not run
    assert msg.splitlines()[0] == "failures: fix FAILED"   # never "failures: none"
    assert "FAILED" in capsys.readouterr().out


def test_failed_units_and_incomplete_targets_are_flagged_not_fatal(fake, capsys):
    _calibrate_all(fake.home)
    fake.state["counts"] = {"structures": {"failed": 3},
                            "dock": {"failed": 0, "incomplete_targets": 2}}
    assert _main(fake, "--from", "uniprot") == 3
    assert _runs(fake.events)[-1] == "dock"
    title, msg, prio = fake.notes[-1]
    assert "failures" in title.lower() and prio == 1
    # the failures line comes first: notify truncates long messages
    assert msg.splitlines()[0] == "failures: structures failed=3; dock incomplete_targets=2"
    out = capsys.readouterr().out
    assert "structures failed=3" in out and "dock incomplete_targets=2" in out


def test_clean_run_says_no_failures(fake):
    _calibrate_all(fake.home)
    assert _main(fake, "--from", "uniprot") == 0
    title, msg, _ = fake.notes[-1]
    assert title == "run done" and "failures: none" in msg


def test_uncalibrated_stage_is_benchmarked_first_except_pockets(fake, capsys):
    c = Calibration.load(fake.home)
    c.update("uniprot", 1.0)      # calibrated: no benchmark
    c.save()
    assert _main(fake, "--from", "uniprot", "--to", "similarity") == 0
    benches = [e[1] for e in fake.events if e[0] == "benchmark"]
    assert benches == [("structures",), ("fix",), ("similarity",)]
    # each benchmark runs after its stage's units() (no benchmark when 0 units are left)
    # and before the stage's ETA / start notification
    for s in ("structures", "fix", "similarity"):
        i = fake.events.index(("benchmark", (s,)))
        assert fake.events.index(("units", s)) < i < fake.events.index(("notify",
                                                                        f"{s} starting"))
    assert "pockets: 2 units (things) - ETA uncalibrated" in capsys.readouterr().out


def test_dry_run_prints_the_estimate_table_and_runs_nothing(fake, capsys, monkeypatch):
    monkeypatch.setattr(run_mod, "estimate_rows",
                        lambda layout, nprocs, calib: ([{"stage": "ids", "kind": "trivial",
                                                         "units": "-", "seconds": 0.0,
                                                         "note": ""}], []))
    assert _main(fake, "--dry-run") == 0
    out = capsys.readouterr().out
    assert out.startswith("stage ")
    assert _runs(fake.events) == [] and fake.notes == []
    assert not (fake.run / "logs").exists()


def test_run_accepts_the_stage_options():
    args = cli.build_parser().parse_args(
        ["run", "--txt", "inputs/ids.txt", "--ph", "7.0", "--worker-timeout", "60",
         "--pm-bin-dir", "/opt/pm", "--per-protein", "3", "--exclude",
         "--library", "inputs/x.sdf", "--seed", "7", "--replicates", "5",
         "--threads-per-job", "2", "--unit-timeout", "100"])
    assert args.txt == "inputs/ids.txt" and args.per_protein == 3 and args.exclude == []
    assert args.library == "inputs/x.sdf" and args.replicates == 5 and args.seed == 7
    assert args.from_stage == "ids" and args.to_stage == "dock"
    # ids' source options stay mutually exclusive
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["run", "--txt", "a", "--xlsx", "b"])


def test_conflicting_option_definitions_are_rejected():
    def a(p):
        p.add_argument("--seed", type=int, default=1)

    def b(p):
        p.add_argument("--seed", type=int, default=2)

    parser = argparse.ArgumentParser()
    with pytest.raises(ValueError, match="--seed"):
        run_mod.merge_stage_options(parser, [("one", a), ("two", b)])


def test_identical_option_definitions_are_merged():
    def a(p):
        p.add_argument("--seed", type=int, default=1)

    parser = argparse.ArgumentParser()
    run_mod.merge_stage_options(parser, [("one", a), ("two", a)])
    assert parser.parse_args(["--seed", "3"]).seed == 3


def test_notify_records_every_notification_in_the_run_log(monkeypatch, tmp_path):
    monkeypatch.delenv("PUSHOVER_APP_TOKEN", raising=False)
    monkeypatch.delenv("PUSHOVER_USER_KEY", raising=False)
    notify_mod.notify("fix starting", "fix: 2 units\nsecond line", priority=1, run_dir=tmp_path)
    lines = (tmp_path / "logs" / "notifications.log").read_text().splitlines()
    assert len(lines) == 1
    stamp, prio, sent, title, msg = lines[0].split("\t")
    assert prio == "1" and sent == "not sent" and title == "fix starting"
    assert msg == "fix: 2 units\\nsecond line"
    assert "T" in stamp


def test_resume_skips_finished_stages_and_estimates_only_remaining_units(fake, capsys):
    _calibrate_all(fake.home)
    (fake.run / "inputs").mkdir()
    (fake.run / "inputs" / "ids.txt").write_text("P25815\n")
    fake.state["fail"] = "fix"
    with pytest.raises(RuntimeError):
        _main(fake, "--txt", "inputs/ids.txt")
    capsys.readouterr()
    fake.events.clear()
    fake.notes.clear()
    fake.state["fail"] = None

    assert _main(fake) == 0     # no ids options: 00_ids/ids.txt exists
    assert _runs(fake.events) == ["fix", "pockets", "annotate", "select",
                                  "similarity", "ligands", "dock"]
    out = capsys.readouterr().out
    # ETA only for the remaining work: fix has 1 unit left
    assert "fix: 1 units (things) - ETA" in out
    starts = [t for t, _, _ in fake.notes if t.endswith(" starting")]
    assert starts == [f"{s} starting" for s in ("fix", "pockets", "annotate", "select",
                                                "similarity", "ligands", "dock")]
    summary = fake.notes[-1][1]
    for s in ("ids", "uniprot", "structures"):
        line = next(l for l in summary.splitlines() if l.startswith(s + " "))
        assert "skipped (done)" in line
    assert "fix" in summary and "failures: none" in summary


def test_stage_without_manifest_and_zero_units_still_runs(fake, capsys):
    _calibrate_all(fake.home)
    fake.state["remaining"]["uniprot"] = 0      # no manifest yet: bookkeeping still to do
    assert _main(fake, "--from", "uniprot", "--to", "uniprot") == 0
    assert _runs(fake.events) == ["uniprot"]


def test_skipped_stage_with_flagged_manifest_counts_is_still_flagged(fake):
    _calibrate_all(fake.home)
    d = fake.run / "07_docking"
    d.mkdir()
    (d / "manifest.json").write_text(json.dumps({"counts": {"failed": 0,
                                                            "incomplete_targets": 1}}))
    fake.state["remaining"]["dock"] = 0
    assert _main(fake, "--from", "dock", "--to", "dock") == 3
    assert _runs(fake.events) == []
    title, msg, prio = fake.notes[-1]
    assert prio == 1 and msg.splitlines()[0] == "failures: dock incomplete_targets=1"


def test_no_benchmark_for_a_stage_with_no_units_left(fake):
    fake.state["remaining"]["dock"] = 0
    d = fake.run / "07_docking"
    d.mkdir()
    (d / "manifest.json").write_text(json.dumps({"counts": {"failed": 0}}))
    assert _main(fake, "--from", "dock") == 0
    assert [e for e in fake.events if e[0] == "benchmark"] == []


def test_sigterm_during_auto_benchmark_ends_in_run_failed_summary(fake):
    fake.state["bench_sigterm"] = True
    previous = signal.getsignal(signal.SIGTERM)
    with pytest.raises(SystemExit) as exc:
        _main(fake, "--from", "uniprot", "--to", "fix")
    assert exc.value.code == 143
    assert _runs(fake.events) == []
    title, msg, prio = fake.notes[-1]
    assert title == "run FAILED at uniprot" and prio == 1
    assert signal.getsignal(signal.SIGTERM) is previous      # handler restored


def test_options_of_stages_outside_the_range_are_warned(fake, capsys):
    _calibrate_all(fake.home)
    assert _main(fake, "--from", "uniprot", "--to", "uniprot", "--per-protein", "3",
                 "--library", "inputs/x.sdf", "--threads-per-job", "1") == 0
    err = capsys.readouterr().err.splitlines()
    warned = [l for l in err if "warning" in l]
    assert any("--per-protein" in l and "select" in l for l in warned)
    assert any("--library" in l and "dock" in l for l in warned)
    assert not any("--threads-per-job" in l for l in warned)   # default value: not "given"
    assert len(warned) == 2


def test_sighup_ends_in_run_failed_summary_exit_129(fake):
    fake.state["bench_sigterm"] = True
    fake.state["bench_signal"] = signal.SIGHUP
    previous = signal.getsignal(signal.SIGHUP)
    with pytest.raises(SystemExit) as exc:
        _main(fake, "--from", "uniprot", "--to", "fix")
    assert exc.value.code == 129
    title, msg, prio = fake.notes[-1]
    assert title == "run FAILED at uniprot" and prio == 1
    assert signal.getsignal(signal.SIGHUP) is previous


# --- item 1: a skipped stage must have unchanged inputs ---------------------------------

def _manifest(fake, stage, **fields):
    d = fake.run / DIRS[stage]
    d.mkdir(parents=True, exist_ok=True)
    m = {"counts": {"failed": 0}, "git_commit": common._git_commit(), **fields}
    (d / "manifest.json").write_text(json.dumps(m))


def test_zero_unit_stage_with_changed_input_runs_through_the_gate(fake, capsys):
    _calibrate_all(fake.home)
    ids = fake.run / "00_ids" / "ids.txt"
    ids.parent.mkdir()
    ids.write_text("P25815\nP53350\n")
    _manifest(fake, "uniprot", inputs=[{"path": str(ids), "sha256": common.sha256_file(ids)}])
    fake.state["remaining"]["uniprot"] = 0
    assert _main(fake, "--from", "uniprot", "--to", "uniprot") == 0
    assert _runs(fake.events) == []                  # unchanged: skipped
    ids.write_text("P25815\n")                      # a subset
    capsys.readouterr()
    assert _main(fake, "--from", "uniprot", "--to", "uniprot") == 0
    assert _runs(fake.events) == ["uniprot"]         # changed: run, not skipped
    assert "inputs changed since its manifest" in capsys.readouterr().out


def test_manifest_input_recorded_under_another_run_dir_path_is_found(fake):
    """Paths recorded in the container (/work/...) resolve inside the run folder."""
    _calibrate_all(fake.home)
    ids = fake.run / "00_ids" / "ids.txt"
    ids.parent.mkdir()
    ids.write_text("P25815\n")
    _manifest(fake, "uniprot", inputs=[{"path": "/work/00_ids/ids.txt",
                                        "sha256": common.sha256_file(ids)}])
    fake.state["remaining"]["uniprot"] = 0
    assert _main(fake, "--from", "uniprot", "--to", "uniprot") == 0
    assert _runs(fake.events) == []


def test_zero_unit_stage_with_missing_input_runs(fake):
    _calibrate_all(fake.home)
    _manifest(fake, "fix", inputs=[{"path": str(fake.run / "02_structures" / "gone.csv"),
                                    "sha256": "0" * 64}])
    fake.state["remaining"]["fix"] = 0
    assert _main(fake, "--from", "fix", "--to", "fix") == 0
    assert _runs(fake.events) == ["fix"]


# --- item 10 / 12: provenance warnings -------------------------------------------------

@pytest.mark.parametrize("commit", ["0123abcd" * 5, "unknown", None])
def test_skipping_a_stage_from_another_code_version_warns(fake, capsys, commit):
    _calibrate_all(fake.home)
    fields = {} if commit is None else {"git_commit": commit}
    d = fake.run / DIRS["dock"]
    d.mkdir()
    (d / "manifest.json").write_text(json.dumps({"counts": {"failed": 0}, **fields}))
    fake.state["remaining"]["dock"] = 0
    assert _main(fake, "--from", "dock") == 0     # a warning only, not a failure
    out = capsys.readouterr().out
    assert "dock: warning: skipped stage was produced by code version" in out
    assert any(t == "dock: skipped stage from another code version" for t, _, _ in fake.notes)


def test_skipping_a_stage_of_the_current_code_version_does_not_warn(fake, capsys):
    _calibrate_all(fake.home)
    _manifest(fake, "dock")
    fake.state["remaining"]["dock"] = 0
    assert _main(fake, "--from", "dock") == 0
    assert "code version" not in capsys.readouterr().out


def test_unknown_code_version_warns_once_at_run_start(fake, capsys, monkeypatch):
    _calibrate_all(fake.home)
    monkeypatch.setattr(common, "_git_commit", lambda: "unknown")
    assert _main(fake, "--from", "uniprot", "--to", "fix") == 0
    out = capsys.readouterr().out
    assert out.count("run: warning: code version unknown") == 1
    titles = [t for t, _, _ in fake.notes]
    assert titles.count("run: code version unknown") == 1
    assert titles.index("run: code version unknown") < titles.index("uniprot starting")


# --- item 8: preflight --------------------------------------------------------------------

def test_preflight_problem_stops_before_the_first_stage(fake, capsys):
    _calibrate_all(fake.home)
    fake.state["preflight"] = ["dock: ligand library /x.sdf not found"]
    assert _main(fake, "--from", "uniprot") == 2
    assert _runs(fake.events) == [] and not [e for e in fake.events if e[0] == "units"]
    assert "ligand library /x.sdf not found" in capsys.readouterr().err


def test_preflight_checks_library_smina_and_pocketmatch(fake, tmp_path, monkeypatch):
    from polyphar.common import RunLayout
    from polyphar.stages import dock

    check = fake.state["real_preflight"]
    layout = RunLayout(fake.run)
    args = argparse.Namespace(library="inputs/missing.sdf", pm_bin_dir=str(tmp_path / "nopm"))
    monkeypatch.setattr(dock, "smina_version", lambda: "unavailable")
    problems = check(layout, args, ["similarity", "select", "dock"])
    text = "\n".join(problems)
    assert "missing.sdf not found" in text
    assert "smina" in text and "not runnable" in text
    assert "PocketMatch" in text
    assert check(layout, args, ["uniprot", "structures"]) == []   # neither in range


# --- item 9: a stage with too many failed units stops the run ---------------------------

def test_too_many_failed_units_stop_the_run(fake, capsys):
    _calibrate_all(fake.home)
    fake.state["counts"]["pockets"] = {"failed": 3, "models": 10}
    assert _main(fake, "--from", "pockets") == 1
    assert _runs(fake.events) == ["pockets"]
    title, msg, prio = fake.notes[-1]
    assert title.startswith("run STOPPED at pockets") and prio == 1
    assert "3/10" in title + msg


def test_keep_going_overrides_the_failed_fraction(fake):
    _calibrate_all(fake.home)
    fake.state["counts"]["pockets"] = {"failed": 3, "models": 10}
    assert _main(fake, "--from", "pockets", "--keep-going") == 3
    assert _runs(fake.events) == ["pockets", "annotate", "select", "similarity",
                                  "ligands", "dock"]


def test_failed_fraction_at_the_limit_continues(fake):
    _calibrate_all(fake.home)
    fake.state["counts"]["fix"] = {"failed": 2, "inputs": 10}
    assert _main(fake, "--from", "fix", "--to", "pockets") == 3
    assert _runs(fake.events) == ["fix", "pockets"]


def test_failed_fraction_excludes_expected_exclusions():
    f = run_mod.failed_fraction
    # unmapped pairs are expected exclusions: not in failed, not in the denominator
    assert f("structures", {"failed": 1, "pdb_entries_listed": 12, "accessions": 3,
                            "unmapped": {"no SIFTS chain for accession": 5}}) == (1, 10)
    assert f("pockets", {"failed": 1, "models": 4, "outside_segment": 50}) == (1, 4)
    assert f("fix", {"failed": 0, "inputs": 0}) is None
    assert f("dock", {"failed": 3}) is None
