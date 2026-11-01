"""`polyphar run`: every stage in order through the same ETA gate (fake stages)."""
import argparse
import sys
import types

import pytest

from polyphar import cli, notify as notify_mod
from polyphar.common import WorkUnits
from polyphar.eta import Calibration
from polyphar.stages import run as run_mod

ORDER = ("ids", "uniprot", "structures", "fix", "pockets", "similarity", "select", "dock")


@pytest.fixture
def fake(monkeypatch, tmp_path):
    """Replace every stage with a fake that keeps the real stage's options. events is
    one ordered log of units / run / notifications / benchmarks / prints."""
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    events, notes = [], []
    state = {"fail": None, "counts": {}}

    for name in ORDER:
        real = cli._load_stage(name)
        mod = types.ModuleType(f"fake_run_{name}")

        def units(layout, args, name=name):
            events.append(("units", name))
            return WorkUnits(2, "things")

        def run(layout, args, progress, name=name):
            events.append(("run", name))
            progress.advance()
            progress.advance()
            if state["fail"] == name:
                raise RuntimeError(f"{name} blew up")
            return {"failed": 0, **state["counts"].get(name, {})}

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
        for s in stages:
            calib.update(s, 1.0)
        calib.save()
        return {s: {"status": "calibrated"} for s in stages}

    monkeypatch.setattr(run_mod, "run_benchmark", fake_bench)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    return types.SimpleNamespace(events=events, notes=notes, state=state, run=run_dir,
                                 home=tmp_path / "home")


def _calibrate_all(home):
    c = Calibration.load(home)
    for s in ("uniprot", "structures", "fix", "pockets", "similarity", "dock"):
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
    assert _runs(fake.events) == list(ORDER)
    # per stage: units -> "<stage> starting" notification -> run (the single-stage gate)
    for name in ORDER:
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
    for name in ORDER:
        assert name in msg
    assert "run summary" in out
    assert (fake.run / "logs" / "run.log").read_text().count("run summary") == 1


def test_from_to_select_a_contiguous_range(fake):
    _calibrate_all(fake.home)
    assert _main(fake, "--from", "fix", "--to", "select") == 0
    assert _runs(fake.events) == ["fix", "pockets", "similarity", "select"]


def test_from_after_to_is_an_error(fake):
    with pytest.raises(SystemExit, match="--from"):
        _main(fake, "--from", "dock", "--to", "fix")
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


def test_ids_options_required_without_ids_txt(fake):
    with pytest.raises(SystemExit, match="--txt or --xlsx"):
        _main(fake, "--to", "uniprot")
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
    assert "failures: fix FAILED" in msg                # never "failures: none"
    assert "FAILED" in capsys.readouterr().out


def test_failed_units_and_incomplete_targets_are_flagged_not_fatal(fake, capsys):
    _calibrate_all(fake.home)
    fake.state["counts"] = {"structures": {"failed": 3},
                            "dock": {"failed": 0, "incomplete_targets": 2}}
    assert _main(fake, "--from", "uniprot") == 0
    assert _runs(fake.events)[-1] == "dock"
    title, msg, prio = fake.notes[-1]
    assert "failures" in title.lower()
    assert "structures failed=3" in msg
    assert "dock incomplete_targets=2" in msg
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
    # each benchmark runs right before its stage's units()
    for s in ("structures", "fix", "similarity"):
        assert fake.events.index(("benchmark", (s,))) < fake.events.index(("units", s))
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
         "--pm-bin-dir", "/opt/pm", "--pmax", "0.5", "--exclude", "--hub-tie-break", "union",
         "--library", "inputs/x.sdf", "--seed", "7", "--chunk-size", "5",
         "--threads-per-job", "2", "--chunk-timeout", "100"])
    assert args.txt == "inputs/ids.txt" and args.pmax == 0.5 and args.exclude == []
    assert args.library == "inputs/x.sdf" and args.chunk_size == 5 and args.seed == 7
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
