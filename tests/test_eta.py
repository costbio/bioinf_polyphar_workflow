import json
from datetime import datetime

import pytest

from polyphar.common import WorkUnits
from polyphar.eta import (PARALLEL_EFFICIENCY, Calibration, Estimate, Progress,
                          estimate_stage, format_eta, polyphar_home, slots_for)


def _calib(tmp_path, **costs):
    c = Calibration.load(tmp_path)
    for stage, cost in costs.items():
        c.update(stage, cost)
    return c


def test_estimate_uses_largest_job_bound(tmp_path):
    calib = _calib(tmp_path, fix=1.0)
    est = estimate_stage("fix", WorkUnits(3, "bytes", [100, 1, 1]), 8, calib,
                         PARALLEL_EFFICIENCY)
    assert est.seconds == 100


def test_estimate_total_over_slots_when_balanced(tmp_path):
    calib = _calib(tmp_path, fix=2.0)
    est = estimate_stage("fix", WorkUnits(8, "bytes", [10] * 8), 4, calib,
                         PARALLEL_EFFICIENCY)
    assert est.seconds == 80 * 2.0 / 4


def test_estimate_count_only(tmp_path):
    calib = _calib(tmp_path, pockets=3.0)
    est = estimate_stage("pockets", WorkUnits(10, "structure"), 5, calib,
                         PARALLEL_EFFICIENCY)
    assert est.seconds == 10 * 3.0 / 5


def test_estimate_uncalibrated_none(tmp_path):
    est = estimate_stage("fix", WorkUnits(5, "bytes"), 4, Calibration.load(tmp_path),
                         PARALLEL_EFFICIENCY)
    assert est.seconds is None
    assert est.basis == "uncalibrated — run polyphar benchmark"
    assert format_eta(est, datetime(2026, 1, 1)) == est.basis


def test_estimate_trivial_stages_and_zero_count(tmp_path):
    calib = Calibration.load(tmp_path)
    assert estimate_stage("ids", WorkUnits(500, "ids"), 4, calib,
                          PARALLEL_EFFICIENCY) == Estimate(0.0, "trivial")
    assert estimate_stage("select", WorkUnits(1, "x"), 4, calib,
                          PARALLEL_EFFICIENCY).seconds == 0.0
    assert estimate_stage("fix", WorkUnits(0, "x"), 4, calib,
                          PARALLEL_EFFICIENCY).seconds == 0.0


@pytest.mark.parametrize("stage,nprocs,tpj,expected", [
    ("fix", 8, 1, 8),
    ("pockets", 8, 1, 8),
    ("similarity", 8, 1, 4),
    ("similarity", 1, 1, 1),        # at least one slot
    ("dock", 16, 4, 4),
    ("dock", 2, 4, 1),              # never below one slot
    ("uniprot", 32, 1, 8),
    ("structures", 4, 1, 4),
    ("ids", 4, 1, None),
])
def test_slots_table(stage, nprocs, tpj, expected):
    assert slots_for(stage, nprocs, tpj) == expected


def test_estimate_dock_uses_threads_per_job(tmp_path):
    calib = _calib(tmp_path, dock=10.0)
    est = estimate_stage("dock", WorkUnits(8, "ligand"), 16, calib,
                         PARALLEL_EFFICIENCY, threads_per_job=4)
    assert est.seconds == 8 * 10.0 / 4


def test_format_eta_shape():
    now = datetime(2026, 10, 31, 10, 45)  # a Saturday
    out = format_eta(Estimate(3 * 3600 + 20 * 60, "x"), now)
    assert out == "~3 h 20 min (finish ≈ Sat 14:05)"
    assert format_eta(Estimate(12 * 60, "x"), now).startswith("~12 min (finish")
    assert format_eta(Estimate(30, "x"), now).startswith("~30 s (finish")
    assert format_eta(Estimate(0.0, "trivial"), now) == "instant (trivial)"
    # far-future finish carries a date
    assert "Nov" in format_eta(Estimate(9 * 86400, "x"), now)


def test_calibration_roundtrip_running_mean(tmp_path):
    home = tmp_path / "new" / "home"
    c = Calibration.load(home)
    assert c.cost("fix") is None
    c.update("fix", 2.0)
    c.update("fix", 4.0)
    c.save()
    raw = json.loads((home / "calibration.json").read_text())
    entry = raw["fix"]
    assert entry["seconds_per_unit"] == 3.0
    assert entry["n_observations"] == 2
    assert entry["nproc"] >= 1 and entry["cpu_model"] and entry["updated"]
    assert Calibration.load(home).cost("fix") == 3.0
    assert not list(home.glob("*.tmp"))  # atomic save leaves no temp file


def test_calibration_ignores_bad_data(tmp_path):
    (tmp_path / "calibration.json").write_text("{not json")
    c = Calibration.load(tmp_path)
    assert c.cost("fix") is None
    c.update("fix", 0.0)
    c.update("fix", float("nan"))
    assert c.cost("fix") is None


def test_polyphar_home(monkeypatch, tmp_path):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path))
    assert polyphar_home() == tmp_path
    monkeypatch.delenv("POLYPHAR_HOME")
    assert polyphar_home().name == ".polyphar"


class Rec:
    def __init__(self):
        self.calls = []

    def __call__(self, title, message, priority=0):
        self.calls.append((title, message, priority))


class Clock:
    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t


def test_progress_milestones_fire_once():
    rec, clock = Rec(), Clock()
    p = Progress("fix", 100, Estimate(100.0, "x"), rec, clock=clock)
    for _ in range(100):
        clock.t += 1
        p.advance()
    assert len(rec.calls) == 4
    assert [c[0] for c in rec.calls] == ["fix 5%", "fix 25%", "fix 50%", "fix 75%"]
    assert "revised ETA" in rec.calls[0][1]
    assert "revised ETA" not in rec.calls[1][1]


def test_progress_big_jump_fires_each_milestone_once():
    rec = Rec()
    p = Progress("fix", 100, Estimate(None, "x"), rec)
    p.advance(60)
    p.advance(40)
    assert len(rec.calls) == 4


def test_progress_live_eta_in_milestone():
    rec, clock = Rec(), Clock()
    p = Progress("fix", 100, Estimate(None, "x"), rec, clock=clock)
    clock.t = 600.0
    p.advance(50)       # 50% after 10 min -> ~10 min remaining
    assert "~10 min" in rec.calls[-1][1]


def test_progress_size_weighted():
    rec = Rec()
    p = Progress("fix", 2, Estimate(None, "x"), rec, total_size=100.0)
    p.advance(1, size=50)
    assert [c[0] for c in rec.calls] == ["fix 5%", "fix 25%", "fix 50%"]


def test_progress_log_throttled_and_written(tmp_path, capsys):
    clock = Clock()
    log = tmp_path / "logs" / "fix.log"
    p = Progress("fix", 1000, Estimate(None, "x"), Rec(), milestones=(),
                 log_path=log, clock=clock)
    for _ in range(10):
        clock.t += 5
        p.advance()          # 50 s total: one periodic line
    clock.t += 20
    p.advance()              # 70 s: second periodic line
    progress_lines = [l for l in log.read_text().splitlines() if "/1000 units" in l]
    assert len(progress_lines) == 2
    assert "units" in capsys.readouterr().out


def test_progress_notifier_errors_do_not_propagate():
    def boom(*a, **k):
        raise RuntimeError("push down")

    p = Progress("fix", 10, Estimate(None, "x"), boom)
    p.advance(10)
    p.done({"ok": 10})
    p.fail(ValueError("x"))


def test_progress_done_updates_calibration(tmp_path):
    rec, clock = Rec(), Clock()
    calib = Calibration.load(tmp_path)
    p = Progress("fix", 10, Estimate(None, "x"), rec, calib=calib, slots=4.0, clock=clock)
    clock.t = 20.0
    p.advance(10)
    rec.calls.clear()
    p.done({"fixed": 9, "failed": 1})
    # 20 s wall * 4 slots / 10 units = 8 s per unit on one process
    assert calib.cost("fix") == 8.0
    assert Calibration.load(tmp_path).cost("fix") == 8.0   # saved to disk
    title, msg, prio = rec.calls[-1]
    assert title == "fix done" and prio == 1   # failed > 0: priority 1 (item 13)
    assert "fixed=9" in msg and "failures: 1" in msg and "estimated n/a" in msg


def test_progress_done_nothing_done_keeps_calibration(tmp_path):
    calib = Calibration.load(tmp_path)
    Progress("fix", 10, Estimate(None, "x"), Rec(), calib=calib).done({})
    assert calib.cost("fix") is None


def test_progress_fail_priority_1_with_last_10_log_lines():
    rec, clock = Rec(), Clock()
    p = Progress("fix", 1000, Estimate(None, "x"), rec, milestones=(), clock=clock)
    for _ in range(30):
        clock.t += 61
        p.advance()          # a log line each time
    p.fail(RuntimeError("kaboom"))
    title, msg, prio = rec.calls[-1]
    assert title == "fix FAILED" and prio == 1
    lines = msg.splitlines()
    assert len(lines) == 10
    assert "RuntimeError: kaboom" in lines[-1]


def test_calibration_entry_without_n_observations(tmp_path):
    (tmp_path / "calibration.json").write_text(json.dumps({
        "fix": {"seconds_per_unit": 2.0},
        "dock": {"seconds_per_unit": 3.0, "n_observations": "many"},
    }))
    c = Calibration.load(tmp_path)
    assert c.cost("fix") == 2.0 and c.cost("dock") == 3.0
    c.update("fix", 4.0)                      # treated as 1 prior observation
    assert c.cost("fix") == 3.0 and c.data["fix"]["n_observations"] == 2


def test_calibration_update_survives_bad_in_memory_entry(tmp_path):
    c = Calibration(tmp_path / "calibration.json", {"fix": {"seconds_per_unit": 2.0}})
    c.update("fix", 4.0)
    assert c.cost("fix") == 3.0


def test_progress_done_succeeds_with_legacy_calibration_file(tmp_path):
    (tmp_path / "calibration.json").write_text(
        json.dumps({"fix": {"seconds_per_unit": 2.0}}))
    clock = Clock()
    p = Progress("fix", 4, Estimate(None, "x"), Rec(),
                 calib=Calibration.load(tmp_path), slots=1.0, clock=clock)
    clock.t = 8.0
    p.advance(4)
    p.done({"ok": 4})                          # must not raise
    assert Calibration.load(tmp_path).cost("fix") == 2.0   # (2*1 + 2)/2


def test_progress_done_never_raises_even_if_bookkeeping_breaks(tmp_path):
    class Broken:
        def update(self, *a):
            raise KeyError("n_observations")

    p = Progress("fix", 1, Estimate(None, "x"), Rec(), calib=Broken())
    p.advance()
    p.done({})                                 # logged warning, no exception
    assert "bookkeeping failed" in p._recent[-1]


def test_progress_size_mode_advance_without_size_uses_mean_unit_size(tmp_path):
    rec, clock = Rec(), Clock()
    calib = Calibration.load(tmp_path)
    p = Progress("fix", 4, Estimate(None, "x"), rec, total_size=400.0,
                 calib=calib, slots=1.0, clock=clock)
    clock.t = 40.0
    p.advance(2)                               # no size: 2 * (400/4) = 200 of 400
    assert p.work_done == 200.0
    assert [c[0] for c in rec.calls] == ["fix 5%", "fix 25%", "fix 50%"]
    p.advance(2)
    p.done({})
    assert calib.cost("fix") == 40.0 / 400.0   # calibration no longer skipped
