"""Ruling 8: the ETA gate wired into cli.main for every per-stage subcommand."""
import subprocess
import sys
import types

import pytest

from polyphar import cli
from polyphar.common import WorkUnits
from polyphar.eta import Calibration, Progress


@pytest.fixture
def gate(monkeypatch, tmp_path):
    """Fake 'fix' stage + recorders. events is the ordered call log."""
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    events, notes = [], []
    state = {"fail": False}

    mod = types.ModuleType("fake_gate_stage")

    def units(layout, args):
        events.append("units")
        return WorkUnits(3, "pdb files")

    def run(layout, args, progress):
        # whatever is on stdout NOW must already contain the ETA line
        events.append("run")
        for _ in range(3):
            progress.advance()
        if state["fail"]:
            raise RuntimeError("stage blew up")
        return {"done": 3, "failed": 0}

    mod.units, mod.run = units, run
    monkeypatch.setitem(sys.modules, "fake_gate_stage", mod)
    monkeypatch.setitem(cli.STAGES, "fix", ("fake_gate_stage", "fake"))
    monkeypatch.setattr(cli, "notify",
                        lambda title, message, priority=0, **kw: notes.append((title, message, priority)) or True)

    real_fail = Progress.fail

    def spy_fail(self, exc):
        events.append("fail")
        return real_fail(self, exc)

    monkeypatch.setattr(Progress, "fail", spy_fail)
    return types.SimpleNamespace(events=events, notes=notes, state=state, home=tmp_path / "home")


def _argv(tmp_path, *extra):
    return ["fix", "--run-dir", str(tmp_path / "run"), "--nprocs", "2", *extra]


def test_eta_line_printed_before_run(gate, tmp_path, capsys, monkeypatch):
    seen = {}
    real_run = sys.modules["fake_gate_stage"].run

    def run(layout, args, progress):
        seen["stdout_at_run"] = capsys.readouterr().out
        return real_run(layout, args, progress)

    sys.modules["fake_gate_stage"].run = run
    assert cli.main(_argv(tmp_path)) == 0
    assert gate.events == ["units", "run"]
    out = seen["stdout_at_run"]
    assert "fix: 3 units (pdb files)" in out
    assert "uncalibrated" in out                     # ETA text, before any work
    # start notification went out before run, then milestones, then done
    assert gate.notes[0][0] == "fix starting"
    assert gate.notes[-1][0] == "fix done"


def test_dry_run_never_calls_run(gate, tmp_path, capsys):
    assert cli.main(_argv(tmp_path, "--dry-run")) == 0
    assert gate.events == ["units"]
    assert "fix: 3 units" in capsys.readouterr().out
    assert [n[0] for n in gate.notes] == ["fix starting"]   # start notify still sent
    assert not (tmp_path / "run").exists()                  # nothing written


def test_dry_run_no_notify(gate, tmp_path):
    assert cli.main(_argv(tmp_path, "--dry-run", "--no-notify")) == 0
    assert gate.notes == []


def test_no_notify_suppresses_all_notifications(gate, tmp_path):
    assert cli.main(_argv(tmp_path, "--no-notify")) == 0
    assert gate.events == ["units", "run"]
    assert gate.notes == []


def test_exception_triggers_fail_and_reraises(gate, tmp_path):
    gate.state["fail"] = True
    with pytest.raises(RuntimeError, match="stage blew up"):
        cli.main(_argv(tmp_path))
    assert gate.events[-1] == "fail"
    title, msg, prio = gate.notes[-1]
    assert title == "fix FAILED" and prio == 1 and "stage blew up" in msg
    assert "fix done" not in [n[0] for n in gate.notes]


def test_exception_exits_nonzero(tmp_path):
    code = (
        "import sys, types\n"
        "from polyphar import cli\n"
        "from polyphar.common import WorkUnits\n"
        "m = types.ModuleType('boom_stage')\n"
        "m.units = lambda l, a: WorkUnits(1, 'x')\n"
        "def run(l, a, p): raise RuntimeError('boom')\n"
        "m.run = run\n"
        "sys.modules['boom_stage'] = m\n"
        "cli.STAGES['fix'] = ('boom_stage', 'x')\n"
        f"sys.exit(cli.main(['fix', '--run-dir', {str(tmp_path / 'run')!r}, '--no-notify']))\n"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       env={**__import__("os").environ, "POLYPHAR_HOME": str(tmp_path / "h")})
    assert r.returncode != 0
    assert "boom" in r.stderr


def test_success_updates_calibration(gate, tmp_path):
    assert cli.main(_argv(tmp_path, "--no-notify")) == 0
    assert Calibration.load(gate.home).cost("fix") is not None
    assert (tmp_path / "run" / "logs" / "fix.log").exists()


def test_calibrated_eta_line(gate, tmp_path, capsys):
    c = Calibration.load(gate.home)
    c.update("fix", 600.0)       # 600 s/unit; 3 units on 2 slots = 900 s
    c.save()
    assert cli.main(_argv(tmp_path, "--dry-run", "--no-notify")) == 0
    assert "~15 min (finish" in capsys.readouterr().out
