"""`bin/colab/fetch.py` — the VM side of the incremental pull.

Before this existed the round trip packed everything once, at the very end, so a reclaimed
Colab VM cost every unit drawn since the session started. These tests pin the selection
rules the host depends on: only finished units travel, only ones the host lacks, and the
marker that certifies a unit enters the archive after the record it certifies.
"""

import json
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

FETCH = Path(__file__).resolve().parents[1] / "bin" / "colab" / "fetch.py"


def _unit(bundle, label, rep, done=True, extras=True):
    """Write one unit into the bundle's results tree, as the worker would."""
    d = bundle / "results" / label / "replicates"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"r{rep}.sdf").write_text(f"{label} r{rep}\n$$$$\n")
    if extras:
        (d / f"r{rep}.log").write_text("engine log\n")
        (d / f"r{rep}.prov.json").write_text(json.dumps({"replicate": rep}))
        (d / f"r{rep}.filtered.tsv").write_text("id\treason\n")
    if done:
        # written last by the worker, and it holds the unit key
        (d / f"r{rep}.done").write_text(f"key-{label}-{rep}\n")
    return f"{label}/r{rep}"


def run_fetch(tmp_path, bundle, have=(), force=False):
    """Run fetch.py against `bundle`; returns (stdout, path to the archive or None)."""
    have_file = tmp_path / "have.txt"
    have_file.write_text("".join(f"{h}\n" for h in have))
    out = tmp_path / "outbox.tgz"
    env = {"FETCH_BUNDLE": str(bundle), "FETCH_HAVE": str(have_file),
           "FETCH_OUT": str(out), "PATH": "/usr/bin:/bin"}
    if force:
        env["FETCH_FORCE"] = "1"
    p = subprocess.run([sys.executable, "-I", str(FETCH)], capture_output=True, text=True,
                       env=env, timeout=120)
    assert p.returncode == 0, p.stderr
    return p.stdout, (out if out.exists() else None)


def members(archive):
    with tarfile.open(archive) as tar:
        return tar.getnames()


@pytest.fixture
def bundle(tmp_path):
    b = tmp_path / "bundle"
    (b / "results").mkdir(parents=True)
    return b


def test_an_empty_bundle_sends_nothing_and_does_not_fail(tmp_path, bundle):
    """pack.py raised SystemExit here, which would abort the wrapper on the first poll."""
    out, archive = run_fetch(tmp_path, bundle)
    assert "FETCH units=0" in out
    assert archive is None


def test_only_units_with_a_marker_travel(tmp_path, bundle):
    """A unit still being written has no .done yet. Sending it would put an uncertified
    record on the host, and the next fetch would have to send it again."""
    _unit(bundle, "A", 0, done=True)
    _unit(bundle, "B", 0, done=False)          # mid-write
    out, archive = run_fetch(tmp_path, bundle)
    assert "FETCH units=1" in out
    names = members(archive)
    assert any("A/replicates/r0.sdf" in n for n in names)
    assert not any("/B/" in n for n in names)


def test_a_unit_the_host_already_holds_is_not_sent_again(tmp_path, bundle):
    _unit(bundle, "A", 0)
    _unit(bundle, "A", 1)
    out, archive = run_fetch(tmp_path, bundle, have=["A/r0"])
    assert "FETCH units=1" in out
    names = members(archive)
    assert any("r1.sdf" in n for n in names)
    assert not any("r0.sdf" in n for n in names)


def test_the_marker_enters_the_archive_after_the_record(tmp_path, bundle):
    """An interrupted extraction must not leave a certificate without its record: the
    host's unit_state would then read DONE for a unit whose SDF is missing."""
    _unit(bundle, "A", 0)
    _out, archive = run_fetch(tmp_path, bundle)
    names = [n for n in members(archive) if "/A/replicates/" in n]
    assert names[-1].endswith("r0.done")
    assert names.index("results/A/replicates/r0.sdf") < names.index("results/A/replicates/r0.done")


def test_every_file_of_a_unit_travels(tmp_path, bundle):
    _unit(bundle, "A", 0)
    _out, archive = run_fetch(tmp_path, bundle)
    got = {Path(n).name for n in members(archive)}
    assert got == {"r0.sdf", "r0.log", "r0.prov.json", "r0.filtered.tsv", "r0.done"}


def test_a_unit_missing_optional_files_still_travels(tmp_path, bundle):
    """filtered.tsv is absent for a unit the engine fully posed; that is not a failure."""
    _unit(bundle, "A", 0, extras=False)
    out, archive = run_fetch(tmp_path, bundle)
    assert "FETCH units=1" in out
    assert {Path(n).name for n in members(archive)} == {"r0.sdf", "r0.done"}


def test_the_worker_log_and_failures_come_home(tmp_path, bundle):
    """pack.py sent only results/, so the merged manifest recorded engine_worker
    "unknown" and the worker's own unit failures never reached the host at all."""
    _unit(bundle, "A", 0)
    (bundle / "worker.json").write_text(json.dumps({"engine_name": "unidock"}))
    (bundle / "failed.tsv").write_text("item\tstage\treason\nA/r9\tdock\tboom\n")
    (bundle / "logs").mkdir()
    (bundle / "logs" / "dock-worker.log").write_text("worker said things\n")
    (bundle / "results" / "discarded.tsv").write_text("unit\treason\n")
    _out, archive = run_fetch(tmp_path, bundle)
    names = members(archive)
    for want in ("worker.json", "failed.tsv", "logs/dock-worker.log",
                 "results/discarded.tsv"):
        assert want in names, want


def test_force_sends_the_extras_even_with_no_new_units(tmp_path, bundle):
    """The final sweep: worker.json only exists once the worker has finished, which is
    after the last unit was already pulled."""
    _unit(bundle, "A", 0)
    (bundle / "worker.json").write_text(json.dumps({"engine_name": "unidock"}))
    out, archive = run_fetch(tmp_path, bundle, have=["A/r0"], force=True)
    assert "FETCH units=0" in out
    assert archive is not None and "worker.json" in members(archive)


def test_a_stale_archive_is_removed_when_there_is_nothing_to_send(tmp_path, bundle):
    """Otherwise the host could download the previous poll's archive a second time."""
    _unit(bundle, "A", 0)
    _out, archive = run_fetch(tmp_path, bundle)
    assert archive is not None and archive.exists()
    out2, archive2 = run_fetch(tmp_path, bundle, have=["A/r0"])
    assert "FETCH units=0" in out2
    assert archive2 is None and not archive.exists()


def test_the_unit_count_is_reported_for_every_new_unit(tmp_path, bundle):
    for i in range(3):
        _unit(bundle, "A", i)
    _unit(bundle, "B", 0)
    out, _archive = run_fetch(tmp_path, bundle, have=["A/r0"])
    assert "FETCH units=3" in out
    assert "A/r1" in out and "B/r0" in out
