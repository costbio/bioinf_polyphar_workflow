"""`bin/polyphar-colab` — the host side of the Colab round trip.

This file exists because that script had no tests at all. It orchestrates ~20 GPU-hours
of work per session and every failure mode it has was found in production, one costly
surprise at a time: a poll client that hung for 8 h, a wrapper that refused to re-attach
to a live worker, and a pull that only happened at the very end so a reclaimed VM cost
everything drawn since the session started.

The `colab` CLI is replaced by `tests/fixtures/fake_colab`, backed by a directory that
plays the VM, and `polyphar-docker` by a recording stub. The real `bin/colab/fetch.py`
runs inside the fake, so the selection logic under test is the shipped one.
"""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "bin" / "polyphar-colab"
FAKE_COLAB = REPO / "tests" / "fixtures" / "fake_colab"

UNITS = [("A", 0), ("A", 1), ("B", 0), ("B", 1)]


def _git(repo, *args):
    subprocess.run(["git", *args], cwd=repo, check=True, capture_output=True)


@pytest.fixture
def world(tmp_path):
    """A run folder, a clean repo the wrapper can pin a commit from, a fake colab CLI and
    a fake polyphar-docker. Returns a namespace of paths plus a `run()` helper."""
    run_dir = tmp_path / "run"
    (run_dir / "07_docking").mkdir(parents=True)

    # A clean git repo: the wrapper refuses to start when polyphar/ is dirty, because the
    # worker installs polyphar from the commit it pins.
    repo = tmp_path / "repo"
    (repo / "polyphar").mkdir(parents=True)
    (repo / "polyphar" / "__init__.py").write_text("")
    shutil.copytree(REPO / "bin" / "colab", repo / "bin" / "colab")
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "t@example.com")
    _git(repo, "config", "user.name", "t")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-qm", "fixture")

    vm = tmp_path / "vm"
    vm.mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    calls = tmp_path / "docker_calls.txt"
    colab_log = tmp_path / "colab_calls.txt"

    # Stands in for bin/polyphar-docker: dock-export writes a bundle, dock-merge counts
    # what arrived. The polyphar internals have their own tests; what is under test here
    # is the orchestration around them.
    (bin_dir / "polyphar-docker").write_text(f"""#!/usr/bin/env python3
import json, os, sys
from pathlib import Path
argv = sys.argv[1:]
open({str(calls)!r}, "a").write(" ".join(argv) + "\\n")
run_dir = Path(os.environ["POLYPHAR_RUN_DIR"])
rel = argv[argv.index("--bundle") + 1]
bundle = run_dir / rel
if argv[0] == "dock-export":
    (bundle / "results").mkdir(parents=True, exist_ok=True)
    (bundle / "receptors").mkdir(exist_ok=True)
    rows = ["label\\treplicate\\tseed\\tkey"]
    for label, rep in {UNITS!r}:
        rows.append(f"{{label}}\\t{{rep}}\\t42\\tkey-{{label}}-{{rep}}")
    (bundle / "units.tsv").write_text("\\n".join(rows) + "\\n")
    (bundle / "run_key.json").write_text(json.dumps({{"units": len({UNITS!r})}}))
    (bundle / "library.sdf").write_text("$$$$\\n")
    sys.exit(0)
if argv[0] == "dock-merge":
    done = sorted(str(p) for p in (bundle / "results").glob("*/replicates/*.done"))
    (run_dir / "merged.json").write_text(json.dumps(
        {{"merged": len(done), "units": [Path(p).parent.parent.name + "/" + Path(p).stem
                                        for p in done]}}))
    sys.exit(0 if len(done) == len({UNITS!r}) else 3)
sys.exit(2)
""")
    (bin_dir / "polyphar-docker").chmod(0o755)
    (bin_dir / "colab").symlink_to(FAKE_COLAB)
    (bin_dir / "docker").write_text("#!/bin/sh\nexit 0\n")      # image inspect passes
    (bin_dir / "docker").chmod(0o755)
    shutil.copyfile(WRAPPER, bin_dir / "polyphar-colab")
    (bin_dir / "polyphar-colab").chmod(0o755)

    def run(_timeout=180, _argv=(), **env):
        e = dict(os.environ)
        e.update({
            "PATH": f"{bin_dir}:{os.environ['PATH']}",
            "POLYPHAR_RUN_DIR": str(run_dir),
            "POLYPHAR_REPO": str(repo),
            "POLYPHAR_HOME": str(tmp_path / "home"),
            "POLYPHAR_VERSION": "test",
            "COLAB_POLL_SECONDS": "0",
            "FAKE_COLAB_VM": str(vm),
            "FAKE_COLAB_LOG": str(colab_log),
        })
        e.update({k: str(v) for k, v in env.items()})
        return subprocess.run([str(bin_dir / "polyphar-colab"), *_argv],
                              capture_output=True, text=True, env=e, timeout=_timeout)

    ns = pytest.importorskip("types").SimpleNamespace(
        run=run, run_dir=run_dir, vm=vm, repo=repo, calls=calls, colab_log=colab_log,
        tmp=tmp_path)
    return ns


def bundles(run_dir):
    return sorted((run_dir / "07_docking" / "export").glob("*/"))


def held(run_dir):
    """Units pulled into the local bundle."""
    out = []
    for b in bundles(run_dir):
        for p in (b / "results").glob("*/replicates/*.done"):
            out.append(f"{p.parent.parent.name}/{p.stem}")
    return sorted(out)


def colab_calls(world):
    try:
        return world.colab_log.read_text().splitlines()
    except OSError:
        return []


# --- the happy path -------------------------------------------------------------------

def test_a_full_round_trip_pulls_every_unit_and_merges(world):
    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert held(world.run_dir) == ["A/r0", "A/r1", "B/r0", "B/r1"], p.stdout + p.stderr
    merged = json.loads((world.run_dir / "merged.json").read_text())
    assert merged["merged"] == 4
    assert (bundles(world.run_dir)[0] / ".merged").exists()


def test_units_arrive_during_the_run_not_only_at_the_end(world):
    """The whole point. Each poll pulls what has finished, so a VM reclaimed mid-session
    costs one interval instead of everything drawn since it started."""
    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    pulls = [l for l in p.stdout.splitlines() if "pulled" in l]
    assert len(pulls) >= 2, p.stdout


def test_a_unit_the_host_already_holds_is_never_fetched_twice(world):
    world.run()
    fetched = []
    for line in colab_calls(world):
        if "fetch.py" in line:
            fetched.append(line)
    # one fetch per poll that had something new, plus the final forced sweep
    assert len(fetched) <= len(UNITS) + 1


def test_a_poll_with_nothing_new_costs_no_transfer(world):
    """When the VM has nothing this host lacks there must be no upload, no exec and no
    download -- that is every poll but the few where a unit actually lands."""
    p = world.run(FAKE_COLAB_UNITS_PER_POLL=0, FAKE_COLAB_GONE_AFTER=3)
    assert p.returncode == 0, p.stdout + p.stderr
    # only the final forced sweep may fetch; the three quiet polls must not
    assert len([l for l in colab_calls(world) if "fetch.py" in l]) == 1
    assert len([l for l in colab_calls(world) if l.startswith("upload have")
                or "/content/have.txt" in l]) == 1


def test_a_lost_session_ends_the_run_instead_of_polling_for_ever(world):
    """Found by these tests: the loop polled a dead session indefinitely. Harmless when
    nothing had been pulled; now that every poll brings results home, the right end to a
    reclaimed VM is to stop and merge what is here."""
    p = world.run(FAKE_COLAB_DIE_AFTER=2, COLAB_MAX_POLL_FAILS=3, _timeout=60)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "giving up after" in p.stdout
    pulled = held(world.run_dir)
    assert pulled, "the units drawn before the VM died should have been pulled"
    merged = json.loads((world.run_dir / "merged.json").read_text())
    assert merged["units"] == pulled          # exactly what was pulled, merged
    for uid in pulled:                        # and intact
        label, name = uid.split("/")
        sdf = bundles(world.run_dir)[0] / "results" / label / "replicates" / f"{name}.sdf"
        assert sdf.read_text().endswith("$$$$\n")


def test_rerunning_after_a_kill_reattaches_to_the_live_worker(world):
    """The wrapper's terminal dies; the detached worker does not. Re-running must pick
    the same bundle back up rather than re-export or refuse."""
    with pytest.raises(subprocess.TimeoutExpired):
        world.run(FAKE_COLAB_HANG_AFTER=2, _timeout=6)   # never reports GONE
    first = held(world.run_dir)
    assert first, "nothing was pulled before the kill"
    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert "re-attaching to bundle" in p.stdout
    assert len(bundles(world.run_dir)) == 1               # nothing re-exported
    assert held(world.run_dir) == ["A/r0", "A/r1", "B/r0", "B/r1"]


def test_rerunning_when_the_session_is_gone_merges_what_was_pulled(world):
    """Since the pull became incremental, "results here, no session" is the ordinary
    state after an interruption. Refusing would strand real results behind a manual
    command."""
    with pytest.raises(subprocess.TimeoutExpired):
        world.run(FAKE_COLAB_HANG_AFTER=2, _timeout=6)
    pulled = held(world.run_dir)
    assert pulled
    (world.vm / ".sessions" / "dock").unlink()            # the VM is reclaimed meanwhile
    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert "merging" in p.stdout
    merged = json.loads((world.run_dir / "merged.json").read_text())
    assert merged["units"] == pulled
    assert (bundles(world.run_dir)[0] / ".merged").exists()


def test_a_failed_download_is_retried_and_loses_nothing(world):
    p = world.run(FAKE_COLAB_FAIL_DOWNLOAD="1")
    assert p.returncode == 0, p.stdout + p.stderr
    assert held(world.run_dir) == ["A/r0", "A/r1", "B/r0", "B/r1"]
    assert "retrying next poll" in p.stdout


def test_a_failed_fetch_exec_is_retried_and_loses_nothing(world):
    p = world.run(FAKE_COLAB_FAIL_EXEC="2")
    assert p.returncode == 0, p.stdout + p.stderr
    assert held(world.run_dir) == ["A/r0", "A/r1", "B/r0", "B/r1"]


def test_only_one_unmerged_bundle_is_allowed(world):
    """Two bundles exported before either is merged hold the same units: both get worked
    and the second merge's draws are discarded, with the GPU time already spent."""
    export = world.run_dir / "07_docking" / "export"
    for name in ("20000101T000000", "20000102T000000"):
        (export / name / "results").mkdir(parents=True)
        (export / name / "units.tsv").write_text("label\treplicate\n")
    p = world.run()
    assert p.returncode != 0
    assert "more than one unmerged bundle" in p.stdout + p.stderr


def test_the_wrapper_refuses_a_dirty_package(world):
    """The worker installs polyphar from the pinned commit, so a dirty tree would run
    different code on the GPU than this host planned the unit keys with."""
    (world.repo / "polyphar" / "__init__.py").write_text("# changed\n")
    p = world.run()
    assert p.returncode != 0
    assert "uncommitted changes" in p.stdout + p.stderr


def test_the_library_reaches_both_export_and_merge(world):
    """`run --library X` used to lose X at the round trip: polyphar-colab called
    dock-export and dock-merge without it, so they fell back to the DEFAULT library.
    Either it does not exist and the run dies at dock, or -- worse -- it does, and the
    GPU docks a different library than the run planned. The library's sha256 is part of
    every unit key, so those results would be keyed to something the run never used."""
    p = world.run(_argv=["--library", "inputs/custom.sdf"])
    assert p.returncode == 0, p.stdout + p.stderr
    calls = world.calls.read_text().splitlines()
    exports = [c for c in calls if c.startswith("dock-export")]
    merges = [c for c in calls if c.startswith("dock-merge")]
    assert exports and all("--library inputs/custom.sdf" in c for c in exports), calls
    assert merges and all("--library inputs/custom.sdf" in c for c in merges), calls


def test_a_refused_launch_is_not_mistaken_for_a_started_worker(world):
    """`colab exec` exits 0 even when the script it ran raised -- the cell executed, which
    is all the CLI reports. Trusting that, the wrapper polled a worker that was never
    started, saw NO-PIDFILE at once, merged zero units and said "done". Found on the real
    Colab path: seven subtype runs reported success having docked nothing."""
    p = world.run(FAKE_COLAB_LAUNCH_REFUSES="1")
    assert p.returncode != 0, p.stdout + p.stderr
    assert "the worker did not start" in p.stdout + p.stderr
    assert not (world.run_dir / "merged.json").exists(), "it merged anyway"


def test_a_failed_engine_install_is_not_mistaken_for_a_ready_vm(world):
    """Same trap one step earlier: without the engine the worker cannot dock at all."""
    p = world.run(FAKE_COLAB_SETUP_REFUSES="1")
    assert p.returncode != 0, p.stdout + p.stderr
    assert "did not install" in p.stdout + p.stderr
