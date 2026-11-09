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
# A unit already recorded in the run folder is never exported again. The real
# dock-export re-plans what is outstanding, and the wrapper's recovery path depends on
# it: after merging a stranded bundle the next export must carry only the rest.
rec_file = run_dir / "recorded.json"
def _recorded():
    return set(json.loads(rec_file.read_text())) if rec_file.is_file() else set()
if argv[0] == "dock-export":
    outstanding = [u for u in {UNITS!r} if f"{{u[0]}}/r{{u[1]}}" not in _recorded()]
    if not outstanding:
        # the real dock-export still writes units.tsv, header only; the wrapper reads the
        # row count to decide "nothing outstanding" and treats a MISSING file as a failure
        bundle.mkdir(parents=True, exist_ok=True)
        (bundle / "units.tsv").write_text("label\\treplicate\\tseed\\tkey\\n")
        print("dock-export: 0 units (no work)")
        sys.exit(0)
    (bundle / "results").mkdir(parents=True, exist_ok=True)
    (bundle / "receptors").mkdir(exist_ok=True)
    rows = ["label\\treplicate\\tseed\\tkey"]
    for label, rep in outstanding:
        rows.append(f"{{label}}\\t{{rep}}\\t42\\tkey-{{label}}-{{rep}}")
    (bundle / "units.tsv").write_text("\\n".join(rows) + "\\n")
    (bundle / "run_key.json").write_text(json.dumps({{"units": len(outstanding)}}))
    (bundle / "library.sdf").write_text("$$$$\\n")
    sys.exit(0)
if argv[0] == "dock-merge":
    done = sorted(str(p) for p in (bundle / "results").glob("*/replicates/*.done"))
    ids = [Path(p).parent.parent.name + "/" + Path(p).stem for p in done]
    recorded = set(ids) | _recorded()
    rec_file.write_text(json.dumps(sorted(recorded)))
    (run_dir / "merged.json").write_text(json.dumps(
        {{"merged": len(done), "units": ids, "recorded": sorted(recorded)}}))
    sys.exit(0 if len(recorded) == len({UNITS!r}) else 3)
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
    """Units the HOST holds, across every bundle.

    A unit can sit in a bundle stranded by an interruption and again in a later one; it
    is still one unit, so this is a set.
    """
    out = set()
    for b in bundles(run_dir):
        for p in (b / "results").glob("*/replicates/*.done"):
            out.add(f"{p.parent.parent.name}/{p.stem}")
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
    # and it says WHY it concluded that: Colab was asked and answered
    assert "giving up: Colab answered" in p.stdout, p.stdout
    assert "is not there" in p.stdout, p.stdout
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
    assert set(pulled) <= set(merged["recorded"]), "the stranded units must be recorded"
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


# --- the one-shot cells: live output, bounded wait --------------------------------------

def test_the_setup_cell_is_bounded_by_its_own_budget_not_by_the_production_one(world):
    """`colab exec`'s --timeout bounds the remote cell, never the client, so the client
    side is the only real bound -- and it was hardcoded at 4200 s.

    On 2026-11-06 a subtype run sat inside that budget for 32 minutes with the VM silent.
    Seventy minutes is the right patience for a production session (a VM you have already
    waited ten minutes for is worth more than a fast failure); it is the wrong patience
    for anything else, and there was no way to say so.
    """
    import time

    t0 = time.monotonic()
    p = world.run(_timeout=90, COLAB_SETUP_TIMEOUT="3", FAKE_COLAB_SETUP_STALLS="60")
    waited = time.monotonic() - t0

    both = p.stdout + p.stderr
    assert p.returncode != 0, both
    assert waited < 45, f"waited {waited:.0f}s for a 3 s budget"
    # and it must SAY so. `set -e` with `pipefail` ends the script the moment `timeout`
    # returns 124, so the first version of this exited 124 with no message at all -- the
    # same silence, just arriving sooner.
    assert "did not install" in both, both
    assert "COLAB_SETUP_TIMEOUT" in both, both
    assert "nothing was lost" in both, both


def test_a_stalled_setup_still_shows_what_the_vm_was_doing(world):
    """The VM is chatty even while it is slow: setup.py prints each install command before
    running it. `setup_out="$(colab exec ...)"` captured that stream and showed none of it
    until the cell ended, so a slow micromamba and a reclaimed VM looked identical for up
    to seventy minutes. colab_cli flushes every chunk as it arrives -- the output was being
    produced and thrown away.
    """
    p = world.run(_timeout=90, COLAB_SETUP_TIMEOUT="3", FAKE_COLAB_SETUP_STALLS="60")
    assert "micromamba install" in p.stdout, p.stdout + p.stderr


def test_the_setup_and_launch_cells_still_reach_the_terminal_on_the_happy_path(world):
    """tee must not eat what the old `printf '%s\\n' "$out"` used to echo."""
    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert "setup: ready" in p.stdout, p.stdout
    assert "worker pid=" in p.stdout, p.stdout


def test_the_remote_cell_gives_up_before_the_client_does(world):
    """If the cell outlived the client the socket would be cut with nothing to show; if it
    were much shorter the client's budget would be a lie. One margin rule, both cells."""
    out = subprocess.run(
        ["bash", "-c",
         f'set -euo pipefail; source <(sed -n "/^cell_timeout() {{/,/^}}/p" {WRAPPER!s}); '
         'for t in 4200 900 300 60; do echo "$t $(cell_timeout $t)"; done'],
        capture_output=True, text=True, check=True).stdout.split()
    pairs = [(int(out[i]), int(out[i + 1])) for i in range(0, len(out), 2)]
    for client, cell in pairs:
        assert 0 < cell < client, (client, cell)
        assert cell >= client * 0.8, f"{cell} is a needlessly small share of {client}"


# --- a lost websocket is not a refusal --------------------------------------------------

def test_a_lost_connection_during_setup_is_asked_again(world):
    """What cost Cluster4 on 2026-11-06.

    `colab exec` exits 0 when the cell RAN and the script raised -- a refusal, and
    definitive. A non-zero exit means the CLI itself gave up: on POSEIDON the websocket
    drops with `RuntimeError: Connection was lost` and the client then hangs until the
    budget kills it, while the kernel carries on installing server-side. setup.py is
    idempotent -- it reports "unidock already installed" and returns in seconds -- so the
    honest answer to a lost connection is to ask again, not to throw the session away.
    """
    p = world.run(FAKE_COLAB_LOSE_SETUP="1")
    assert p.returncode == 0, p.stdout + p.stderr
    assert "asking again" in p.stdout, p.stdout
    assert held(world.run_dir) == ["A/r0", "A/r1", "B/r0", "B/r1"], p.stdout + p.stderr


def test_a_refused_setup_is_never_retried(world):
    """A refusal is the cell's own answer and asking twice only wastes the budget."""
    p = world.run(_timeout=120, FAKE_COLAB_SETUP_REFUSES="1")
    assert p.returncode != 0, p.stdout
    assert "asking again" not in p.stdout, p.stdout
    setups = [l for l in colab_calls(world) if "setup.py" in l]
    assert len(setups) == 1, setups


def test_a_lost_connection_during_launch_asks_the_vm_instead_of_relaunching(world):
    """launch.py must never be re-run blind.

    If the first call DID start the worker -- and on a lost connection it usually did, the
    client just never saw it say so -- a second call is refused with "a dock-worker is
    already running" and the wrapper dies with a worker on the VM burning GPU units that
    nobody is watching. The poll already answers the only question that matters.
    """
    p = world.run(FAKE_COLAB_LOSE_LAUNCH="1")
    assert p.returncode == 0, p.stdout + p.stderr
    assert "asking the VM whether it started" in p.stdout, p.stdout
    launches = [l for l in colab_calls(world) if "launch.py" in l]
    assert len(launches) == 1, f"launch.py must not be re-run: {launches}"
    assert held(world.run_dir) == ["A/r0", "A/r1", "B/r0", "B/r1"], p.stdout + p.stderr


def test_the_dock_options_reach_both_ends_of_the_round_trip(world):
    """`--replicates 5` must mean five everywhere, not three on the GPU host.

    The wrapper forwards them after a `--`; both dock-export and dock-merge have to see
    the same ones, because the first decides which units are outstanding and the second
    has to agree about what it is folding back in.
    """
    p = world.run(_argv=("--", "--replicates", "5", "--seed", "7"))
    assert p.returncode == 0, p.stdout + p.stderr
    calls = world.calls.read_text().splitlines()
    for stage in ("dock-export", "dock-merge"):
        line = next(l for l in calls if l.startswith(stage))
        assert "--replicates 5" in line, (stage, line)
        assert "--seed 7" in line, (stage, line)


def test_the_planned_engine_string_matches_what_the_binary_calls_itself(world):
    """dock-merge warns when the engine that DREW the units is not the one their keys were
    planned under. That warning is the only thing standing between a silently changed
    engine on the VM and a set of results keyed to a version that never ran.

    The wrapper planned under 'Uni-Dock 1.2.0' while the binary reports 'Uni-Dock v1.2.0',
    so the warning fired on every merge and meant nothing -- a check that always fires is
    worse than no check, because people learn to scroll past it.

    The string must stay pinned here rather than read from the engine: this host has no
    GPU to ask, and a unit key cannot be assigned after the work it identifies.
    """
    import re as _re

    text = WRAPPER.read_text()
    assert _re.search(r'engine_version="\$\{POLYPHAR_ENGINE_VERSION:-([^}]+)\}"', text), \
        "the planned engine string must be one named, overridable value"
    planned = _re.search(r'engine_version="\$\{POLYPHAR_ENGINE_VERSION:-([^}]+)\}"',
                         text).group(1)
    # what bin/colab/setup.py installs and the binary then reports of itself
    assert planned == "Uni-Dock v1.2.0", planned
    assert "--engine-version \"$engine_version\"" in text


def test_a_refused_launch_is_still_fatal(world):
    """The probe only covers a lost connection. A refusal exits 0 with a reason, and the
    worker really is not running -- that must still stop the run."""
    p = world.run(_timeout=120, FAKE_COLAB_LAUNCH_REFUSES="1")
    assert p.returncode != 0, p.stdout
    assert "the worker did not start" in p.stdout + p.stderr, p.stdout + p.stderr
    assert "asking the VM" not in p.stdout, p.stdout


def test_the_session_name_is_reachable_without_touching_the_argv(world):
    """Concurrent subtypes need one session each, and the only caller is polyphar-docker,
    which has no --session to forward.

    `bin/polyphar-docker dock` builds the round-trip argv itself; a user cannot get
    `--session` into it. Every run folder therefore reused the session called "dock", so a
    second subtype started while the first was working would re-attach to the first one's
    VM, find a worker already running and refuse. An environment variable is reachable from
    anywhere, which is what the run folders need.
    """
    p = world.run(COLAB_SESSION="dock2", COLAB_GPU="A100")
    assert p.returncode == 0, p.stdout + p.stderr
    calls = colab_calls(world)
    assert calls, "the fake recorded nothing"
    assert all("-s dock2" in c or "dock2" in c for c in calls if " -s " in c), calls
    assert any("--gpu A100" in c for c in calls), calls
    assert not any(" -s dock " in f" {c} " for c in calls), \
        f"the default session must not be used when COLAB_SESSION is set: {calls}"


def test_an_explicit_session_argument_still_wins_over_the_variable(world):
    """The flag is the more specific statement of intent."""
    p = world.run(_argv=("--session", "fromflag"), COLAB_SESSION="fromenv")
    assert p.returncode == 0, p.stdout + p.stderr
    calls = [c for c in colab_calls(world) if " -s " in c]
    assert calls and all("fromflag" in c for c in calls), calls


def test_a_recovery_merge_carries_on_instead_of_telling_you_to_rerun(world):
    """A stranded bundle leaves the run exactly where a fresh invocation would start.

    On 2026-11-06 the first command of a session found a bundle left by a run killed
    hours earlier, merged its 5 units and exited -- correct, but it read as "did nothing"
    and the docking the user asked for had not begun. There is no decision in that step,
    so the wrapper takes it: merge the stranded bundle, then export the next one.
    """
    # strand a bundle: pull some units, then lose the VM
    with pytest.raises(subprocess.TimeoutExpired):
        world.run(FAKE_COLAB_HANG_AFTER=2, _timeout=6)
    stranded = held(world.run_dir)
    assert stranded, "nothing was pulled before the interruption"
    (world.vm / ".sessions" / "dock").unlink()

    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert "merging" in p.stdout
    assert "carrying on to the next one" in p.stdout, p.stdout
    # it really exported again rather than stopping at the merge
    assert len(bundles(world.run_dir)) == 2, [b.name for b in bundles(world.run_dir)]
    assert held(world.run_dir) == ["A/r0", "A/r1", "B/r0", "B/r1"]


def test_the_continuation_happens_once_and_cannot_spin(world):
    """A run that keeps stranding bundles must stop and be looked at, not loop."""
    with pytest.raises(subprocess.TimeoutExpired):
        world.run(FAKE_COLAB_HANG_AFTER=2, _timeout=6)
    (world.vm / ".sessions" / "dock").unlink()

    p = world.run(POLYPHAR_COLAB_CONTINUED="1")
    assert p.returncode == 0, p.stdout + p.stderr
    assert "merging" in p.stdout
    assert "carrying on to the next one" not in p.stdout, p.stdout
    assert len(bundles(world.run_dir)) == 1, "it must not export again when already continued"


def test_an_ordinary_run_does_not_announce_a_continuation(world):
    """The continuation belongs to the recovery path only: a normal round trip ends when
    its own bundle is merged."""
    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert "carrying on to the next one" not in p.stdout, p.stdout
    assert len(bundles(world.run_dir)) == 1


# --- a failed poll is not a dead session ------------------------------------------------

def test_polls_failing_while_the_session_lives_does_not_abandon_it(world):
    """What happened on 2026-11-07.

    POSEIDON lost its network. The colab CLI could not refresh its token and fell back to
    an interactive `Enter the authorization code:` prompt nobody could answer, so every
    poll failed. After ten the wrapper announced "the session is gone", merged, wrote
    `.merged` and exited -- while that VM was still docking at 52 %, and the closed bundle
    then blocked the re-attach that would have recovered it.

    A failed poll says the CLIENT could not ask. Ask Colab before concluding.
    """
    # it must NOT end: a live session with a broken link is a reason to keep polling, so
    # the run is bounded here rather than by the wrapper
    with pytest.raises(subprocess.TimeoutExpired) as exc:
        world.run(_timeout=25, FAKE_COLAB_POLL_DOWN="1", COLAB_MAX_POLL_FAILS="2",
                  COLAB_POLL_SECONDS="0")
    both = (exc.value.stdout or b"").decode() + (exc.value.stderr or b"").decode()
    assert "the session is still there" in both, both
    assert "is not there" not in both, both
    assert "giving up" not in both, both
    assert not (bundles(world.run_dir)[0] / ".merged").exists(), \
        "it must not close the bundle while the session is alive"


def test_when_colab_itself_is_unreachable_the_bundle_is_left_reattachable(world):
    """If `colab sessions` fails too, the session's fate is UNKNOWN. Merge what is home --
    always safe -- but do not close the bundle, because closing it is what made net3
    unrecoverable without moving the marker aside by hand."""
    p = world.run(_timeout=120, FAKE_COLAB_LINK_DOWN="1", COLAB_MAX_POLL_FAILS="2",
                  COLAB_POLL_SECONDS="0")
    both = p.stdout + p.stderr
    assert "cannot reach Colab at all" in both, both
    assert "left OPEN because Colab was unreachable" in both, both
    assert "colab stop -s" in both, both       # the user is told it may still be billing
    b = bundles(world.run_dir)[0]
    assert not (b / ".merged").exists(), "the bundle must stay re-attachable"


def test_a_session_colab_confirms_is_gone_is_closed_normally(world):
    """The third case: Colab answers and the session really is absent. That is a genuine
    loss, the bundle is closed, and the message says Colab confirmed it."""
    with pytest.raises(subprocess.TimeoutExpired):
        world.run(FAKE_COLAB_HANG_AFTER=2, _timeout=6)
    (world.vm / ".sessions" / "dock").unlink()          # really gone
    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert "cannot reach Colab at all" not in p.stdout
    assert (bundles(world.run_dir)[0] / ".merged").exists()


def test_an_empty_stranded_bundle_is_dropped_not_a_dead_end(world):
    """An export whose session never started leaves a bundle with nothing in it.

    On 2026-11-09 Colab answered `Service Unavailable` to net1's L4 assignment request.
    The export had already written a 199-unit bundle, so the run folder was left with an
    unmerged bundle holding zero results -- and every subsequent `polyphar dock` died on
    the unmerged-bundle guard, telling the user to delete it by hand. There is nothing in
    it to lose and export re-plans from what is recorded, so it is dropped, exactly as the
    "nothing outstanding" case already does.
    """
    world.run()                                     # a normal round trip, bundle merged
    first = held(world.run_dir)
    assert first

    # a later export that never got a VM: a bundle with units.tsv and no results
    stranded = world.run_dir / "07_docking" / "export" / "20990101T000000"
    (stranded / "results").mkdir(parents=True)
    (stranded / "units.tsv").write_text("label\treplicate\tseed\tkey\nZ\t0\t42\tkey-Z-0\n")
    (stranded / "run_key.json").write_text('{"units": 1}')

    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert "has no results (its session never started); dropping it" in p.stdout, p.stdout
    assert not stranded.exists(), "the empty bundle must be gone"
    assert held(world.run_dir) == first, "nothing recorded may be disturbed"


def test_a_stranded_bundle_WITH_results_is_still_never_dropped(world):
    """The drop is for empty bundles only: one holding real draws is merged, never
    deleted, because those draws cost GPU time and cannot be re-created."""
    with pytest.raises(subprocess.TimeoutExpired):
        world.run(FAKE_COLAB_HANG_AFTER=2, _timeout=6)
    pulled = held(world.run_dir)
    assert pulled
    (world.vm / ".sessions" / "dock").unlink()
    p = world.run()
    assert p.returncode == 0, p.stdout + p.stderr
    assert "dropping it" not in p.stdout, p.stdout
    assert set(pulled) <= set(held(world.run_dir))
