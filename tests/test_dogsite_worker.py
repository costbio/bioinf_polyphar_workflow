"""Host-side DoGSiteScorer worker (bin/polyphar-dogsite-worker) against the fake binary
tests/fixtures/fake_dogsite, plus the wrapper's worker lifecycle (bin/polyphar-docker)."""
import ast
import json
import os
import socket
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import pytest

from polyphar.common import RunLayout, sha256_file
from polyphar.stages import pockets

REPO = Path(__file__).resolve().parent.parent
WORKER = REPO / "bin" / "polyphar-dogsite-worker"
WRAPPER = REPO / "bin" / "polyphar-docker"
FAKE = REPO / "tests" / "fixtures" / "fake_dogsite"
BENCH = REPO / "polyphar" / "data" / "bench" / "2mjw_B.pdb"
HINT = "license/installation problem — see README §DoGSiteScorer"
INFO_COLUMNS = ["UniProt ID", "Gene Name", "PDB", "Chain ID", "AlphaFold", "Resolution",
                "Structure method", "SIFTS coverage", "File"]


def _make_run(root, names, acc="P25815"):
    """Run folder with fixed PDB chains `names` (e.g. "2mjw_B") listed in
    Protein_info.csv, and their jobs enqueued. Returns the layout."""
    layout = RunLayout(root)
    rows = []
    for name in names:
        pdb, chain = name.split("_", 1)
        for base in (layout.structures, layout.pdbfixer):
            p = base / acc / f"{name}.pdb"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(f"ATOM {name}\n")
        rows.append([acc, "S100P_HUMAN", pdb.upper(), chain, False, "-", "-", "-",
                     f"02_structures/{acc}/{name}.pdb"])
    layout.structures.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=INFO_COLUMNS).to_csv(layout.protein_info, index=False)
    pockets.enqueue(layout)
    return layout


def _worker_cmd(run_dir, dogsite=FAKE, *extra):
    return [sys.executable, str(WORKER), "--run-dir", str(run_dir),
            "--dogsite", str(dogsite), *extra]


def _env(**kw):
    env = {k: v for k, v in os.environ.items()
           if k not in ("FAKE_DOGSITE_SLEEP", "FAKE_DOGSITE_LOG", "DOGSITE_BIN")}
    env.update({k: str(v) for k, v in kw.items()})
    return env


def _queue(layout):
    return sorted(p.name for p in (layout.pockets / "queue").iterdir())


def _logged_models(log):
    """Models from the fake's invocation log, without the preflight runs."""
    lines = Path(log).read_text().splitlines() if Path(log).exists() else []
    return [l for l in lines if Path(l).resolve() != BENCH.resolve()]


def test_worker_runs_queue_with_stub(tmp_path):
    layout = _make_run(tmp_path, ["2mjw_B", "2mjw_C", "1abc_A"])
    r = subprocess.run(_worker_cmd(tmp_path, FAKE, "--nprocs", "2"),
                       env=_env(), capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr
    assert _queue(layout) == ["P25815_1abc_A.done", "P25815_2mjw_B.done",
                              "P25815_2mjw_C.done"]   # no heartbeat
    for name in ("2mjw_B", "2mjw_C", "1abc_A"):
        done = json.loads((layout.pockets / "queue" / f"P25815_{name}.done").read_text())
        assert set(done) == {"exit", "seconds", "dogsite_sha256"}
        assert done["exit"] == 0 and done["seconds"] >= 0
        assert done["dogsite_sha256"] == sha256_file(FAKE)
        out = layout.pockets / "P25815" / f"P25815_{name}.pdb"
        assert Path(f"{out}_desc.txt").is_file()
        assert Path(f"{out}_res_P_1.pdb").read_text().startswith("HEADER")
        assert (layout.pdbfixer / "P25815" / f"{name}.pdb").read_text() == f"ATOM {name}\n"
    df = pockets.parse_pockets(layout.pockets, layout.protein_info)
    assert list(df["Pocket Filename"]) == [f"P25815_{n}.pdb_res_P_1.pdb"
                                           for n in ("2mjw_B", "2mjw_C", "1abc_A")]
    assert (df["drugScore"] >= 0.5).all()     # the fake pocket survives the 0.5 cutoff


def test_worker_failure_recorded(tmp_path):
    layout = _make_run(tmp_path, ["FAIL_A", "2mjw_B"])
    src = layout.pdbfixer / "P25815" / "FAIL_A.pdb"
    before = src.read_bytes()
    r = subprocess.run(_worker_cmd(tmp_path), env=_env(), capture_output=True, text=True,
                       timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr
    assert _queue(layout) == ["P25815_2mjw_B.done", "P25815_FAIL_A.failed"]
    failed = json.loads((layout.pockets / "queue" / "P25815_FAIL_A.failed").read_text())
    assert failed["exit"] == 3
    assert "FAIL" in failed["stderr_tail"]
    assert src.read_bytes() == before
    assert not Path(f"{layout.pockets}/P25815/P25815_FAIL_A.pdb_desc.txt").exists()


def test_two_workers_never_double_claim(tmp_path):
    names = [f"{i}xyz_A" for i in range(40)]
    layout = _make_run(tmp_path, names)
    log = tmp_path / "fake.log"
    env = _env(FAKE_DOGSITE_SLEEP="0.15", FAKE_DOGSITE_LOG=log)
    q = layout.pockets / "queue"
    procs = [subprocess.Popen(_worker_cmd(tmp_path, FAKE, "--nprocs", "3"), env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for _ in range(2)]
    claimers = set()
    try:
        deadline = time.monotonic() + 120
        while any(p.poll() is None for p in procs) and time.monotonic() < deadline:
            for name in os.listdir(q):
                if ".running." in name:
                    claimers.add(int(name.rsplit(".", 1)[1]))
            time.sleep(0.02)
    finally:
        outs = [p.communicate(timeout=60) for p in procs]
    for p, (out, err) in zip(procs, outs):
        assert p.returncode == 0, out + err
    # both workers really worked on the queue at the same time
    assert claimers == {p.pid for p in procs}
    assert _queue(layout) == sorted(f"P25815_{n}.done" for n in names)
    models = _logged_models(log)
    assert len(models) == 40 and len(set(models)) == 40


def test_stale_running_reclaimed(tmp_path):
    layout = _make_run(tmp_path, ["2mjw_B"])
    dead = subprocess.Popen(["true"])
    dead.wait()
    q = layout.pockets / "queue"
    host = socket.gethostname()
    (q / "P25815_2mjw_B.job").rename(q / f"P25815_2mjw_B.running.{host}.{dead.pid}")
    # a claim by a live process on another host is left alone
    (q / "other.running.elsewhere.1").write_text(
        json.dumps({"model": "03_pdbfixer/P25815/2mjw_B.pdb",
                    "out": "04_pockets/P25815/other.pdb", "size": 1}))
    r = subprocess.run(_worker_cmd(tmp_path), env=_env(), capture_output=True, text=True,
                       timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr
    assert _queue(layout) == ["P25815_2mjw_B.done", "other.running.elsewhere.1"]
    assert Path(f"{layout.pockets}/P25815/P25815_2mjw_B.pdb_desc.txt").is_file()


def test_preflight_failure_claims_nothing(tmp_path):
    layout = _make_run(tmp_path, ["2mjw_B", "2mjw_C"])
    bad = tmp_path / "bad_dogsite"
    bad.write_text("#!/bin/sh\necho 'license expired for host' >&2\nexit 1\n")
    bad.chmod(0o755)
    r = subprocess.run(_worker_cmd(tmp_path, bad, "--watch"), env=_env(),
                       capture_output=True, text=True, timeout=60)
    assert r.returncode != 0
    assert "license expired for host" in r.stderr
    assert HINT in r.stderr
    assert _queue(layout) == ["P25815_2mjw_B.job", "P25815_2mjw_C.job"]


def test_worker_missing_binary_is_usage_error(tmp_path):
    _make_run(tmp_path, ["2mjw_B"])
    r = subprocess.run(_worker_cmd(tmp_path, tmp_path / "nope"), env=_env(),
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 2
    assert "nope" in r.stderr


def test_worker_rejects_paths_outside_run_dir(tmp_path):
    layout = _make_run(tmp_path, ["2mjw_B"])
    q = layout.pockets / "queue"
    (q / "evil.job").write_text(json.dumps(
        {"model": "../outside.pdb", "out": "04_pockets/P25815/evil.pdb", "size": 1}))
    r = subprocess.run(_worker_cmd(tmp_path), env=_env(), capture_output=True, text=True,
                       timeout=60)
    assert r.returncode == 0, r.stderr
    assert "outside the run folder" in json.loads((q / "evil.failed").read_text())["stderr_tail"]


def test_watch_stops_at_closed_and_sigterm_unclaims(tmp_path):
    layout = _make_run(tmp_path, ["2mjw_B", "2mjw_C"])
    q = layout.pockets / "queue"
    env = _env(FAKE_DOGSITE_SLEEP="3")
    p = subprocess.Popen(_worker_cmd(tmp_path, FAKE, "--watch"), env=env,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        deadline = time.monotonic() + 30
        while sum(".running." in n for n in os.listdir(q)) < 2:
            assert time.monotonic() < deadline and p.poll() is None
            time.sleep(0.05)
        assert any(n.startswith("WORKER.") for n in os.listdir(q))
        p.terminate()
        p.communicate(timeout=30)
    finally:
        if p.poll() is None:
            p.kill()
    # interrupted jobs go back to the queue; no output that would look done
    assert _queue(layout) == ["P25815_2mjw_B.job", "P25815_2mjw_C.job"]
    assert not list(layout.pockets.glob("P25815/*_desc.txt"))
    out = subprocess.run(["pgrep", "-f", str(FAKE)], capture_output=True, text=True)
    assert out.stdout.strip() == "", "fake dogsite children left running"

    # --watch: runs the queue, then exits only once CLOSED appears
    p = subprocess.Popen(_worker_cmd(tmp_path, FAKE, "--watch"), env=_env(),
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        deadline = time.monotonic() + 30
        while sum(n.endswith(".done") for n in os.listdir(q)) < 2:
            assert time.monotonic() < deadline and p.poll() is None
            time.sleep(0.05)
        time.sleep(0.5)
        assert p.poll() is None              # queue empty, but not closed yet
        (q / "CLOSED").write_text("")
        p.communicate(timeout=30)
    finally:
        if p.poll() is None:
            p.kill()
    assert p.returncode == 0


def test_benchmark_writes_pockets_calibration(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    (home / "calibration.json").write_text(json.dumps(
        {"fix": {"seconds_per_unit": 1e-5, "n_observations": 3, "nproc": 8,
                 "cpu_model": "x", "updated": "2026-01-01T00:00:00+00:00"}}))
    r = subprocess.run([sys.executable, str(WORKER), "--benchmark", "--dogsite", str(FAKE)],
                       env=_env(POLYPHAR_HOME=home, FAKE_DOGSITE_SLEEP="0.2"),
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr
    cal = json.loads((home / "calibration.json").read_text())
    assert cal["fix"]["n_observations"] == 3            # other entries kept
    entry = cal["pockets"]
    assert set(entry) == {"seconds_per_unit", "n_observations", "nproc", "cpu_model", "updated"}
    assert entry["n_observations"] == 1
    # per byte, single process: >= the fake's 0.2 s sleep over the bench bytes
    bench_bytes = sum(p.stat().st_size for p in BENCH.parent.glob("*.pdb"))
    assert 0.2 / bench_bytes <= entry["seconds_per_unit"] < 5.0 / bench_bytes
    assert [p.name for p in home.iterdir()] == ["calibration.json"]   # no temp left
    # the polyphar ETA engine reads it
    from polyphar.eta import Calibration
    assert Calibration.load(home).cost("pockets") == entry["seconds_per_unit"]


def test_worker_is_stdlib_python38_and_standalone():
    src = WORKER.read_text()
    assert src.startswith("#!/usr/bin/env python3\n")
    assert os.access(WORKER, os.X_OK) and os.access(FAKE, os.X_OK)
    tree = ast.parse(src, feature_version=(3, 8))
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            mods.add(node.module.split(".")[0])
    assert mods <= set(sys.stdlib_module_names), mods - set(sys.stdlib_module_names)
    assert "polyphar" not in mods
    for newer in (".removeprefix(", ".removesuffix(", "cancel_futures", "match ", " | None"):
        assert newer not in src, newer
    r = subprocess.run([sys.executable, str(WORKER), "--help"], capture_output=True, text=True)
    assert r.returncode == 0 and "--run-dir" in r.stdout and "--benchmark" in r.stdout


# --- wrapper lifecycle (bin/polyphar-docker) ----------------------------------------

FAKE_DOCKER = """#!/bin/sh
# records whether a DoGSiteScorer worker heartbeat exists while "the container" runs
q="$RUN/04_pockets/queue"
n=$(ls "$q" 2>/dev/null | grep -c '^WORKER\\.' || true)
echo "$n" > "$RUN/docker_saw_workers"
exit "${FAKE_DOCKER_EXIT:-0}"
"""


@pytest.mark.parametrize("argv, worker", [
    (["pockets"], 1),
    (["run"], 1),
    (["run", "--from", "fix", "--to", "pockets"], 1),
    (["run", "--from=pockets"], 1),
    (["run", "--to", "fix"], 0),
    (["run", "--to=structures"], 0),
    (["run", "--from", "similarity"], 0),
    (["pockets", "--dry-run"], 0),
    (["ids"], 0),
])
def test_wrapper_worker_only_when_pockets_in_range(tmp_path, argv, worker):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "docker").write_text(FAKE_DOCKER)
    (bindir / "docker").chmod(0o755)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    home = tmp_path / "home"
    home.mkdir()
    (home / "calibration.json").write_text(json.dumps(
        {"pockets": {"seconds_per_unit": 1e-6, "n_observations": 1, "nproc": 1,
                     "cpu_model": "x", "updated": "x"}}))
    env = _env(PATH=f"{bindir}:{os.environ['PATH']}", RUN=run_dir, POLYPHAR_RUN_DIR=run_dir,
               POLYPHAR_HOME=home, DOGSITE_BIN=FAKE, FAKE_DOCKER_EXIT=7)
    r = subprocess.run([str(WRAPPER), *argv], env=env, capture_output=True, text=True,
                       timeout=60)
    assert r.returncode == 7, r.stdout + r.stderr    # the container's exit status
    assert (run_dir / "docker_saw_workers").read_text().strip() == str(worker)
    left = subprocess.run(["pgrep", "-f", f"polyphar-dogsite-worker --run-dir {run_dir}"],
                          capture_output=True, text=True)
    assert left.stdout.strip() == ""
    q = run_dir / "04_pockets" / "queue"
    assert not q.exists() or not [n for n in os.listdir(q) if n.startswith("WORKER.")]


def test_wrapper_requires_dogsite_bin_for_pockets(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    env = _env(POLYPHAR_RUN_DIR=run_dir, POLYPHAR_HOME=tmp_path / "home",
               PATH="/usr/bin:/bin")
    r = subprocess.run([str(WRAPPER), "pockets"], env=env, capture_output=True, text=True)
    assert r.returncode == 2
    assert "DOGSITE_BIN" in r.stderr


def test_wrapper_benchmarks_pockets_before_run_when_uncalibrated(tmp_path):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "docker").write_text(FAKE_DOCKER)
    (bindir / "docker").chmod(0o755)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    home = tmp_path / "home"
    env = _env(PATH=f"{bindir}:{os.environ['PATH']}", RUN=run_dir, POLYPHAR_RUN_DIR=run_dir,
               POLYPHAR_HOME=home, DOGSITE_BIN=FAKE)
    r = subprocess.run([str(WRAPPER), "run"], env=env, capture_output=True, text=True,
                       timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "pockets" in json.loads((home / "calibration.json").read_text())


@pytest.mark.docker
def test_wrapper_starts_and_stops_worker(tmp_path):
    r = subprocess.run(["docker", "image", "inspect", "polyphar:dev"], capture_output=True)
    if r.returncode != 0:
        pytest.skip("image polyphar:dev not built (see tests/test_docker.py)")
    run_dir = tmp_path / "run"
    layout = RunLayout(run_dir)
    names = ["2mjw_B", "1abc_A"]
    rows = []
    for name in names:
        for base in (layout.structures, layout.pdbfixer):
            p = base / "P25815" / f"{name}.pdb"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(BENCH.read_text())
        pdb, chain = name.split("_")
        rows.append(["P25815", "S100P_HUMAN", pdb.upper(), chain, False, "-", "-", "-",
                     f"02_structures/P25815/{name}.pdb"])
    pd.DataFrame(rows, columns=INFO_COLUMNS).to_csv(layout.protein_info, index=False)
    env = _env(POLYPHAR_RUN_DIR=run_dir, POLYPHAR_VERSION="dev",
               POLYPHAR_HOME=tmp_path / "home", DOGSITE_BIN=FAKE)
    r = subprocess.run([str(WRAPPER), "pockets", "--no-notify"], env=env,
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    df = pd.read_csv(layout.drugscore_csv, index_col=0)
    assert sorted(df["Pocket Filename"]) == sorted(f"P25815_{n}.pdb_res_P_1.pdb" for n in names)
    man = json.loads((layout.pockets / "manifest.json").read_text())
    assert man["versions"]["dogsite_sha256"] == sha256_file(FAKE)
    left = subprocess.run(["pgrep", "-f", f"polyphar-dogsite-worker --run-dir {run_dir}"],
                          capture_output=True, text=True)
    assert left.stdout.strip() == "", "worker still running after the wrapper exited"
    q = layout.pockets / "queue"
    assert not [n for n in os.listdir(q) if n.startswith("WORKER.")]
    assert (q / "CLOSED").exists()
