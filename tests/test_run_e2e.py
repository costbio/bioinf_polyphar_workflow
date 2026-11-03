"""End to end: `polyphar run` on two proteins through bin/polyphar-docker.

Needs internet (UniProt, RCSB, AlphaFold DB), the polyphar image (see tests/test_docker.py for
the build command; tag from $POLYPHAR_VERSION, default dev) and $DOGSITE_BIN (the host
DoGSiteScorer, or tests/fixtures/fake_dogsite for a local run: Ruling 5). PocketMatch and
smina run inside the image. The 3-ligand bench library (typical small molecules) keeps the
docking stage short. POLYPHAR_HOME is a fresh temp dir, so every stage is uncalibrated and
`run` benchmarks it first (pockets: the wrapper's worker --benchmark).

$POLYPHAR_E2E_PMAX overrides select's --pmax (fake pockets may leave no pair >= 0.7).
"""
import os
import shutil
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
WRAPPER = REPO / "bin" / "polyphar-docker"
LIBRARY = REPO / "polyphar" / "data" / "bench" / "edrug3d_bench3.sdf"
# benchmark.ORDER, paired with DIRS below (benchmark.STAGE_DIRS).
ORDER = ("ids", "uniprot", "structures", "fix", "pockets", "annotate", "select",
         "similarity", "ligands", "dock")
# File mtimes come from the kernel's coarse clock (one tick, <= 10 ms behind the clock that
# stamped the notification), so a file written right after the notification can carry an
# mtime a few ms before it.
MTIME_SLACK = timedelta(milliseconds=20)
DIRS = ("00_ids", "01_uniprot", "02_structures", "03_pdbfixer", "04_pockets",
        "04_pockets/annotation", "06_targets", "05_similarity", "07_docking/library",
        "07_docking")
# Subdirectories a stage dir contains that belong to ANOTHER stage (or to a queue), so the
# "no output before the ETA" check does not see the earlier stage's files.
NESTED = {"04_pockets": ("queue", "annotation"), "07_docking": ("library",)}


def _notifications(run_dir):
    """[(timestamp, priority, title, message)] from logs/notifications.log, in order."""
    out = []
    for line in (run_dir / "logs" / "notifications.log").read_text().splitlines():
        stamp, prio, _sent, title, msg = line.split("\t")
        out.append((datetime.fromisoformat(stamp), int(prio), title, msg))
    return out


def _first_output_mtime(stage_dir: Path, skip=()):
    times = [p.stat().st_mtime for p in stage_dir.rglob("*")
             if p.is_file() and not any(part in skip for part in p.relative_to(stage_dir).parts)]
    return datetime.fromtimestamp(min(times), tz=timezone.utc) if times else None


@pytest.mark.network
@pytest.mark.dogsite
@pytest.mark.docker
def test_e2e_two_proteins(tmp_path):
    version = os.environ.get("POLYPHAR_VERSION", "dev")
    if subprocess.run(["docker", "image", "inspect", f"polyphar:{version}"],
                      capture_output=True).returncode != 0:
        pytest.skip(f"image polyphar:{version} not built")
    run_dir = tmp_path / "e2e_run"
    (run_dir / "inputs").mkdir(parents=True)
    (run_dir / "inputs" / "ids.txt").write_text("P25815\nP53350\n")
    shutil.copyfile(LIBRARY, run_dir / "inputs" / "bench3.sdf")
    env = {k: v for k, v in os.environ.items() if not k.startswith("PUSHOVER_")}
    env.update(POLYPHAR_RUN_DIR=str(run_dir), POLYPHAR_VERSION=version,
               POLYPHAR_HOME=str(tmp_path / "home"))
    argv = [str(WRAPPER), "run", "--txt", "inputs/ids.txt", "--library", "inputs/bench3.sdf"]
    if os.environ.get("POLYPHAR_E2E_PMAX"):
        argv += ["--pmax", os.environ["POLYPHAR_E2E_PMAX"]]
    r = subprocess.run(argv, cwd=run_dir, env=env, capture_output=True, text=True,
                       timeout=4 * 3600)
    (tmp_path / "run.out").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]

    for d in DIRS:
        assert (run_dir / d / "manifest.json").is_file(), f"{d} has no manifest"
    targets = (run_dir / "06_targets" / "final_targets.csv").read_text().splitlines()
    assert len(targets) >= 2, "final_targets.csv has no target rows"
    results = list((run_dir / "07_docking").glob("df_results_*_fixed.csv"))
    assert results and all(len(p.read_text().splitlines()) > 1 for p in results)

    # one start/ETA notification per stage, before that stage's first unit
    notes = _notifications(run_dir)
    for stage, d in zip(ORDER, DIRS):
        starts = [n for n in notes if n[2] == f"{stage} starting"]
        assert len(starts) == 1, f"{stage}: {len(starts)} start notifications"
        t_start, _, _, msg = starts[0]
        assert msg.startswith(f"{stage}: ") and " - ETA " in msg, msg
        later = [n for n in notes if n[2].startswith(f"{stage} ") and n is not starts[0]]
        assert later and all(n[0] > t_start for n in later), f"{stage}: notification order"
        first_out = _first_output_mtime(run_dir / d, skip=NESTED.get(d, ("queue",)))
        assert first_out is not None and first_out >= t_start - MTIME_SLACK, \
            f"{stage}: output at {first_out} before the ETA at {t_start}"
        assert f"{stage}: " in r.stdout
    title = notes[-1][2]
    assert title.startswith("run done"), title
    assert "run summary" in r.stdout and "run summary" in notes[-1][3]
    print(r.stdout[r.stdout.index("run summary"):])
