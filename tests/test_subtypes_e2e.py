"""The whole pipeline, once per disease subtype, on that subtype's own accessions.

Every stage had tests before this; the chain between them had one, on two proteins. The
three defects that reached live data in November 2026 all lived in the joins -- a receptor
paired with another structure's pocket, a ligand the engine silently never posed, a pocket
PocketMatch could not type -- and each was trivial once seen. A per-stage test cannot see
any of them.

Runs `ids -> uniprot -> structures -> fix -> pockets -> annotate -> select -> similarity
-> ligands -> dock` through bin/polyphar-docker, over the real UniProt, RCSB, SIFTS and
AlphaFold, with fake DoGSiteScorer, the 3-ligand bench library and smina. Two accessions
per subtype, taken from that subtype's real input list.

    DOGSITE_BIN=$PWD/tests/fixtures/fake_dogsite pytest -m subtype_e2e

`--engine smina` is not a detail: on a GPU-less host the default (Uni-Dock) makes
polyphar-docker hand the docking to a REAL Colab session, which a test must never do.

POLYPHAR_HOME is shared across the seven, so only the first pays for the benchmark.
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
WRAPPER = REPO / "bin" / "polyphar-docker"
LIBRARY = REPO / "polyphar" / "data" / "bench" / "edrug3d_bench3.sdf"

# Every stage dir that must carry a manifest when the run is done.
DIRS = ("00_ids", "01_uniprot", "02_structures", "03_pdbfixer", "04_pockets",
        "04_pockets/annotation", "06_targets", "05_similarity", "07_docking/library",
        "07_docking")

# Two reviewed accessions from each subtype's own 00_ids/ids.txt on POSEIDON.
SUBTYPES = {
    "net1": ("P02814", "O75638"),
    "net2": ("P43355", "O14746"),
    "net3": ("O60829", "P01037"),
    "Cluster1": ("P68366", "Q8TBB1"),
    "Cluster2": ("O60356", "Q13323"),
    "Cluster3": ("O60885", "Q99750"),
    "Cluster4": ("Q96RK0", "O94761"),
}


@pytest.fixture(scope="module")
def home(tmp_path_factory):
    """One calibration directory for all seven, so only the first run benchmarks."""
    return tmp_path_factory.mktemp("polyphar_home")


def _run(run_dir, home, accessions, extra=()):
    (run_dir / "inputs").mkdir(parents=True)
    (run_dir / "inputs" / "ids.txt").write_text("".join(f"{a}\n" for a in accessions))
    shutil.copyfile(LIBRARY, run_dir / "inputs" / "bench3.sdf")
    env = {k: v for k, v in os.environ.items() if not k.startswith("PUSHOVER_")}
    # no SMINA_BIN: the image ships the real smina, and a host path would not resolve
    # inside the container. The 3-ligand library is what keeps the docking short.
    env.update(POLYPHAR_RUN_DIR=str(run_dir), POLYPHAR_HOME=str(home),
               POLYPHAR_VERSION=os.environ.get("POLYPHAR_VERSION", "dev"))
    argv = [str(WRAPPER), "run", "--txt", "inputs/ids.txt",
            "--library", "inputs/bench3.sdf", "--engine", "smina", *extra]
    return subprocess.run(argv, cwd=run_dir, env=env, capture_output=True, text=True,
                          timeout=30 * 60)


def _image_or_skip():
    version = os.environ.get("POLYPHAR_VERSION", "dev")
    if subprocess.run(["docker", "image", "inspect", f"polyphar:{version}"],
                      capture_output=True).returncode != 0:
        pytest.skip(f"image polyphar:{version} not built")


@pytest.mark.subtype_e2e
@pytest.mark.network
@pytest.mark.dogsite
@pytest.mark.docker
@pytest.mark.parametrize("subtype", sorted(SUBTYPES))
def test_a_subtype_runs_end_to_end(subtype, tmp_path, home):
    _image_or_skip()
    accessions = SUBTYPES[subtype]
    r = _run(tmp_path / subtype, home, accessions)
    assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-4000:]
    run_dir = tmp_path / subtype

    for d in DIRS:
        assert (run_dir / d / "manifest.json").is_file(), f"{subtype}: {d} has no manifest"

    # the panel is non-empty and every target pairs a pocket with ITS OWN receptor --
    # the invariant that was violated for 83 of net1's 479 rows
    targets = (run_dir / "06_targets" / "final_targets.csv").read_text().splitlines()
    assert len(targets) >= 2, f"{subtype}: no target rows"

    # docking produced a summary per complete target, over the ligands that survived
    results = list((run_dir / "07_docking").glob("df_results_*.csv"))
    assert results, f"{subtype}: no df_results"
    assert all(len(p.read_text().splitlines()) > 1 for p in results)


@pytest.mark.subtype_e2e
@pytest.mark.network
@pytest.mark.dogsite
@pytest.mark.docker
def test_a_run_never_hands_docking_to_colab_without_a_gpu(tmp_path, home):
    """A test that reached the real Colab path would spend GPU credit and need a session.
    `--engine smina` must keep the docking local; this pins that the wrapper agrees."""
    _image_or_skip()
    r = _run(tmp_path / "local", home, SUBTYPES["net1"], extra=("--dry-run",))
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert "no local GPU" not in r.stdout, r.stdout[-2000:]
    assert "colab" not in r.stdout.lower(), r.stdout[-2000:]
