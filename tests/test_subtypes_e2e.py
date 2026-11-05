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

    DOGSITE_BIN=$PWD/tests/fixtures/fake_dogsite2 pytest -m subtype_e2e

Note the `2`: the fake answers `--help` as DoGSite3 or as DoGSiteScorer 2 depending on the
name it is invoked under, and the default --pocket-tool is dogsite2. Pointing DOGSITE_BIN
at the plain `fake_dogsite` name fails the worker's tool check before any stage runs.

**The docking is the production path: Uni-Dock on a Colab Pro+ GPU.** POSEIDON has no
GPU, so `polyphar-docker` hands `dock` to `bin/polyphar-colab`, which exports a bundle,
installs the pinned commit on the VM, polls, pulls results every poll and merges. Testing
smina instead would exercise an engine the study never runs -- the defects that reached
live data were all in the joins, and the Uni-Dock join is the one in use.

That means these tests need a Colab session and spend a little GPU credit. The 3-ligand
library keeps it small: a unit is seconds rather than the ~6 GPU-minutes a 1,960-ligand
unit costs. POLYPHAR_COLAB_KEEP=1 holds one session across all seven, so the engine
install is paid once instead of seven times.

Prerequisites, all of them hard:
  * the commit under test must be PUSHED -- the VM pip-installs polyphar from it;
  * `colab` must be authenticated on this host;
  * ~/.polyphar/pushover.env must exist, or ALLOW_SILENT=1.

POLYPHAR_HOME is shared across the seven, so only the first pays for the benchmark.
"""
import json
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
    env.update(POLYPHAR_RUN_DIR=str(run_dir), POLYPHAR_HOME=str(home),
               POLYPHAR_COLAB_KEEP="1",      # one session for all seven
               # A unit of 3 ligands finishes in seconds, so the default 120 s poll is
               # the whole cost of the round trip. Production keeps 120 s: there a unit
               # is ~6 GPU-minutes and polling harder would just burn API calls.
               COLAB_POLL_SECONDS=os.environ.get("COLAB_POLL_SECONDS", "15"),
               # POLYPHAR_HOME is a fresh temp dir, so it holds no pushover.env and
               # launch.py would refuse to start a worker that cannot notify. A test does
               # not want notifications; a real run does, which is why the refusal exists.
               ALLOW_SILENT="1",
               # A fresh VM installs micromamba, Uni-Dock and the pinned polyphar in
               # about nine minutes. The production budget is 70, because there a VM you
               # have already waited ten minutes for is worth more than a fast failure.
               # On 2026-11-06 a Colab VM went quiet mid-install and this test sat in that
               # 70-minute budget for 37 minutes before anyone looked; twice the real cost
               # is all the patience a test needs.
               COLAB_SETUP_TIMEOUT=os.environ.get("COLAB_SETUP_TIMEOUT", "1200"),
               POLYPHAR_VERSION=os.environ.get("POLYPHAR_VERSION", "dev"))
    # no --engine: the default is Uni-Dock, and on a GPU-less host that is exactly the
    # production path -- polyphar-docker hands dock to the Colab round trip.
    argv = [str(WRAPPER), "run", "--txt", "inputs/ids.txt",
            "--library", "inputs/bench3.sdf", *extra]
    return subprocess.run(argv, cwd=run_dir, env=env, capture_output=True, text=True,
                          timeout=60 * 60)


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
    # 3 is EXIT_FLAGGED: a stage finished with flagged counts, and for dock that normally
    # means incomplete_targets -- Uni-Dock does not always pose every ligand (issue #133),
    # so a replicate can legitimately come back short and the target waits for a redraw.
    # It is the ordinary state until everything is docked, not a failure. Anything else is.
    assert r.returncode in (0, 3), r.stdout[-4000:] + r.stderr[-4000:]
    run_dir = tmp_path / subtype

    # ... but an incomplete target must be NAMED, never silently absent. This is the
    # habit that let six boron drugs make every target permanently incomplete without a
    # word in any log.
    if r.returncode == 3:
        failed = (run_dir / "07_docking" / "failed.tsv").read_text().splitlines()
        assert len(failed) > 1, f"{subtype}: exit 3 with an empty failed.tsv"

    for d in DIRS:
        assert (run_dir / d / "manifest.json").is_file(), f"{subtype}: {d} has no manifest"

    # the panel is non-empty and every target pairs a pocket with ITS OWN receptor --
    # the invariant that was violated for 83 of net1's 479 rows
    targets = (run_dir / "06_targets" / "final_targets.csv").read_text().splitlines()
    assert len(targets) >= 2, f"{subtype}: no target rows"

    # docking produced a summary per complete target, over the ligands that survived
    results = list((run_dir / "07_docking").glob("df_results_*.csv"))
    assert results, f"{subtype}: no df_results"   # at least one target came out whole
    assert all(len(p.read_text().splitlines()) > 1 for p in results)

    # Uni-Dock on a GPU, not a local fallback. A run that quietly fell back to CPU would
    # still produce results, so the provenance is what tells the production path from a
    # lookalike. (Checked here rather than in its own test: tmp_path is per-test, so a
    # second test cannot see this run folder at all -- it would skip, always.)
    provs = list((run_dir / "07_docking").glob("*/replicates/*.prov.json"))
    assert provs, f"{subtype}: no unit provenance"
    for prov in provs:
        info = json.loads(prov.read_text())
        assert info["engine_name"] == "unidock", info
        assert info.get("gpu_model"), f"{subtype}: no GPU recorded in {prov.name}"
