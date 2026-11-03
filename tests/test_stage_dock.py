import json
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

openbabel = pytest.importorskip(
    "openbabel",
    reason="dock needs openbabel for receptor PDBQT preparation; without it "
           "these tests FAILED rather than skipped, which hid real regressions "
           "behind a missing dependency")

from polyphar.stages import dock

FIX = Path(__file__).parent / "fixtures"
D2MJW = FIX / "dock_2mjw"
RECEPTOR = D2MJW / "2mjw_A.pdb"
POCKET = D2MJW / "P25815_2mjw_A.pdb_res_P_1.pdb"
HEAD20 = D2MJW / "edrug3d_head20.sdf"
LEGACY = D2MJW / "legacy_pose1_head20.csv"
LEGACY_PDBQT = D2MJW / "2mjw_A_legacy2023.pdbqt"   # the 2023 receptor PDBQT
AF_MODEL = FIX / "dock_af" / "AF-P01350-F1.pdb"
FAKE_SMINA = FIX / "fake_smina"
LEGACY_SEED = -496630900   # "Using random seed" in the legacy screen_2mjw_eDrug.log


def _subset(src: Path, ids, dst: Path) -> Path:
    """The library records with these IDs, in the given order, titled by ID."""
    recs = {str(r.fields["ID"]): r for r in dock._records(src.read_text())}
    dst.write_text("".join(dock._titled(recs[str(i)]) for i in ids))
    return dst


# --- legacy equivalence and determinism (real smina) ---------------------------------

@pytest.mark.docking
def test_reproduces_legacy_2mjw_scores(tmp_path, record_property):
    """Agreement, not identity, with the 2023 screen: the 2023 receptor PDBQT and the
    logged seed with today's smina. (Bit-exact reproduction is impossible: the 2023 PDBQT
    came from a random hydrogen placement and an older smina build.) Fails only if the
    Spearman rho of the pose-1 scores drops below 0.7."""
    lig = tmp_path / "lig.sdf"
    dock.write_titled_library(HEAD20, lig)
    out = tmp_path / "out.sdf"
    center, size = dock.pocket_box(POCKET)
    dock.run_smina(LEGACY_PDBQT, lig, center, size, out, tmp_path / "out.log",
                   seed=LEGACY_SEED, cpu=12, exhaustiveness=12, num_modes=12)
    got = dock.parse_poses(out, dock.read_library(HEAD20), "2mjw_P_1")
    got = got[got.Pose == 1].reset_index(drop=True)
    legacy = pd.read_csv(LEGACY, index_col=0).reset_index(drop=True)
    assert got.Ligand.tolist() == legacy.Ligand.tolist()
    rho = got.SMINA_Score.corr(legacy.SMINA_Score, method="spearman")
    mad = (got.SMINA_Score.round(1) - legacy.SMINA_Score).abs().mean()   # log has 1 decimal
    same = int((got.SMINA_Score.round(1) == legacy.SMINA_Score).sum())
    record_property("spearman_rho", round(float(rho), 4))
    record_property("mean_abs_delta", round(float(mad), 4))
    print(f"\nlegacy 2mjw agreement: spearman_rho={rho:.4f} mean_abs_delta={mad:.3f} "
          f"identical={same}/20")
    assert rho >= 0.7


# --- receptor preparation ---------------------------------------------------------------

def _atoms(path):
    return [l for l in path.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]


def _polar_h(pdb: Path) -> set:
    """Coordinates of the hydrogens bonded (< 1.25 A) to N, O or S."""
    rows = [(l[76:78].strip(), (float(l[30:38]), float(l[38:46]), float(l[46:54])))
            for l in _atoms(pdb)]
    heavy = [x for e, x in rows if e in ("N", "O", "S")]
    out = set()
    for e, x in rows:
        if e == "H" and min(sum((a - b) ** 2 for a, b in zip(x, y)) for y in heavy) < 1.25 ** 2:
            out.add(tuple(round(c, 3) for c in x))
    return out


def test_receptor_pdbqt_deterministic_and_keeps_pdbfixer_hydrogens(tmp_path):
    a, b = tmp_path / "a" / "r.pdbqt", tmp_path / "b" / "r.pdbqt"
    a.parent.mkdir(), b.parent.mkdir()
    dock.receptor_to_pdbqt(RECEPTOR, a)
    dock.receptor_to_pdbqt(RECEPTOR, b)
    assert a.read_bytes() == b.read_bytes()   # also independent of the output folder
    lines = _atoms(a)
    assert not any(l.startswith(("ROOT", "BRANCH")) for l in a.read_text().splitlines())
    hd = {tuple(round(float(l[k:k + 8]), 3) for k in (30, 38, 46))
          for l in lines if l.split()[-1] == "HD"}
    # no hydrogen added, removed or moved: every polar H of the pdbfixer receptor is an HD
    # at its input coordinates (non-polar H are merged into their carbon, as AutoDock does)
    assert hd == _polar_h(RECEPTOR)
    heavy_in = sum(1 for l in _atoms(RECEPTOR) if l[76:78].strip() != "H")
    assert len(lines) - len(hd) == heavy_in
    assert "REMARK  Name = 2mjw_A.pdb" in a.read_text()


def test_alphafold_receptor_gets_hydrogens_deterministically(tmp_path):
    assert not any(l[76:78].strip() == "H" for l in _atoms(AF_MODEL))
    out = []
    for k in (1, 2):
        h, q = tmp_path / f"h{k}.pdb", tmp_path / f"h{k}.pdbqt"
        dock.protonate_alphafold(AF_MODEL, h)
        dock.receptor_to_pdbqt(h, q)
        out.append((h.read_bytes(), q.read_text().split("\n", 1)[1]))
    assert out[0] == out[1]
    h_atoms = [l for l in _atoms(tmp_path / "h1.pdb") if l[76:78].strip() == "H"]
    heavy = [l for l in _atoms(tmp_path / "h1.pdb") if l[76:78].strip() != "H"]
    assert len(h_atoms) > 0 and len(heavy) == len(_atoms(AF_MODEL))   # hydrogens only
    assert b"CREATED WITH OPENMM" not in out[0][0]


# --- library -----------------------------------------------------------------------------

def test_library_ids_unique_and_complete(tmp_path):
    lib = dock.read_library(HEAD20)
    assert list(lib.columns) == ["ID", "name", "cas", "status", "record_index"]
    assert lib.ID.tolist() == [str(i) for i in range(1, 21)]
    assert lib.record_index.tolist() == list(range(20))
    assert lib.name[3] == "BUTABARBITAL"   # value stripped
    recs = dock._records(HEAD20.read_text())
    dup = tmp_path / "dup.sdf"
    dup.write_text("".join("".join(r.lines) for r in recs[:2] + recs[:1]))
    with pytest.raises(ValueError, match="duplicate ID"):
        dock.read_library(dup)
    noid = tmp_path / "noid.sdf"
    noid.write_text("".join(recs[0].lines).replace("> <ID>", "> <IDX>"))
    with pytest.raises(ValueError, match="without an ID"):
        dock.read_library(noid)
    titled = tmp_path / "t.sdf"
    dock.write_titled_library(HEAD20, titled)
    t = dock._records(titled.read_text())
    assert [r.lines[0].strip() for r in t] == [f"EDRUG3D_{i}" for i in range(1, 21)]
    assert all(a.lines[1:] == b.lines[1:] for a, b in zip(t, recs))   # rest verbatim


def test_placeholder_titles_resolved_by_id(tmp_path):
    recs = dock._records(HEAD20.read_text())
    src = tmp_path / "mol1.sdf"
    src.write_text("".join("Molecule-1\n" + "".join(r.lines[1:]) for r in recs[2:4]))
    lib = dock.read_library(src)
    assert lib.ID.tolist() == ["3", "4"] and lib.name.tolist() == ["HISTAMINE", "BUTABARBITAL"]
    titled = tmp_path / "t.sdf"
    dock.write_titled_library(src, titled)
    out = tmp_path / "out.sdf"
    subprocess.run([sys.executable, str(FAKE_SMINA), "-l", str(titled), "--out", str(out),
                    "--log", str(tmp_path / "l"), "--seed", "1", "--cpu", "1"], check=True)
    df = dock.parse_poses(out, lib, "X_P_1")
    assert df.columns.tolist() == ["Pose", "SMINA_Score", "Ligand", "ID", "Pocket"]
    assert df.groupby("ID").Ligand.first().to_dict() == {"3": "HISTAMINE", "4": "BUTABARBITAL"}
    assert df.Pose.tolist() == [1, 2, 1, 2]


# --- the stage (fake smina, real receptor preparation) -----------------------------------------
class Prog:
    def __init__(self):
        self.calls = []

    def advance(self, n=1, size=None):
        self.calls.append(n)


def _args(**kw):
    base = dict(nprocs=2, library="inputs/lib.sdf", seed=42, replicates=3,
                threads_per_job=1, dry_run=False)
    base.update(kw)
    return SimpleNamespace(**base)


@pytest.fixture
def run_folder(tmp_path, monkeypatch):
    """A run folder with one PDB target (2mjw) and one AlphaFold target, a 20-ligand
    library, and the fake smina."""
    from polyphar.common import RunLayout

    monkeypatch.setenv("SMINA_BIN", str(FAKE_SMINA))
    calls = tmp_path / "calls.txt"
    monkeypatch.setenv("FAKE_SMINA_CALLS", str(calls))
    monkeypatch.delenv("FAKE_SMINA_SKIP_IDS", raising=False)
    monkeypatch.delenv("FAKE_SMINA_JITTER", raising=False)
    run = tmp_path / "run"
    layout = RunLayout(run)
    files = {"03_pdbfixer/P25815/2mjw_A.pdb": RECEPTOR,
             "04_pockets/P25815/P25815_2mjw_A.pdb_res_P_1.pdb": POCKET,
             "02_structures/P01350/AF-P01350-F1.pdb": AF_MODEL,
             "04_pockets/P01350/AF-P01350-F1.pdb_res_P_1.pdb": POCKET,
             "inputs/lib.sdf": HEAD20}
    for rel, src in files.items():
        (run / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, run / rel)
    layout.targets.mkdir(parents=True)
    pd.DataFrame({
        "Node": ["n1", "n2"], "UniProt_ID": ["P25815", "P01350"], "Gene_Name": ["S100P", "GAST"],
        "Organism": ["Homo sapiens"] * 2, "Structure_ID": ["2mjw", "AF-P01350-F1"],
        "Similarity": [1.0, 0.8], "Protein": ["2mjw_A", "AF-P01350-F1"], "Chain": ["A", "A"],
        "Pocket": ["P_1", "P_1"],
        "Pocket_Path": ["04_pockets/P25815/P25815_2mjw_A.pdb_res_P_1.pdb",
                        "04_pockets/P01350/AF-P01350-F1.pdb_res_P_1.pdb"],
        "Receptor_Path": ["03_pdbfixer/P25815/2mjw_A.pdb", "02_structures/P01350/AF-P01350-F1.pdb"],
        "Is_Hub": [True, False]}).to_csv(layout.final_targets_csv, index=False)
    return layout, calls


LABELS = ("2mjw_A_P_1", "AF-P01350-F1_P_1")


def _plan_of(layout, args):
    lib = dock._library_path(layout, args, copy=False)
    return dock._plan(layout, args, lib, dock.smina_version())


def _discards(layout):
    path = layout.docking / dock.DISCARD_LOG
    if not path.is_file():
        return pd.DataFrame(columns=["label", "replicate", "reason"])
    return pd.read_csv(path, sep="\t", dtype=str)


# --- the unit is (target, replicate) ---------------------------------------------------

def test_two_replicates_of_one_target_are_distinct_units(run_folder):
    """The replicate index and its seed are part of a unit's identity, so replicate 1 is
    not a re-run of replicate 0: different key, different seed, different result file."""
    layout, _ = run_folder
    plan = _plan_of(layout, _args(replicates=3))
    mine = sorted((u for u in plan.units if u.label == LABELS[0]), key=lambda u: u.replicate)
    assert [u.replicate for u in mine] == [0, 1, 2]
    assert [u.seed for u in mine] == [42, 43, 44]          # base_seed + i, recorded
    assert len({u.key for u in mine}) == 3                  # the index is IN the key
    assert len({dock.unit_paths(layout, u)[0] for u in mine}) == 3

    # the only thing that differs between them is the replicate index and its seed
    same = [u for u in _plan_of(layout, _args(replicates=3)).units if u.label == LABELS[0]]
    assert {u.key for u in same} == {u.key for u in mine}   # planning twice is stable
    # ... and raising --replicates adds units without touching the existing keys
    more = _plan_of(layout, _args(replicates=5)).units
    assert len(more) == 2 * 5
    assert {u.key for u in mine} <= {u.key for u in more}


def test_run_writes_one_result_per_replicate(run_folder):
    layout, calls = run_folder
    args = _args()
    assert dock.units(layout, args).count == 120   # 2 targets x 3 replicates x 20 ligands
    prog = Prog()
    counts = dock.run(layout, args, prog)
    assert sum(prog.calls) == 120
    assert counts["units"] == 6 and counts["units_docked"] == 6 and counts["failed"] == 0
    assert counts["replicates"] == 3 and counts["complete"] == 2
    assert counts["incomplete_targets"] == 0 and counts["discarded"] == 0
    for label in LABELS:
        for i in range(3):
            assert (layout.docking / label / "replicates" / f"r{i}.sdf").is_file()
            assert (layout.docking / label / "replicates" / f"r{i}.done").is_file()
        df = pd.read_csv(layout.docking / f"df_results_{label}_fixed.csv", dtype={"ID": str})
        assert df.ID.tolist() == [str(i) for i in range(1, 21) for _ in (1, 2)]
        assert (layout.docking / label / "poses.sdf").is_file()
    # one engine call per unit, over the WHOLE library, with that replicate's seed
    lines = calls.read_text().splitlines()
    assert len(lines) == 6
    assert {l.split()[2] for l in lines} == {"n=20"}
    assert sorted(l.split()[4] for l in lines) == ["seed=42"] * 2 + ["seed=43"] * 2 + \
        ["seed=44"] * 2
    assert {l.split()[0] for l in lines} == {
        str(layout.docking / dock.LIB_DIRNAME / dock.TITLED_LIBRARY)}
    m = json.loads((layout.docking / "manifest.json").read_text())
    assert m["base_seed"] == 42 and m["engine"] == "smina fake 0.0"
    assert m["args"]["replicates"] == 3 and m["args"]["unit_timeout"] == 43200
    assert len(m["units"]) == 6
    assert {(u["replicate"], u["seed"]) for u in m["units"]} == {(0, 42), (1, 43), (2, 44)}
    assert {u["state"] for u in m["units"]} == {"done"}
    assert (layout.docking / "failed.tsv").read_text() == "item\tstage\treason\n"
    assert not (layout.docking / dock.DISCARD_LOG).exists()   # nothing was thrown away


# --- resume: skip, never redo ----------------------------------------------------------

def test_resume_skips_a_completed_replicate_and_never_reruns_it(run_folder):
    """A .done marker certifies that the draw was made and recorded, NOT what it scored.
    Re-running it would produce a different sample, so resume must leave it alone."""
    layout, calls = run_folder
    dock.run(layout, _args(), Prog())
    (layout.docking / "manifest.json").unlink()          # as if killed before the end
    # one unit loses its marker entirely: that one, and only that one, is outstanding
    (layout.docking / LABELS[0] / "replicates" / "r1.done").unlink()
    before = {label: {i: (layout.docking / label / "replicates" / f"r{i}.sdf").read_bytes()
                      for i in range(3)} for label in LABELS}
    assert dock.units(layout, _args()).count == 20       # one unit, whole library

    calls.write_text("")
    prog = Prog()
    counts = dock.run(layout, _args(), prog)
    assert prog.calls == [20]
    assert len(calls.read_text().splitlines()) == 1      # exactly one engine call
    assert counts["units_docked"] == 1 and counts["skipped_existing"] == 5
    # the lost marker left an UNCERTIFIED result file behind; replacing it is itself a
    # discarded draw, so it is logged rather than quietly overwritten
    assert counts["discarded"] == 1 and counts["complete"] == 2
    assert "uncertified" in _discards(layout).reason[0]
    # every OTHER replicate's record is byte-identical: untouched, not re-drawn
    after = {label: {i: (layout.docking / label / "replicates" / f"r{i}.sdf").read_bytes()
                     for i in range(3)} for label in LABELS}
    assert after[LABELS[1]] == before[LABELS[1]]
    assert {i: after[LABELS[0]][i] for i in (0, 2)} == {i: before[LABELS[0]][i] for i in (0, 2)}

    # finished again: a second pass runs nothing at all
    calls.write_text("")
    assert dock.units(layout, _args()).count == 0
    dock.run(layout, _args(), Prog())
    assert calls.read_text() == ""


def test_a_recorded_replicate_is_not_redrawn_even_with_a_jittering_engine(run_folder,
                                                                         monkeypatch):
    """With a non-deterministic engine, re-running a unit would silently replace one
    sample with another. Resume must make that impossible."""
    layout, _ = run_folder
    monkeypatch.setenv("FAKE_SMINA_JITTER", "1")
    dock.run(layout, _args(), Prog())
    scores = {i: (layout.docking / LABELS[0] / "replicates" / f"r{i}.sdf").read_text()
              for i in range(3)}
    assert len({v for v in scores.values()}) == 3        # the replicates ARE different draws
    (layout.docking / "manifest.json").unlink()
    dock.run(layout, _args(), Prog())
    assert {i: (layout.docking / LABELS[0] / "replicates" / f"r{i}.sdf").read_text()
            for i in range(3)} == scores


def test_a_fixed_valued_engine_gives_identical_replicates_for_any_nprocs(run_folder):
    """The counterpart: with a stubbed engine that returns fixed values the records are
    exact and independent of nprocs. Nothing here asserts a real engine's scores."""
    layout, _ = run_folder
    dock.run(layout, _args(nprocs=1), Prog())
    a = [(layout.docking / LABELS[0] / "replicates" / f"r{i}.sdf").read_bytes()
         for i in range(3)]
    assert a[0] == a[1] == a[2]
    shutil.rmtree(layout.docking)
    dock.run(layout, _args(nprocs=4), Prog())
    assert [(layout.docking / LABELS[0] / "replicates" / f"r{i}.sdf").read_bytes()
            for i in range(3)] == a


# --- discard and re-draw ---------------------------------------------------------------

def test_a_corrupted_unit_is_discarded_redrawn_as_the_same_replicate_and_logged(run_folder,
                                                                                calls=None):
    """The one legitimate repeat. The discarded draw was a real sample, so dropping it is
    a (small) selection effect and has to leave a record."""
    layout, calls = run_folder
    dock.run(layout, _args(), Prog())
    (layout.docking / "manifest.json").unlink()
    victim = layout.docking / LABELS[0] / "replicates" / "r1.sdf"
    marker = layout.docking / LABELS[0] / "replicates" / "r1.done"
    key_before = marker.read_text()
    victim.write_text("EDRUG3D_1\n truncated, no terminator\n")   # a half-copied record

    plan = _plan_of(layout, _args())
    u = next(u for u in plan.units if u.label == LABELS[0] and u.replicate == 1)
    assert dock.unit_state(layout, u) == dock.STATE_CORRUPT
    assert dock.units(layout, _args()).count == 20

    calls.write_text("")
    counts = dock.run(layout, _args(), Prog())
    assert counts["discarded"] == 1 and counts["units_docked"] == 1
    assert counts["skipped_existing"] == 5 and counts["complete"] == 2
    # re-drawn under the SAME replicate index and the same seed, not as a new replicate
    assert [l.split()[4] for l in calls.read_text().splitlines()] == ["seed=43"]
    assert sorted(p.name for p in (layout.docking / LABELS[0] / "replicates").glob("r*.sdf")) \
        == ["r0.sdf", "r1.sdf", "r2.sdf"]
    assert marker.read_text() == key_before          # same identity, a different sample
    assert dock.unit_state(layout, u) == dock.STATE_DONE

    log = _discards(layout)
    assert len(log) == 1
    assert log.label[0] == LABELS[0] and log.replicate[0] == "1" and log.seed[0] == "43"
    assert "re-drawn under the same replicate index" in log.reason[0]
    assert log.unit_key[0] == key_before.strip()


def test_a_missing_result_file_is_a_corrupt_unit_not_a_done_one(run_folder):
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())
    (layout.docking / "manifest.json").unlink()
    (layout.docking / LABELS[1] / "replicates" / "r0.sdf").unlink()   # marker survives
    plan = _plan_of(layout, _args())
    u = next(u for u in plan.units if u.label == LABELS[1] and u.replicate == 0)
    assert dock.unit_state(layout, u) == dock.STATE_CORRUPT
    counts = dock.run(layout, _args(), Prog())
    assert counts["discarded"] == 1 and counts["complete"] == 2


def test_the_discard_log_is_append_only_across_runs(run_folder):
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())
    for rep in (0, 1):
        (layout.docking / "manifest.json").unlink()
        (layout.docking / LABELS[0] / "replicates" / f"r{rep}.sdf").unlink()
        dock.run(layout, _args(), Prog())
    log = _discards(layout)
    assert log.replicate.tolist() == ["0", "1"]          # the first row survived the second run
    assert (layout.docking / "failed.tsv").read_text() == "item\tstage\treason\n"   # that one IS reset


def test_a_superseded_record_is_logged_before_it_is_overwritten(run_folder):
    """Changing the base seed makes every unit a different draw. The old records are
    replaced -- which is still a discarded sample, so it is still logged."""
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())
    counts = dock.run(layout, _args(seed=100), Prog())
    assert counts["discarded"] == 6 and counts["units_docked"] == 6
    log = _discards(layout)
    assert len(log) == 6 and set(log.reason.str.startswith("superseded")) == {True}


# --- orphans ----------------------------------------------------------------------------

def test_orphans_never_sweep_replicate_results(run_folder):
    """The hazard the replicate layout introduces: _orphans decides on whole top-level
    labels, so a current target's replicate records must be invisible to it however many
    of them there are."""
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())
    plan = _plan_of(layout, _args())
    assert dock._orphans(layout, plan) == []
    dock._obsolete(layout, dock._orphans(layout, plan))
    for label in LABELS:
        assert sorted(p.name for p in (layout.docking / label / "replicates").glob("*.sdf")) \
            == ["r0.sdf", "r1.sdf", "r2.sdf"]

    # a replicate count that no longer matches the plan is still not an orphan
    plan5 = _plan_of(layout, _args(replicates=5))
    assert dock._orphans(layout, plan5) == []
    plan1 = _plan_of(layout, _args(replicates=1))
    assert dock._orphans(layout, plan1) == []

    # ... but a whole label `select` dropped still is, and it is MOVED, not deleted
    ft = pd.read_csv(layout.final_targets_csv)
    ft[ft.Protein == "2mjw_A"].to_csv(layout.final_targets_csv, index=False)
    counts = dock.run(layout, _args(), Prog())
    assert counts["obsoleted"] == 2 and counts["complete"] == 1
    assert not (layout.docking / LABELS[1]).exists()
    moved = {o["to"] for o in json.loads(
        (layout.docking / "manifest.json").read_text())["obsoleted"]}
    assert any((layout.docking / to / "replicates" / "r2.sdf").is_file()
               for to in moved if (layout.docking / to).is_dir())
    # the surviving target kept every draw
    assert sorted(p.name for p in (layout.docking / LABELS[0] / "replicates").glob("*.sdf")) \
        == ["r0.sdf", "r1.sdf", "r2.sdf"]


def test_orphans_ignores_nested_directories_entirely(tmp_path):
    """Structural counterpart: whatever a current label contains, _orphans looks only at
    07_docking's own entries."""
    import types

    from polyphar.common import RunLayout

    layout = RunLayout(tmp_path / "run")
    d = layout.docking
    for sub in ("T_P_1/replicates", "T_P_1/replicates", "GONE_P_9/replicates",
                "ligands", "obsolete", "library", "export", "discarded"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    (d / "T_P_1" / "replicates" / "r7.sdf").write_text("$$$$\n")
    plan = types.SimpleNamespace(targets=[types.SimpleNamespace(label="T_P_1")])
    assert sorted(p.name for p in dock._orphans(layout, plan)) == ["GONE_P_9"]


def test_orphans_never_sweep_the_export_or_discard_directories(tmp_path):
    """`dock_export`'s default bundle directory is 07_docking/export/<ts>/ and the
    archive of discarded draws is 07_docking/discarded/. Both are inside the stage
    directory and neither is a target label, so before they were reserved the orphan
    sweep -- which runs on every `dock` AND every `dock-merge` -- relocated a
    worked-but-unmerged bundle, and the archived draws, into obsolete/."""
    import types

    from polyphar.common import RunLayout
    from polyphar.stages import dock_export

    layout = RunLayout(tmp_path / "run")
    d = layout.docking
    for sub in (f"{dock_export.EXPORT_DIRNAME}/20260101T000000", "discarded/20260101T000000"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    (d / dock_export.EXPORT_DIRNAME / "20260101T000000" / "units.tsv").write_text("x\n")
    plan = types.SimpleNamespace(targets=[])
    assert dock._orphans(layout, plan) == []
    dock._obsolete(layout, dock._orphans(layout, plan))
    assert (d / dock_export.EXPORT_DIRNAME / "20260101T000000" / "units.tsv").is_file()
    assert (d / "discarded" / "20260101T000000").is_dir()
    assert dock_export.EXPORT_DIRNAME in dock.RESERVED_DIRS
    assert dock.DISCARD_DIRNAME in dock.RESERVED_DIRS


# --- failures ---------------------------------------------------------------------------

def test_incomplete_target_flagged(run_folder, monkeypatch):
    layout, _ = run_folder
    monkeypatch.setenv("FAKE_SMINA_SKIP_IDS", "3,5,7")   # r0 loses 3 ligands
    counts = dock.run(layout, _args(), Prog())
    assert counts["complete"] == 0 and counts["incomplete_targets"] == 2
    assert counts["failed"] == 0   # no unit failed; the ligands just have no poses
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str)
    rows = list(failed.set_index("item").loc[LABELS[0], "reason"])
    # EVERY replicate is checked for ligand coverage, not only r0: a target whose r1 holds
    # a fraction of the library must not pass as finished (the mean over replicates would
    # average it in silently).
    assert rows == [f"incomplete: 3 ligand ID(s) without poses in r{i}: 3,5,7; "
                    f"log: 07_docking/2mjw_A_P_1/replicates/r{i}.log" for i in range(3)]
    assert not (layout.docking / f"df_results_{LABELS[0]}_fixed.csv").exists()
    assert not (layout.docking / LABELS[0] / "poses.sdf").exists()


def test_a_failed_unit_leaves_the_target_incomplete(run_folder, monkeypatch, tmp_path):
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())   # complete first: stale results must go
    bad = tmp_path / "bad_smina"
    bad.write_text("#!/bin/sh\ncase \"$*\" in *--version*) echo 'smina fake 0.0';; "
                   "*r1.sdf*) echo boom >&2; exit 3;; *) exec " + str(FAKE_SMINA) +
                   " \"$@\";; esac\n")
    bad.chmod(0o755)
    monkeypatch.setenv("SMINA_BIN", str(bad))
    for f in layout.docking.glob("*/replicates/r1.done"):
        f.unlink()
    counts = dock.run(layout, _args(), Prog())
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str).set_index("item")
    assert failed.loc[f"{LABELS[0]}/r1", "reason"] == (
        "smina exited with 3: boom; log: 07_docking/2mjw_A_P_1/replicates/r1.log")
    assert failed.loc[LABELS[0], "reason"] == "replicate(s) not recorded: r1"
    assert "boom" in (layout.docking / LABELS[0] / "replicates" / "r1.log").read_text()
    assert counts["incomplete_targets"] == 2 and counts["failed"] == 2
    assert not list(layout.docking.glob("df_results_*"))


def test_missing_receptor_file_is_a_failed_target(run_folder):
    layout, _ = run_folder
    (layout.run_dir / "02_structures/P01350/AF-P01350-F1.pdb").unlink()
    counts = dock.run(layout, _args(), Prog())
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str).set_index("item")
    assert "missing file: 02_structures/P01350/AF-P01350-F1.pdb" in failed.loc[
        LABELS[1], "reason"]
    assert counts["complete"] == 1 and counts["incomplete_targets"] == 1


def test_unit_timeout_kills_the_engine_and_fails_the_unit(run_folder, monkeypatch):
    import time

    layout, _ = run_folder
    monkeypatch.setenv("FAKE_SMINA_SLEEP_ON", "r1.sdf")   # sleeps 60 s
    t0 = time.monotonic()
    counts = dock.run(layout, _args(unit_timeout=1), Prog())
    assert time.monotonic() - t0 < 40
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str).set_index("item")
    for label in LABELS:
        assert failed.loc[f"{label}/r1", "reason"] == (
            f"timeout after 1 s; log: 07_docking/{label}/replicates/r1.log")
        log = (layout.docking / label / "replicates" / "r1.log").read_text()
        assert "Using random seed" in log and "timeout after 1 s" in log
        assert not (layout.docking / label / "replicates" / "r1.done").exists()
        assert (layout.docking / label / "replicates" / "r2.done").exists()   # others ran
    assert counts["failed"] == 2 and counts["incomplete_targets"] == 2
    assert not list(layout.docking.glob("df_results_*"))
    # scoped to THIS test's run folder (it is in every engine call's argv): a bare
    # `pgrep -f <FAKE_SMINA>` also matches a concurrent pytest's engine children.
    assert subprocess.run(["pgrep", "-f", f"{FAKE_SMINA}.*{layout.run_dir}"],
                          capture_output=True).returncode == 1


def test_cli_dry_run_and_library_outside_inputs(run_folder, tmp_path, capsys):
    from polyphar import cli

    layout, _ = run_folder
    outside = tmp_path / "elsewhere.sdf"
    shutil.copyfile(HEAD20, outside)
    argv = ["dock", "--run-dir", str(layout.run_dir), "--nprocs", "4", "--no-notify",
            "--library", str(outside), "--replicates", "2"]
    assert cli.main(argv + ["--dry-run"]) == 0
    assert "dock: 80 units (ligand)" in capsys.readouterr().out
    assert not (layout.inputs / "elsewhere.sdf").exists()   # dry-run writes nothing
    assert cli.main(argv) == 0
    assert (layout.inputs / "elsewhere.sdf").read_bytes() == HEAD20.read_bytes()
    m = json.loads((layout.docking / "manifest.json").read_text())
    assert m["args"]["library"] == "inputs/elsewhere.sdf" and m["args"]["replicates"] == 2


def _pose(lid, score=-5.0):
    return (f"EDRUG3D_{lid}\n x\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
            f"> <minimizedAffinity>\n{score}\n\n$$$$\n")


def test_parse_poses_requires_contiguous_runs_of_at_most_num_modes(tmp_path):
    lib = dock.read_library(HEAD20)
    f = tmp_path / "p.sdf"
    f.write_text(_pose(1) + _pose(1) + _pose(2))
    assert dock.parse_poses(f, lib, "L").Pose.tolist() == [1, 2, 1]
    f.write_text(_pose(1) + _pose(2) + _pose(1))
    with pytest.raises(ValueError, match="ligand ID 1 are not contiguous"):
        dock.parse_poses(f, lib, "L")
    f.write_text(_pose(3) * 13)
    with pytest.raises(ValueError, match="ligand ID 3 has more than 12 poses"):
        dock.parse_poses(f, lib, "L")
    f.write_text(_pose(3) * 12)
    assert len(dock.parse_poses(f, lib, "L")) == 12


def test_result_intact_rejects_an_unterminated_or_empty_record(tmp_path):
    f = tmp_path / "r.sdf"
    assert not dock.result_intact(f)
    f.write_text("")
    assert not dock.result_intact(f)
    f.write_text(_pose(1)[:-6])
    assert not dock.result_intact(f)
    f.write_text(_pose(1))
    assert dock.result_intact(f)


# --- the split stage: export / worker / merge --------------------------------------------

def _bundle_args(bundle=None, **kw):
    base = dict(bundle=bundle, max_units=None, engine_version=None, unit_timeout=None)
    base.update(kw)
    return _args(**base)


def test_export_worker_merge_round_trip(run_folder, tmp_path):
    """The three-way split must reach exactly the state `dock` would, without the merging
    host ever running the engine."""
    from polyphar.stages import dock_export, dock_merge, dock_worker

    layout, calls = run_folder
    bundle = tmp_path / "bundle"
    args = _bundle_args(str(bundle))
    assert dock_export.units(layout, args).count == 120
    out = dock_export.run(layout, args, Prog())
    assert out["units"] == 6
    assert not calls.exists() or calls.read_text() == ""   # export runs no engine
    units = dock.read_bundle_units(bundle)
    assert {(u.unit.label, u.unit.replicate) for u in units} == {
        (l, i) for l in LABELS for i in range(3)}
    # the keys in the bundle are the plan's keys, not recomputed ones
    plan = _plan_of(layout, args)
    assert {u.unit.key for u in units} == {u.key for u in plan.units}
    assert (bundle / dock.BUNDLE_LIBRARY).is_file()
    assert json.loads((bundle / dock.BUNDLE_KEY).read_text())["n_ligands"] == 20

    assert dock_worker.units(layout, args).count == 120
    wc = dock_worker.run(layout, args, Prog())
    assert wc["units_docked"] == 6 and wc["failed"] == 0
    assert len(calls.read_text().splitlines()) == 6
    assert not list(layout.docking.glob("*/replicates/*.sdf"))   # nothing merged yet

    assert dock_merge.units(layout, args).count == 6
    mc = dock_merge.run(layout, args, Prog())
    assert mc["merged"] == 6 and mc["rejected"] == 0
    assert mc["complete"] == 2 and mc["incomplete_targets"] == 0
    for label in LABELS:
        assert (layout.docking / f"df_results_{label}_fixed.csv").is_file()
        for i in range(3):
            assert (layout.docking / label / "replicates" / f"r{i}.done").is_file()
    assert dock.units(layout, _args()).count == 0        # the run folder is now finished


def test_export_skips_recorded_replicates_and_caps_with_max_units(run_folder, tmp_path):
    from polyphar.stages import dock_export

    layout, _ = run_folder
    dock.run(layout, _args(replicates=1), Prog())        # r0 of both targets recorded
    bundle = tmp_path / "b"
    out = dock_export.run(layout, _bundle_args(str(bundle), max_units=3), Prog())
    assert out["units"] == 3 and out["outstanding"] == 4   # r1, r2 of two targets
    got = [(u.unit.label, u.unit.replicate) for u in dock.read_bundle_units(bundle)]
    assert 0 not in {r for _l, r in got}                 # never a recorded draw
    assert [r for _l, r in got] == [1, 1, 2]             # replicate-major


def test_merge_refuses_units_whose_key_is_not_this_runs(run_folder, tmp_path):
    """A bundle made under other settings describes different draws. Merging it would mix
    two populations under one target."""
    from polyphar.stages import dock_export, dock_merge, dock_worker

    layout, _ = run_folder
    bundle = tmp_path / "b"
    args = _bundle_args(str(bundle))
    dock_export.run(layout, args, Prog())
    dock_worker.run(layout, args, Prog())
    # the run folder now wants a different base seed: every key moves
    other = _bundle_args(str(bundle), seed=900)
    counts = dock_merge.run(layout, other, Prog())
    assert counts["merged"] == 0 and counts["rejected"] == 6
    assert not list(layout.docking.glob("*/replicates/*.sdf"))
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str)
    rejected = failed[failed.stage == "dock-merge"]
    assert len(rejected) == 6
    assert rejected.reason.str.contains("unit key differs").all()


def test_merge_never_replaces_a_recorded_draw(run_folder, tmp_path):
    from polyphar.stages import dock_export, dock_merge, dock_worker

    layout, _ = run_folder
    bundle = tmp_path / "b"
    args = _bundle_args(str(bundle))
    dock_export.run(layout, args, Prog())
    dock_worker.run(layout, args, Prog())
    dock_merge.run(layout, args, Prog())
    kept = (layout.docking / LABELS[0] / "replicates" / "r0.sdf").read_bytes()
    # the bundle's copy is tampered with; a second merge must not import it over the record
    (bundle / "results" / LABELS[0] / "replicates" / "r0.sdf").write_text("EDRUG3D_1\n$$$$\n")
    counts = dock_merge.run(layout, args, Prog())
    assert counts["merged"] == 0
    assert (layout.docking / LABELS[0] / "replicates" / "r0.sdf").read_bytes() == kept


def test_worker_skips_done_units_and_redraws_corrupt_ones(run_folder, tmp_path):
    from polyphar.stages import dock_export, dock_worker

    layout, calls = run_folder
    bundle = tmp_path / "b"
    args = _bundle_args(str(bundle))
    dock_export.run(layout, args, Prog())
    dock_worker.run(layout, args, Prog())
    calls.write_text("")
    assert dock_worker.run(layout, args, Prog())["units_docked"] == 0   # all skipped
    assert calls.read_text() == ""
    (bundle / "results" / LABELS[1] / "replicates" / "r2.sdf").write_text("truncated\n")
    counts = dock_worker.run(layout, args, Prog())
    assert counts["units_docked"] == 1 and counts["skipped_existing"] == 5
    log = pd.read_csv(bundle / "results" / dock.DISCARD_LOG, sep="\t", dtype=str)
    assert log.label.tolist() == [LABELS[1]] and log.replicate.tolist() == ["2"]


def test_cli_wires_the_split_subcommands():
    from polyphar.cli import STAGES, build_parser

    assert STAGES["dock-export"][0] == "polyphar.stages.dock_export"
    assert STAGES["dock-worker"][0] == "polyphar.stages.dock_worker"
    assert STAGES["dock-merge"][0] == "polyphar.stages.dock_merge"
    args = build_parser().parse_args(["dock-export", "--max-units", "7", "--replicates", "5"])
    assert args.max_units == 7 and args.replicates == 5
    # the split commands are NOT pipeline stages: `polyphar run` still runs `dock`
    from polyphar.benchmark import ORDER

    assert [s for s in ORDER if s.startswith("dock")] == ["dock"]


@pytest.mark.docking
def test_real_smina_version_is_available(monkeypatch):
    monkeypatch.delenv("SMINA_BIN", raising=False)
    v = dock.smina_version()
    assert v != "unavailable" and v.startswith("smina")


# --- regressions: a draw must never be destroyed by a pass that cannot replace it -------

def test_merge_of_an_unworked_bundle_destroys_nothing(run_folder, tmp_path):
    """P1, the data-destroying bug. `dock-merge` plans the run with execute=False; the
    planner used to DISCARD every outstanding unit before the task list was filtered, so a
    fully docked run folder whose bundle was exported under another engine string (the
    normal remote path: the engine is on the other machine) lost all six draws on a merge
    that ran no engine and merged nothing."""
    from polyphar.stages import dock_export, dock_merge

    layout, calls = run_folder
    dock.run(layout, _args(), Prog())
    before = {p.relative_to(layout.docking).as_posix(): p.read_bytes()
              for p in sorted(layout.docking.glob("*/replicates/r*.sdf"))}
    assert len(before) == 6
    bundle = tmp_path / "b"
    # exported (and merged) under the REMOTE engine's version, as dock-export is meant to be
    args = _bundle_args(str(bundle), engine_version="Uni-Dock 1.2.0")
    dock_export.run(layout, args, Prog())
    calls.write_text("")
    counts = dock_merge.run(layout, _bundle_args(str(bundle)), Prog())   # worker never ran
    assert counts["merged"] == 0 and counts["discarded"] == 0
    assert calls.read_text() == ""                       # no engine call to replace them
    assert {p.relative_to(layout.docking).as_posix(): p.read_bytes()
            for p in sorted(layout.docking.glob("*/replicates/r*.sdf"))} == before
    assert not (layout.docking / dock.DISCARD_LOG).exists()


def test_a_unit_whose_receptor_preparation_failed_keeps_its_draws(run_folder):
    """The same class: the task list was filtered by `u.label in prepared`, but the
    discard had already happened for every outstanding unit."""
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())
    (layout.docking / LABELS[0] / "replicates" / "r1.done").unlink()   # outstanding again
    kept = (layout.docking / LABELS[0] / "replicates" / "r1.sdf").read_bytes()
    (layout.run_dir / "04_pockets/P25815/P25815_2mjw_A.pdb_res_P_1.pdb").write_text("x\n")
    counts = dock.run(layout, _args(), Prog())       # pocket_box now raises for this target
    assert counts["discarded"] == 0
    assert (layout.docking / LABELS[0] / "replicates" / "r1.sdf").read_bytes() == kept


def test_a_discarded_record_is_archived_not_deleted(run_folder):
    """`_obsolete` moves orphan labels, "never deleted". A per-unit draw is a smaller
    thing to lose, not a less valuable one: it moves to 07_docking/discarded/ and
    discarded.tsv says where."""
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())
    kept = {label: (layout.docking / label / "replicates" / "r0.sdf").read_bytes()
            for label in LABELS}
    dock.run(layout, _args(seed=100), Prog())        # every key changes: all six replaced
    log = _discards(layout)
    assert len(log) == 6 and log.archived.notna().all()
    archived = [layout.docking / a for a in log.archived]
    assert all(p.is_file() for p in archived)
    for label in LABELS:
        assert any(p.read_bytes() == kept[label] for p in archived)


def test_a_bundle_draw_dropped_for_one_already_on_record_is_logged(run_folder, tmp_path):
    """P3. Export only excludes units already ON RECORD, so two bundles exported before
    either is merged hold the same units. The first merge wins; the loser was dropped with
    no row anywhere, which is an unlogged selection effect."""
    from polyphar.stages import dock_export, dock_merge, dock_worker

    layout, _ = run_folder
    b1, b2 = tmp_path / "b1", tmp_path / "b2"
    dock_export.run(layout, _bundle_args(str(b1)), Prog())
    dock_export.run(layout, _bundle_args(str(b2)), Prog())   # nothing merged yet: identical
    for b in (b1, b2):
        dock_worker.run(layout, _bundle_args(str(b)), Prog())
    assert dock_merge.run(layout, _bundle_args(str(b1)), Prog())["merged"] == 6
    counts = dock_merge.run(layout, _bundle_args(str(b2)), Prog())
    assert counts["merged"] == 0 and counts["dropped"] == 6
    log = _discards(layout)
    assert len(log) == 6
    assert all("dropped on merge" in r for r in log.reason)
    assert all(str(b2) in r for r in log.reason)
    assert all(Path(a).is_file() for a in log.archived)   # the draw is still in the bundle
    # merging the same bundle again logs the same drop once, not twice
    dock_merge.run(layout, _bundle_args(str(b2)), Prog())
    assert len(_discards(layout)) == 6


def test_merged_manifest_names_the_engine_that_actually_drew_the_results(run_folder,
                                                                        tmp_path):
    """P4. The keys are planned under the bundle's run_key.json engine (they must be, or
    they would not be the keys the worker was given), but the Methods source has to name
    the engine that RAN. worker.json held it and nothing read it."""
    from polyphar.stages import dock_export, dock_merge, dock_worker

    layout, _ = run_folder
    bundle = tmp_path / "b"
    dock_export.run(layout, _bundle_args(str(bundle), engine_version="smina fake 0.0"),
                    Prog())
    shim = tmp_path / "unidock"
    shim.write_text('#!/bin/sh\ncase "$*" in *--version*) echo "Uni-Dock 1.2.0";; '
                    f'*) exec {FAKE_SMINA} "$@";; esac\n')
    shim.chmod(0o755)
    import os as _os

    old = _os.environ["SMINA_BIN"]
    _os.environ["SMINA_BIN"] = str(shim)
    try:
        dock_worker.run(layout, _bundle_args(str(bundle)), Prog())
    finally:
        _os.environ["SMINA_BIN"] = old
    counts = dock_merge.run(layout, _bundle_args(str(bundle)), Prog())
    m = json.loads((layout.docking / "manifest.json").read_text())
    assert m["engine"] == "smina fake 0.0"            # what the keys were planned under
    assert m["engine_worker"] == "Uni-Dock 1.2.0"     # what actually docked
    assert m["merged_from"] == str(bundle)
    assert m["counts"]["merged"] == counts["merged"] == 6
    assert m["counts"]["rejected"] == 0


def test_every_replicate_is_checked_for_ligand_coverage(run_folder):
    """P5. `result_intact` only reads the last line, so a record truncated on a clean
    `$$$$` boundary passes it, and only r0 was ever parsed. A 15%-complete r1 stayed
    `done`, with no failed row -- and would have entered the mean over replicates."""
    layout, _ = run_folder
    assert dock.run(layout, _args(), Prog())["complete"] == 2
    r1 = layout.docking / LABELS[0] / "replicates" / "r1.sdf"
    records = r1.read_text().split("$$$$\n")
    r1.write_text("$$$$\n".join(records[:6]) + "$$$$\n")   # 6 of 40 pose records
    assert dock.result_intact(r1)                          # still "intact"
    assert dock.unit_state(layout, _plan_of(layout, _args()).units_of(LABELS[0])[1]) == \
        dock.STATE_DONE
    (layout.docking / "manifest.json").unlink()
    counts = dock.run(layout, _args(), Prog())
    assert counts["complete"] == 1 and counts["incomplete_targets"] == 1
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str)
    assert any("without poses in r1" in r for r in failed.reason)


def test_lowering_replicates_keeps_the_surplus_records_and_counts_them(run_folder):
    """A recorded draw is never thrown away for being surplus -- but nothing used to say
    it was still there."""
    layout, _ = run_folder
    dock.run(layout, _args(replicates=5), Prog())
    counts = dock.run(layout, _args(replicates=3), Prog())
    assert counts["units_present_outside_plan"] == 4     # r3, r4 of both targets
    assert len(list(layout.docking.glob("*/replicates/r*.sdf"))) == 10


def test_the_receptor_key_stamp_follows_the_openbabel_version(run_folder, monkeypatch):
    """P7. The PDBQT the engine reads is deliberately not in the unit key (a Babel upgrade
    must not invalidate every draw), so the stamp has to carry the version or a resumed
    folder and a fresh one would feed different receptors under identical keys."""
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())
    target = dock.read_targets(layout)[0]
    stamp = layout.docking / LABELS[0] / "receptor.key"
    before, src_before = stamp.read_text(), dock._src_key(target)
    monkeypatch.setattr(dock, "_openbabel_version", lambda: "9.9.9")
    assert dock._src_key(target) == src_before           # unit keys are NOT disturbed
    assert dock._prep_stamp(target) != before.strip()    # but the stamp is
    dock._prepare_receptor(layout, target)               # so the receptor re-prepares
    assert stamp.read_text() != before


def test_a_changed_receptor_pdbqt_is_warned_about_against_the_last_manifest(run_folder,
                                                                           capsys):
    """The other half of P7: the PDBQT is out of the key on purpose, so a change in it has
    to be reported instead of silently absorbed."""
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())
    pdbqt = layout.docking / LABELS[0] / "receptor.pdbqt"
    pdbqt.write_text(pdbqt.read_text() + "REMARK tampered\n")
    capsys.readouterr()
    dock.run(layout, _args(), Prog())
    assert "receptor PDBQT changed" in capsys.readouterr().err
    assert LABELS[0] in json.loads(
        (layout.docking / "manifest.json").read_text())["receptors"]
