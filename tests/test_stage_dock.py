import json
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

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


@pytest.mark.docking
def test_cpu_and_chunk_invariance(tmp_path):
    """One job at cpu=12 equals 2 chunks at cpu=1, the second in reversed order: a
    ligand's result depends neither on --cpu nor on its position or chunk mates."""
    ids = [3, 7, 8, 9]   # small ligands: histamine, guanidine, menadione, theophylline
    rec = tmp_path / "rec.pdbqt"
    dock.receptor_to_pdbqt(RECEPTOR, rec)
    center, size = dock.pocket_box(POCKET)
    lib = dock.read_library(HEAD20)

    def table(ligs: Path, cpu: int):
        out = ligs.with_suffix(".out.sdf")
        dock.run_smina(rec, ligs, center, size, out, ligs.with_suffix(".log"), seed=42,
                       cpu=cpu)
        return dock.parse_poses(out, lib, "L")

    one = table(_subset(HEAD20, ids, tmp_path / "all.sdf"), 12)
    two = pd.concat([table(_subset(HEAD20, ids[:2], tmp_path / "c0.sdf"), 1),
                     table(_subset(HEAD20, ids[:1:-1], tmp_path / "c1.sdf"), 1)])
    key = ["ID", "Pose"]
    one = one.sort_values(key).reset_index(drop=True)
    two = two.sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(one, two)


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
    base = dict(nprocs=2, library="inputs/lib.sdf", seed=42, chunk_size=8, threads_per_job=1,
                dry_run=False)
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


def test_run_writes_results_per_target(run_folder):
    layout, calls = run_folder
    args = _args()
    assert dock.units(layout, args).count == 40   # 2 targets x 20 ligands
    prog = Prog()
    counts = dock.run(layout, args, prog)
    assert sum(prog.calls) == 40
    assert counts["complete"] == 2 and counts["failed"] == 0 and counts["chunks"] == 6
    for label in ("2mjw_A_P_1", "AF-P01350-F1_P_1"):
        df = pd.read_csv(layout.docking / f"df_results_{label}_fixed.csv", dtype={"ID": str})
        assert df.columns.tolist() == ["Pose", "SMINA_Score", "Ligand", "ID", "Pocket"]
        assert df.ID.tolist() == [str(i) for i in range(1, 21) for _ in (1, 2)]
        assert set(df.Pocket) == {label}
        assert df.Ligand[0] == "SULFAPYRIDINE"
        assert (layout.docking / label / "poses.sdf").is_file()
    assert (layout.docking / "AF-P01350-F1_P_1" / "receptor_h.pdb").is_file()
    assert not (layout.docking / "2mjw_A_P_1" / "receptor_h.pdb").exists()
    assert {l.split()[1] for l in calls.read_text().splitlines()} == {"cpu=1"}
    m = json.loads((layout.docking / "manifest.json").read_text())
    assert m["seed"] == 42 and m["smina"] == "smina fake 0.0"
    assert m["library_md5"] == dock._md5(HEAD20)
    rec = m["receptors"]["2mjw_A_P_1"]
    assert rec["receptor"] == "03_pdbfixer/P25815/2mjw_A.pdb"
    assert len(rec["pdbqt_sha256"]) == 64 and "pdbfixer pH 7.0" in rec["preparation"]
    assert "addMissingHydrogens" in m["receptors"]["AF-P01350-F1_P_1"]["preparation"]
    assert (layout.docking / "failed.tsv").read_text() == "item\tstage\treason\n"

    # finished on the same inputs: nothing to do, nothing changes (any nprocs)
    before = calls.read_text()
    assert dock.units(layout, _args(nprocs=5)).count == 0
    assert dock.run(layout, _args(nprocs=5), Prog()) == counts
    assert calls.read_text() == before
    # a different seed redocks everything
    assert dock.units(layout, _args(seed=7)).count == 40


def test_results_identical_for_any_nprocs_and_chunking(run_folder, tmp_path):
    layout, _ = run_folder
    dock.run(layout, _args(nprocs=1, chunk_size=20), Prog())
    a = (layout.docking / "df_results_2mjw_A_P_1_fixed.csv").read_bytes()
    shutil.rmtree(layout.docking)
    dock.run(layout, _args(nprocs=3, chunk_size=3), Prog())
    assert (layout.docking / "df_results_2mjw_A_P_1_fixed.csv").read_bytes() == a


def test_incomplete_target_flagged(run_folder, monkeypatch):
    layout, _ = run_folder
    monkeypatch.setenv("FAKE_SMINA_SKIP_IDS", "3,5,7")   # one chunk loses 3 ligands
    counts = dock.run(layout, _args(), Prog())
    assert counts["complete"] == 0 and counts["incomplete"] == 2
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str)
    row = failed.set_index("item").loc["2mjw_A_P_1", "reason"]
    assert row.endswith("3 ligand ID(s) without poses: 3,5,7")
    assert not (layout.docking / "df_results_2mjw_A_P_1_fixed.csv").exists()
    assert not (layout.docking / "2mjw_A_P_1" / "poses.sdf").exists()


def test_failed_chunk_marks_target_incomplete(run_folder, monkeypatch, tmp_path):
    layout, _ = run_folder
    dock.run(layout, _args(), Prog())   # complete first: stale results must go
    bad = tmp_path / "bad_smina"
    bad.write_text("#!/bin/sh\ncase \"$*\" in *--version*) echo 'smina fake 0.0';; "
                   "*c8_00001*) echo boom >&2; exit 3;; *) exec " + str(FAKE_SMINA) + " \"$@\";; esac\n")
    bad.chmod(0o755)
    monkeypatch.setenv("SMINA_BIN", str(bad))
    for f in layout.docking.glob("*/chunks/c8_00001.done"):
        f.unlink()
    counts = dock.run(layout, _args(), Prog())
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str).set_index("item")
    assert "smina exited with 3: boom" in failed.loc["2mjw_A_P_1/c8_00001", "reason"]
    assert failed.loc["2mjw_A_P_1", "reason"].endswith(
        "8 ligand ID(s) without poses: 9,10,11,12,13,14,15,16")
    assert counts["incomplete"] == 2 and counts["failed"] == 4
    assert not list(layout.docking.glob("df_results_*"))


def test_resume_skips_done_chunks(run_folder):
    layout, calls = run_folder
    dock.run(layout, _args(), Prog())
    (layout.docking / "manifest.json").unlink()            # as if killed before the end
    (layout.docking / "2mjw_A_P_1" / "chunks" / "c8_00002.done").unlink()   # 4 ligands
    (layout.docking / "AF-P01350-F1_P_1" / "chunks" / "c8_00000.done").write_text("stale\n")
    assert dock.units(layout, _args()).count == 4 + 8
    calls.write_text("")
    prog = Prog()
    counts = dock.run(layout, _args(), prog)
    assert sorted(prog.calls) == [4, 8]
    assert len(calls.read_text().splitlines()) == 2
    assert counts["skipped_existing"] == 4 and counts["complete"] == 2


def test_missing_receptor_file_is_a_failed_target(run_folder):
    layout, _ = run_folder
    (layout.run_dir / "02_structures/P01350/AF-P01350-F1.pdb").unlink()
    counts = dock.run(layout, _args(), Prog())
    failed = pd.read_csv(layout.docking / "failed.tsv", sep="\t", dtype=str).set_index("item")
    assert "missing file: 02_structures/P01350/AF-P01350-F1.pdb" in failed.loc[
        "AF-P01350-F1_P_1", "reason"]
    assert counts["complete"] == 1 and counts["incomplete"] == 1


def test_cli_dry_run_and_library_outside_inputs(run_folder, tmp_path, capsys):
    from polyphar import cli

    layout, _ = run_folder
    outside = tmp_path / "elsewhere.sdf"
    shutil.copyfile(HEAD20, outside)
    argv = ["dock", "--run-dir", str(layout.run_dir), "--nprocs", "4", "--no-notify",
            "--library", str(outside), "--chunk-size", "5"]
    assert cli.main(argv + ["--dry-run"]) == 0
    assert "dock: 40 units (ligand)" in capsys.readouterr().out
    assert not (layout.inputs / "elsewhere.sdf").exists()   # dry-run writes nothing
    assert cli.main(argv) == 0
    assert (layout.inputs / "elsewhere.sdf").read_bytes() == HEAD20.read_bytes()
    m = json.loads((layout.docking / "manifest.json").read_text())
    assert m["args"]["library"] == "inputs/elsewhere.sdf" and m["args"]["chunk_size"] == 5
