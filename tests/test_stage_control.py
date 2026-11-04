"""`polyphar control`: the positive control the stats stage already enforces.

The readout maths (self-rank denominator, pose-in-pocket geometry, the JSON `stats`
accepts) is tested on constructed inputs with no engine and no Open Babel, so it runs in
any environment. The end-to-end stage test needs Open Babel (the co-crystal SDF) and uses
tests/fixtures/fake_smina_posed, a stub that writes poses with coordinates.
"""

import csv
import json
import math
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

from polyphar.common import RunLayout
from polyphar.stages import control, stats

FIX = Path(__file__).parent / "fixtures"
RECEPTOR = FIX / "dock_2mjw" / "2mjw_A.pdb"
POCKET = FIX / "dock_2mjw" / "P25815_2mjw_A.pdb_res_P_1.pdb"
AF_MODEL = FIX / "dock_af" / "AF-P01350-F1.pdb"
FAKE_POSED = FIX / "fake_smina_posed"


class Prog:
    def __init__(self):
        self.calls = []

    def advance(self, n=1, size=None):
        self.calls.append(n)


def _args(**kw):
    base = dict(nprocs=2, seed=42, control_replicates=1, threads_per_job=1,
                unit_timeout=600, engine="smina", search_mode=None, batch_size=None,
                dry_run=False, no_notify=True,
                positive_control_out=control.INPUTS_POSITIVE_CONTROL)
    base.update(kw)
    return SimpleNamespace(**base)


# --- readout 1: the self-rank denominator ------------------------------------------------

def test_self_rank_is_against_the_whole_panel_not_a_subset():
    """The home pocket's rank is a percentile of EVERY target the ligand was docked
    against. Best of 100 is 0.005 (top 1%); the same score against 10 targets is 0.05.
    Shrinking the denominator is exactly how a self-rank is made to look good."""
    big = {f"t{i}": -1.0 - i * 0.01 for i in range(100)}
    big["home"] = -9.0                       # better than all 100
    assert control.self_rank_percentile("home", big) == pytest.approx(0.5 / 101)

    small = {f"t{i}": -1.0 - i * 0.01 for i in range(10)}
    small["home"] = -9.0
    assert control.self_rank_percentile("home", small) == pytest.approx(0.5 / 11)
    # the SAME home score, ranked against 11 instead of 101 targets, is 10x worse
    assert control.self_rank_percentile("home", small) > 0.04


def test_self_rank_of_the_worst_and_of_chance():
    scores = {f"t{i}": -float(i) for i in range(1, 101)}
    scores["home"] = 0.0                     # worse than every one of them
    assert control.self_rank_percentile("home", scores) == pytest.approx(100.5 / 101)
    # a ligand tied with everything sits at chance, not at the top
    tied = {f"t{i}": -5.0 for i in range(100)}
    tied["home"] = -5.0
    assert control.self_rank_percentile("home", tied) == pytest.approx(0.5)


def test_a_ligand_missing_from_any_panel_target_gets_no_self_rank():
    """Not a smaller denominator: no rank at all, with the reason recorded."""
    labels = [control.Label("A_P_1", "P1", "04_pockets/P1/a.pdb", "lig1", "LIG", "A", "1",
                            20, 0.9, Path("x"), ("LIG", "A", "1")),
              control.Label("B_P_1", "P2", "04_pockets/P2/b.pdb", "lig2", "LIG", "A", "1",
                            20, 0.9, Path("x"), ("LIG", "A", "1"))]
    scores = {"A_P_1": {"lig1": -9.0, "lig2": -1.0},
              "B_P_1": {"lig1": -1.0},          # lig2 missing here
              "C_P_1": {"lig1": -1.0, "lig2": -1.0}}
    rows, reasons = control.self_ranks(labels, scores)
    assert [r["Ligand_ID"] for r in rows] == ["lig1"]
    assert rows[0]["Panel_Targets"] == 3
    assert reasons and reasons[0][0] == "lig2" and "no self-rank" in reasons[0][1]


# --- readout 2: pose-in-pocket geometry ---------------------------------------------------

def _shell(radius=5.0, centre=(0.0, 0.0, 0.0)):
    """A cavity lining: atoms on a sphere around `centre`."""
    out = []
    for sx in (-1, 1):
        for k in range(3):
            p = [centre[0], centre[1], centre[2]]
            p[k] += sx * radius
            out.append(tuple(p))
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                r = radius / math.sqrt(3)
                out.append((centre[0] + sx * r, centre[1] + sy * r, centre[2] + sz * r))
    return out


def _box(xyz):
    lo = [min(c[k] for c in xyz) for k in range(3)]
    hi = [max(c[k] for c in xyz) for k in range(3)]
    return (tuple((hi[k] + lo[k]) / 2 for k in range(3)),
            tuple(hi[k] - lo[k] for k in range(3)))


def test_a_pose_in_the_middle_of_the_cavity_is_in_the_pocket():
    lining = _shell()
    ok, n_near, asym, why = control.pose_in_pocket((0.0, 0.0, 0.0), lining, *_box(lining))
    assert ok and why == ""
    assert n_near == len(lining)
    assert asym == pytest.approx(0.0, abs=1e-9)


def test_a_pose_far_from_the_cavity_is_not_in_the_pocket():
    lining = _shell()
    ok, n_near, _asym, why = control.pose_in_pocket((60.0, 0.0, 0.0), lining, *_box(lining))
    assert not ok
    assert n_near == 0
    assert "open space" in why


@pytest.mark.parametrize("standoff", [3.2, 3.6, 4.0, 4.5])
def test_a_pose_on_a_real_surface_patch_is_not_in_the_pocket(standoff):
    """All of the lining on one side of the centroid: a molecule lying ON the surface,
    which is the failure mode the asymmetry criterion exists to catch.

    On REAL lining, at a REALISTIC standoff. The original form of this test used a flat
    synthetic 15-atom patch with the centroid 6.0 A directly above its centre, and no
    ligand sits 6 A off a surface with nothing else nearby: it proved the criterion rejects
    a case that cannot occur. On that synthetic geometry the asymmetry is exactly
    standoff/SHELL_RADIUS, so it does NOT clear MAX_ASYMMETRY below 4 A -- the flat patch
    never established anything about van der Waals contact. Real lining is not a plane: a
    centroid 3.2-4.5 A outside the outermost lining atom of the 2mjw pocket sees its
    neighbours at asymmetry 0.67-0.79 and is rejected across the whole range."""
    lining = control.pocket_atoms(POCKET, hetatm=True)
    centre, size = __import__("polyphar.stages.dock", fromlist=["x"]).pocket_box(POCKET)
    mid = control.centroid(lining)
    far = max(lining, key=lambda p: math.dist(p, mid))
    d = math.dist(far, mid)
    probe = tuple(far[k] + (far[k] - mid[k]) / d * standoff for k in range(3))
    ok, n_near, asym, why = control.pose_in_pocket(probe, lining, centre, size)
    assert not ok, f"asymmetry {asym} over {n_near} lining atoms"
    assert asym > control.MAX_ASYMMETRY
    assert "surface patch" in why


def test_a_pose_pressed_against_a_box_face_is_not_in_the_pocket():
    lining = _shell()
    centre, size = _box(lining)
    # enclosed enough to pass the asymmetry test, but hard against the +x face
    on_face = (size[0] / 2 - 0.1, 0.0, 0.0)
    ok, _n, asym, why = control.pose_in_pocket(on_face, lining, centre, size)
    assert not ok and "box edge" in why
    assert asym <= control.MAX_ASYMMETRY      # the EDGE test is what rejected it


def test_the_best_pose_is_chosen_on_score_not_on_record_order(tmp_path):
    sdf = tmp_path / "poses.sdf"
    def rec(lid, score, x):
        return (f"EDRUG3D_{lid}\n t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                f"{x:10.4f}{0.0:10.4f}{0.0:10.4f} C   0  0  0  0  0  0\nM  END\n"
                f"> <minimizedAffinity>\n{score}\n\n$$$$\n")
    sdf.write_text(rec("L1", -1.0, 50.0) + rec("L1", -7.0, 1.0))
    best = control.best_poses(sdf, "minimizedAffinity")
    assert best["L1"][0] == -7.0
    assert best["L1"][1] == [(1.0, 0.0, 0.0, "C")]
    assert best["L1"][2] == ""


# --- the JSON `stats` accepts --------------------------------------------------------------

def _rank_rows(n, percentile):
    return [{"Ligand_ID": f"l{i}", "Ligand": "LIG", "Home_Target": f"t{i}",
             "Home_Score": -9.0, "Best_Score": -9.0, "Panel_Targets": 200,
             "Self_Rank": 1, "Self_Rank_Percentile": percentile} for i in range(n)]


def _pose_rows(n_in, n_out, self_pairs=0, rmsd=0.5):
    """`self_pairs` of the in-pocket rows are ALSO self pairs carrying an RMSD, which is
    what readout 2 gates on; the rest are the non-cognate pairs, which carry none."""
    rows = [{"Target": f"t{i}", "Ligand_ID": f"l{i}", "Score": -9.0, "In_Pocket": True,
             "Contact_Fraction": 0.9,
             "Self_Pair": i < self_pairs,
             "RMSD": rmsd if i < self_pairs else "",
             "Redock_Success": (rmsd <= control.REDOCK_RMSD_MAX) if i < self_pairs else ""}
            for i in range(n_in)]
    rows += [{"Target": f"u{i}", "Ligand_ID": f"l{i}", "Score": -1.0, "In_Pocket": False,
              "Contact_Fraction": 0.1, "Self_Pair": False, "RMSD": "",
              "Redock_Success": ""} for i in range(n_out)]
    return rows


def _lib_rows(n, pct=0.01):
    return [{"Ligand_ID": f"l{i}", "Ligand": "LIG", "Home_Target": f"t{i}",
             "Control_Score": -9.0, "Library_Ligands": 1960, "Better_In_Library": 10,
             "Percentile": pct, "Top_Decile": pct <= control.TOP_DECILE, "Reason": ""}
            for i in range(n)]


def _result(n=20, percentile=0.004, n_in=85, n_out=15, self_pairs=20, rmsd=0.5,
            lib_pct=0.01, n_panel=200, n_total=200):
    return control.positive_control_json(
        _rank_rows(n, percentile), _pose_rows(n_in, n_out, self_pairs, rmsd), [],
        n_panel, "polyphar control: a test",
        lib_rows=_lib_rows(n, lib_pct), n_panel_total=n_total)


def test_a_passing_control_validates_against_what_stats_accepts(tmp_path):
    result = _result()
    path = tmp_path / "positive_control.json"
    path.write_text(json.dumps(result))
    verdict, detail = stats.positive_control_verdict(path)
    assert verdict == "yes"
    # the loop the task asked to close: a passing control that states no `source` makes
    # `stats` write a failed.tsv line. This one states one.
    assert detail["source"].startswith("polyphar control")
    assert result["n_redocked"] == 20
    assert result["median_self_rank_percentile"] == pytest.approx(0.004)
    assert result["redock_rmsd_success_fraction"] == pytest.approx(1.0)
    assert result["library_top_decile_fraction"] == pytest.approx(1.0)
    assert result["pose_in_pocket_fraction"] == pytest.approx(0.85)


def test_a_failing_control_is_uninterpretable_not_a_near_miss(tmp_path):
    result = _result(percentile=0.4, n_in=50, n_out=50, rmsd=9.0, lib_pct=0.8)
    path = tmp_path / "positive_control.json"
    path.write_text(json.dumps(result))
    verdict, _ = stats.positive_control_verdict(path)
    assert verdict.startswith("no:")
    assert "UNINTERPRETABLE" in verdict


def test_every_readout_the_json_reports_is_a_fraction_in_range():
    """stats range-checks both gated readouts and refuses anything outside [0, 1]; a
    percentile or a fraction produced here can never be out of range."""
    result = _result(n=7, percentile=0.03, n_in=3, n_out=4, self_pairs=3)
    for key in ("median_self_rank_percentile", "redock_rmsd_success_fraction",
                "library_top_decile_fraction", "pose_in_pocket_fraction"):
        assert 0.0 <= result[key] <= 1.0 and not isinstance(result[key], bool)
    sr = result["self_rank"]
    assert 0.0 <= sr["fraction_top_1pct"] <= 1.0
    assert 0.0 <= sr["fraction_top_5pct"] <= 1.0 and sr["fraction_top_5pct"] == 1.0


def test_the_caveat_travels_with_the_number():
    """It is quoted by people who did not read the code, so it is IN the result, not only
    in a docstring: the JSON and the summary both carry it."""
    result = _result(n=5, percentile=0.01, n_in=9, n_out=1, self_pairs=5)
    for phrase in ("solved apo", "LARGE", "WELL-STUDIED", "NOT AN OPTIMISATION TARGET",
                   "HOLO", "EASIEST BENCHMARK", "pseudouridine", "WHAT A PASS LICENSES",
                   "+0.33 kcal/mol smina offset, which cancels exactly"):
        assert phrase in result["caveat"]
    md = control.summary_markdown(result, {"verdict": "yes", **stats.THRESHOLDS})
    assert result["caveat"] in md
    assert "## Caveat" in md


# --- labelling ------------------------------------------------------------------------------

def _pdb(lines):
    return "\n".join(lines) + "\nEND\n"


def _het(resn, chain, num, coords, serial=1):
    return [f"HETATM{serial + i:5d}  C{i:<2} {resn:>3} {chain}{int(num):4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
            for i, (x, y, z) in enumerate(coords)]


def test_solvent_ions_and_small_fragments_are_not_a_co_crystal_label(tmp_path):
    """rg_calibrate.py's SKIP list and the 12-60 heavy-atom window, unchanged: a glycerol,
    a zinc and a 3-atom fragment are not evidence of a drug-binding site."""
    coords = [(float(i) * 0.5, 0.0, 0.0) for i in range(20)]
    chain = tmp_path / "chain.pdb"
    chain.write_text(_pdb(_het("GOL", "A", 1, coords)
                          + _het("ZN", "A", 2, coords[:1], serial=100)
                          + _het("FRG", "A", 3, coords[:3], serial=200)))
    groups = control.het_groups(chain)
    assert ("GOL", "A", "1") not in groups and ("ZN", "A", "2") not in groups
    assert ("FRG", "A", "3") in groups        # not solvent, but too small for the window
    sized = {k: v for k, v in groups.items()
             if control.MIN_HEAVY <= len(v) <= control.MAX_HEAVY}
    assert sized == {}


def test_a_ligand_in_the_pocket_is_labelled_and_one_beside_it_is_not(tmp_path):
    """The label is the contact criterion, not proximity to the file: 70% of the ligand's
    heavy atoms within 4.5 A of the pocket lining."""
    lining = [(float(i), 0.0, 0.0) for i in range(15)]
    inside = [(float(i), 1.0, 0.0) for i in range(14)]        # all within 4.5 A
    far = [(float(i), 40.0, 0.0) for i in range(14)]
    assert control.contact_fraction(inside, lining) == 1.0
    assert control.contact_fraction(far, lining) == 0.0
    assert control.contact_fraction(inside[:4] + far[:10], lining) == pytest.approx(4 / 14)


def test_alphafold_targets_carry_no_co_crystal_label(tmp_path):
    layout = RunLayout(tmp_path)
    from polyphar.stages import dock

    t = dock.Target("AF-P01350-F1_P_1", tmp_path / "x", tmp_path / "y",
                    "02_structures/P01350/AF-P01350-F1.pdb",
                    "04_pockets/P01350/p.pdb", alphafold=True)
    assert control.chain_pdb_of(layout, t) is None
    label, why = control.label_target(layout, t)
    assert label is None and "AlphaFold" in why


def test_the_extracted_chain_is_found_through_receptor_path_not_the_pocket_name(tmp_path):
    """A chimeric chain's pocket file is `<accession>_<stem>_res_P_1.pdb`, so splitting the
    pocket name on `_res_P` yields `<accession>_<stem>` and looks for a chain file that does
    not exist (rg_calibrate.py does this and silently finds no ligand). Receptor_Path does
    not have the problem, and it is also spelling-independent across DoGSite tools."""
    from polyphar.stages import dock

    layout = RunLayout(tmp_path)
    t = dock.Target("2mjw_A_P_1", tmp_path / "r", tmp_path / "p",
                    "03_pdbfixer/P25815/2mjw_A.pdb",
                    "04_pockets/P25815/P25815_2mjw_A_res_P_1.pdb", alphafold=False)
    assert control.chain_pdb_of(layout, t) == tmp_path / "02_structures/P25815/2mjw_A.pdb"


# --- the stage ---------------------------------------------------------------------------

def _co_crystal_chain(src: Path, pocket: Path, dst: Path, n_atoms: int = 16) -> None:
    """`src` with a drug-like HETATM ligand planted at the centroid of `pocket`."""
    xyz = control.pocket_atoms(pocket)
    c = control.centroid(xyz)
    # a compact cluster at the cavity centre: every atom is within 4.5 A of the lining
    coords = [(c[0] + (i % 4) * 0.6 - 0.9, c[1] + (i // 4) * 0.6 - 0.9, c[2])
              for i in range(n_atoms)]
    body = [l for l in src.read_text().splitlines() if not l.startswith("END")]
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(_pdb(body + _het("LIG", "A", 900, coords, serial=9000)))


def test_pdbfixer_leaves_the_extracted_chains_coordinates_where_they_are(tmp_path):
    """The assumption the whole stage rests on: the co-crystal ligand is read from the
    EXTRACTED CHAIN (02_structures) and the pocket from DoGSiteScorer output made from the
    PDBFIXER chain, so the two must share a frame. pockets.py's `_author_atoms`/`_lookup`
    allow 0.001 A for it; here it is measured, and the tolerance is never even exercised --
    every heavy atom comes back byte-identical in the PDB's 3 decimals."""
    pytest.importorskip("pdbfixer", reason="pdbfixer not importable")
    from polyphar.stages.fix import fix_one

    src = tmp_path / "st" / "P25815" / "2mjw_A.pdb"
    src.parent.mkdir(parents=True)
    shutil.copyfile(RECEPTOR, src)
    fix_one(src, tmp_path / "st", tmp_path / "out")

    def heavy(p):
        """(residue, atom name, x, y, z) per heavy atom of the first model, in file
        order. NOT a set of coordinate triples: `set(before) <= set(after)` is containment,
        which two different atoms at the same position would satisfy and which says nothing
        about any particular atom. This is a per-atom identity check."""
        out = []
        for line in p.read_text().splitlines():
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")) and control._is_heavy(line):
                out.append((line[17:20].strip(), line[22:27].strip(), line[12:16].strip(),
                            line[30:38], line[38:46], line[46:54]))
        return out

    before, after = heavy(src), heavy(tmp_path / "out" / "P25815" / "2mjw_A.pdb")
    assert before, "the fixture has no heavy atoms"
    assert len(before) == 788
    # every atom, in order, byte-identical in the PDB's 3 decimals: the 0.001 A tolerance
    # pockets.py allows is never exercised.
    assert before == after
    # Measured the same way out of tree on 2mjw_B (not a fixture: 1.1 MB): all 719 heavy
    # atoms of its first model come back identically, but only once the CHAIN ID is
    # ignored -- pdbfixer renames the chain, and it collapses this NMR entry's 10 models to
    # the first. Neither matters here (the pocket and the ligand are read from the same
    # first model, and the chain ID is not part of any coordinate), but a comparison that
    # includes the chain ID will fail on B and must not be written.


@pytest.mark.openbabel
def test_the_ligand_sdf_does_not_carry_open_babels_clock(tmp_path):
    """Open Babel stamps the molblock's program line with MMDDYYHHmm. Left in, the control
    library's sha256 changes at every minute boundary, which changes dock's unit keys and
    re-docks the whole control on a rerun that happens to cross one."""
    lines = [f"HETATM{9000 + i:5d}  C{i:<2} LIG A 900    {i * 1.4:8.3f}{0.0:8.3f}"
             f"{0.0:8.3f}  1.00  0.00           C" for i in range(14)]
    rec = control.ligand_record(lines, "ID1", "LIG")
    assert "OpenBabel" not in rec
    assert rec.split("\n")[1] == control.PROGRAM_LINE
    assert rec == control.ligand_record(lines, "ID1", "LIG")


@pytest.fixture
def run_folder(tmp_path, monkeypatch):
    """A run folder with one experimental target carrying a planted co-crystal ligand and
    one AlphaFold target that cannot carry one."""
    import pandas as pd

    monkeypatch.setenv("SMINA_BIN", str(FAKE_POSED))
    monkeypatch.setenv("FAKE_POSE_BEST_RECEPTOR", "2mjw")
    # echo mode: the pose is the input ligand's own atom block, which is what a real engine
    # returns and the only thing a redocking RMSD can be computed against.
    monkeypatch.setenv("FAKE_POSE_ECHO", "1")
    monkeypatch.delenv("FAKE_POSE_RMSD", raising=False)
    monkeypatch.delenv("FAKE_POSE_OUTSIDE_IDS", raising=False)
    monkeypatch.delenv("FAKE_POSE_STRONG_IDS", raising=False)
    run = tmp_path / "run"
    layout = RunLayout(run)
    for rel, src in {"03_pdbfixer/P25815/2mjw_A.pdb": RECEPTOR,
                     "04_pockets/P25815/P25815_2mjw_A.pdb_res_P_1.pdb": POCKET,
                     "02_structures/P01350/AF-P01350-F1.pdb": AF_MODEL,
                     "04_pockets/P01350/AF-P01350-F1.pdb_res_P_1.pdb": POCKET}.items():
        (run / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, run / rel)
    _co_crystal_chain(RECEPTOR, POCKET, run / "02_structures/P25815/2mjw_A.pdb")
    layout.targets.mkdir(parents=True)
    pd.DataFrame({
        "Node": ["n1", "n2"], "UniProt_ID": ["P25815", "P01350"],
        "Gene_Name": ["S100P", "GAST"], "Organism": ["Homo sapiens"] * 2,
        "Structure_ID": ["2mjw", "AF-P01350-F1"], "Similarity": [1.0, 0.8],
        "Protein": ["2mjw_A", "AF-P01350-F1"], "Chain": ["A", "A"],
        "Pocket": ["P_1", "P_1"],
        "Pocket_Path": ["04_pockets/P25815/P25815_2mjw_A.pdb_res_P_1.pdb",
                        "04_pockets/P01350/AF-P01350-F1.pdb_res_P_1.pdb"],
        "Receptor_Path": ["03_pdbfixer/P25815/2mjw_A.pdb",
                          "02_structures/P01350/AF-P01350-F1.pdb"],
        "Is_Hub": [True, False]}).to_csv(layout.final_targets_csv, index=False)
    # a production library, for readout 3 (the control's score against the library's at the
    # same target). `production_docking` runs it; the bare fixture deliberately does not,
    # so the "readout 3 is unavailable" path is exercised too.
    layout.inputs.mkdir(parents=True, exist_ok=True)
    (layout.inputs / "e-Drug3D_2083.sdf").write_text("".join(
        f"MOL{i}\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        f"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n"
        f"> <ID>\n{i + 1}\n\n> <name>\ndrug{i}\n\n$$$$\n" for i in range(12)))
    return layout


def production_docking(layout):
    """Run `dock` on the fixture's tiny production library, so 07_docking exists."""
    from polyphar.stages import dock

    dock.run(layout, SimpleNamespace(
        nprocs=1, seed=42, replicates=1, threads_per_job=1, unit_timeout=600,
        engine="smina", search_mode=None, batch_size=None, poses=dock.POSES_FIRST,
        library=dock.DEFAULT_LIBRARY, dry_run=False, no_notify=True), Prog())


@pytest.fixture
def docked_run_folder(run_folder):
    production_docking(run_folder)
    return run_folder


def test_labels_only_the_experimental_target(run_folder):
    from polyphar.stages import dock

    layout = run_folder
    labels, unlabelled = control.label_targets(layout, dock.read_targets(layout))
    assert [l.target for l in labels] == ["2mjw_A_P_1"]
    assert labels[0].resname == "LIG" and labels[0].heavy_atoms == 16
    assert labels[0].contact == 1.0
    assert [t for t, _ in unlabelled] == ["AF-P01350-F1_P_1"]


@pytest.mark.openbabel
def test_the_stage_redocks_the_label_and_writes_the_file_stats_reads(run_folder):
    layout = run_folder
    counts = control.run(layout, _args(), Prog())
    d = control.stage_dir(layout)

    assert counts["labelled"] == 1 and counts["ligands"] == 1
    assert counts["panel_targets"] == 2 and counts["panel_targets_scored"] == 2
    assert counts["self_ranked"] == 1
    # one best pose per (target, ligand): the whole panel, not only the home target
    assert counts["poses"] == 2

    # the docking went to 08_control/docking, NOT to 07_docking
    assert (d / "docking" / "manifest.json").is_file()
    assert not layout.docking.exists()
    assert (layout.run_dir / control.LIBRARY).is_file()

    installed = layout.run_dir / control.INPUTS_POSITIVE_CONTROL
    result = json.loads(installed.read_text())
    assert result == json.loads((d / control.POSITIVE_CONTROL).read_text())
    assert result["n_redocked"] == 1
    assert result["source"].startswith("polyphar control")
    assert result["caveat"] == control.CAVEAT
    # the home target scores best (FAKE_POSE_BEST_RECEPTOR=2mjw), so the self-rank is the
    # best available on a two-target panel
    assert result["median_self_rank_percentile"] == pytest.approx(0.25)
    assert result["pose_in_pocket_fraction"] == 1.0
    # the stub echoes the ligand's own coordinates, so the self pair is recovered exactly
    assert result["redock_rmsd_success_fraction"] == 1.0
    assert result["n_redock_rmsd"] == 1 and result["redock_rmsd"]["median_rmsd"] == 0.0
    # the whole panel is scored, and both numbers are in the JSON, not only the manifest
    assert result["panel_targets_scored"] == 2 and result["panel_targets_total"] == 2
    assert result["panel_scored_fraction"] == 1.0
    # ... and `stats` reads it, which is the loop this stage closes. On ONE ligand it is
    # `unknown`, not `no` and certainly not `yes`: n = 1 is not evidence either way.
    verdict, detail = stats.positive_control_verdict(installed)
    assert verdict.startswith("unknown:")
    assert f"minimum of {stats.MIN_SELF_RANK_N}" in verdict
    assert "source" in detail
    assert (d / "summary.md").is_file() and control.CAVEAT in (d / "summary.md").read_text()
    assert (d / "manifest.json").is_file()    # written last


@pytest.mark.openbabel
def test_poses_outside_the_cavity_lower_the_fraction(run_folder, monkeypatch):
    layout = run_folder
    labels, _ = control.label_targets(layout, __import__(
        "polyphar.stages.dock", fromlist=["x"]).read_targets(layout))
    monkeypatch.setenv("FAKE_POSE_OUTSIDE_IDS", labels[0].ligand_id)
    counts = control.run(layout, _args(), Prog())
    assert counts["poses"] == 2
    assert counts["pose_in_pocket_fraction"] == 0.0
    result = json.loads((layout.run_dir
                         / control.INPUTS_POSITIVE_CONTROL).read_text())
    assert result["pose_in_pocket"]["n_in_pocket"] == 0
    # a pose 1000 A away is also a failed redocking
    assert result["redock_rmsd_success_fraction"] == 0.0


@pytest.mark.openbabel
def test_a_pose_beyond_two_angstrom_is_not_a_recovered_pose(run_folder, monkeypatch):
    """The gate, exercised on both sides. Pose-in-pocket cannot tell these apart -- a
    displacement of 1.9 A and one of 2.1 A are both well inside the cavity -- which is
    why it was not the right readout to gate on."""
    layout, run = run_folder, run_folder.run_dir
    import shutil as sh

    monkeypatch.setenv("FAKE_POSE_RMSD", "1.9")
    control.run(layout, _args(), Prog())
    near = json.loads((run / control.INPUTS_POSITIVE_CONTROL).read_text())
    assert near["redock_rmsd_success_fraction"] == 1.0
    assert near["redock_rmsd"]["median_rmsd"] == pytest.approx(1.9)
    assert near["pose_in_pocket_fraction"] == 1.0

    sh.rmtree(control.stage_dir(layout))
    (run / control.LIBRARY).unlink()
    monkeypatch.setenv("FAKE_POSE_RMSD", "2.1")
    control.run(layout, _args(), Prog())
    far = json.loads((run / control.INPUTS_POSITIVE_CONTROL).read_text())
    assert far["redock_rmsd_success_fraction"] == 0.0
    assert far["redock_rmsd"]["median_rmsd"] == pytest.approx(2.1)
    # ... and the descriptive readout is blind to the difference
    assert far["pose_in_pocket_fraction"] == near["pose_in_pocket_fraction"] == 1.0


@pytest.mark.openbabel
def test_readout_three_compares_against_the_production_library_at_the_same_target(
        docked_run_folder, monkeypatch):
    """`M` is within-target across ligands; the self-rank is its transpose. The third
    readout is the axis `M` is built on, and it is free: same receptor, same box, same
    engine, only the ligand set differs."""
    layout = docked_run_folder
    labels, _ = control.label_targets(layout, __import__(
        "polyphar.stages.dock", fromlist=["x"]).read_targets(layout))
    monkeypatch.setenv("FAKE_POSE_STRONG_IDS", labels[0].ligand_id)
    counts = control.run(layout, _args(), Prog())
    assert counts["library_percentiles"] == 1
    result = json.loads((layout.run_dir
                         / control.INPUTS_POSITIVE_CONTROL).read_text())
    lp = result["library_percentile"]
    assert lp["n_ligands"] == 1 and lp["n_top_decile"] == 1
    assert result["library_top_decile_fraction"] == 1.0
    rows = list(csv.DictReader(
        (control.stage_dir(layout) / "library_percentile.csv").open()))
    assert rows[0]["Home_Target"] == "2mjw_A_P_1"
    assert int(rows[0]["Library_Ligands"]) == 12 and rows[0]["Top_Decile"] == "True"


@pytest.mark.openbabel
def test_a_control_ligand_the_library_beats_at_home_fails_readout_three(docked_run_folder):
    """Without the bonus the control ligand is mid-pack among the library at its own
    pocket, and that is a failing readout 3 -- the number `M` needs."""
    layout = docked_run_folder
    control.run(layout, _args(), Prog())
    result = json.loads((layout.run_dir
                         / control.INPUTS_POSITIVE_CONTROL).read_text())
    assert result["library_top_decile_fraction"] == 0.0
    assert result["library_percentile"]["median_percentile"] > control.TOP_DECILE


@pytest.mark.openbabel
def test_readout_three_is_unavailable_without_a_production_docking(run_folder):
    """No 07_docking: the readout is not silently 0 of 0 reported as a pass. n is 0, and
    `stats` says `unknown`."""
    layout = run_folder
    counts = control.run(layout, _args(), Prog())
    assert counts["library_percentiles"] == 0
    installed = layout.run_dir / control.INPUTS_POSITIVE_CONTROL
    assert json.loads(installed.read_text())["n_library_percentile"] == 0
    verdict, _ = stats.positive_control_verdict(installed)
    assert verdict.startswith("unknown:") and "library percentile" in verdict
    failed = (control.stage_dir(layout) / "failed.tsv").read_text()
    assert "df_results" in failed


@pytest.mark.openbabel
def test_resume_reports_zero_units_and_rewrites_nothing(run_folder):
    layout = run_folder
    args = _args()
    assert control.units(layout, args).count > 0
    control.run(layout, args, Prog())
    assert control.units(layout, args).count == 0

    d = control.stage_dir(layout)
    before = {p.name: p.read_bytes() for p in d.glob("*") if p.is_file()}
    prog = Prog()
    counts = control.run(layout, args, prog)
    assert prog.calls == []                   # nothing was docked again
    assert counts["self_ranked"] == 1         # the finished counts come back
    assert {p.name: p.read_bytes() for p in d.glob("*") if p.is_file()} == before


@pytest.mark.openbabel
def test_nprocs_does_not_change_the_result(run_folder, tmp_path):
    import os

    layout = run_folder
    control.run(layout, _args(nprocs=1), Prog())
    one = (layout.run_dir / control.INPUTS_POSITIVE_CONTROL).read_text()

    other = tmp_path / "run4"
    shutil.copytree(layout.run_dir, other)
    for stale in (control.stage_dir(RunLayout(other)),):
        shutil.rmtree(stale)
    os.remove(other / control.LIBRARY)
    control.run(RunLayout(other), _args(nprocs=4), Prog())
    assert (other / control.INPUTS_POSITIVE_CONTROL).read_text() == one


@pytest.mark.openbabel
def test_a_panel_with_no_co_crystal_ligand_writes_no_passing_control(run_folder):
    """A control that cannot be produced must not leave behind something `stats` can read
    as a pass: no JSON at all, so stats says "unknown: not run", which is the truth."""
    layout = run_folder
    (layout.run_dir / "02_structures/P25815/2mjw_A.pdb").write_text(
        RECEPTOR.read_text())                 # the planted ligand is gone
    counts = control.run(layout, _args(), Prog())
    assert counts["labelled"] == 0 and counts["failed"] == 1
    installed = layout.run_dir / control.INPUTS_POSITIVE_CONTROL
    assert not installed.exists()
    verdict, _ = stats.positive_control_verdict(installed)
    assert verdict.startswith("unknown:")
    failed = (control.stage_dir(layout) / "failed.tsv").read_text()
    assert "not a control that passed" in failed


@pytest.mark.openbabel
def test_a_stale_passing_control_is_removed_when_the_control_cannot_be_produced(run_folder):
    layout = run_folder
    control.run(layout, _args(), Prog())
    installed = layout.run_dir / control.INPUTS_POSITIVE_CONTROL
    assert installed.is_file()
    shutil.rmtree(control.stage_dir(layout))
    (layout.run_dir / control.LIBRARY).unlink()
    (layout.run_dir / "02_structures/P25815/2mjw_A.pdb").write_text(RECEPTOR.read_text())
    control.run(layout, _args(), Prog())
    assert not installed.exists()


# --- wiring -------------------------------------------------------------------------------

def test_control_is_a_stage_between_dock_and_stats():
    """It must be in ORDER, or `polyphar run --to stats` would skip the control the stats
    stage then reports as "not run"."""
    from polyphar.benchmark import ORDER, STAGE_DIRS

    assert ORDER.index("dock") < ORDER.index("control") < ORDER.index("stats")
    assert ORDER[-1] == "stats"
    assert STAGE_DIRS["control"] == control.DIRNAME


def test_control_is_a_cli_subcommand_with_docks_engine_options():
    from polyphar.cli import build_parser

    args = build_parser().parse_args(["control", "--run-dir", "/tmp/x", "--engine", "smina",
                                      "--control-replicates", "2", "--seed", "7"])
    assert args.control_replicates == 2 and args.seed == 7
    assert control.dock_args(args).replicates == 2


def test_run_merges_dock_and_control_options_without_a_clash():
    """`polyphar run` refuses two stages that define one option differently; control shares
    dock's engine options by calling the same code, so the merge must succeed."""
    from polyphar.cli import build_parser

    args = build_parser().parse_args(["run", "--run-dir", "/tmp/x", "--to", "stats"])
    assert args.replicates == 3 and args.control_replicates == 1

def test_the_control_estimate_uses_docks_calibrated_rate():
    from polyphar.eta import CALIBRATION_ALIAS

    assert CALIBRATION_ALIAS["control"] == "dock"


def test_the_control_forces_replicate_zero_poses():
    """Readout 2 reads replicate 0's poses, so the pose policy is not the caller's to
    lower."""
    from polyphar.stages import dock

    assert control.dock_args(_args(poses="none")).poses == dock.POSES_FIRST


# --- P2: a non-finite score can never produce a passing rank -------------------------------

def test_a_non_finite_score_is_refused_not_ranked_as_the_best():
    """NaN compares false against everything: `better` counts 0 and `ties` counts only
    itself, so a NaN HOME score used to come out at percentile 0.0 -- better than the best
    achievable 0.5/N -- and an all-NaN matrix read as a perfect control, verdict `yes`."""
    nan = float("nan")
    scores = {f"T{i}": float(i) for i in range(100)}
    scores["T0"] = nan
    with pytest.raises(ValueError, match="non-finite"):
        control.self_rank_percentile("T0", scores)
    with pytest.raises(ValueError, match="non-finite"):
        control.self_rank_percentile("T0", {f"T{i}": nan for i in range(100)})
    # a non-finite COMPETITOR is refused too: it corrupts the ranking around it
    with pytest.raises(ValueError, match="non-finite"):
        control.self_rank_percentile("T0", {"T0": -9.0, "T1": nan, "T2": -1.0})


def test_the_control_reads_df_results_through_the_one_reader(run_folder):
    """There is one df_results reader in this repo and it already rejects non-finite
    scores, saying so in its docstring. The control used to do its own pd.read_csv +
    float()."""
    from polyphar.stages import stats

    layout = run_folder
    ctrl = control.ControlLayout(layout.run_dir)
    ctrl.docking.mkdir(parents=True, exist_ok=True)
    (ctrl.docking / "df_results_2mjw_A_P_1.csv").write_text(
        "ID,score_mean\nbad,nan\ngood,-7.0\n")
    (ctrl.docking / "df_results_AF-P01350-F1_P_1.csv").write_text(
        "ID,score_mean\nbad,-1.0\ngood,-1.0\n")
    from polyphar.stages import dock

    scores, reasons = control._summary_scores(layout, dock.read_targets(layout))
    assert scores["2mjw_A_P_1"] == {"good": -7.0}        # the NaN row is gone
    assert any("not finite" in why for _item, why in reasons)
    # ... and the ligand that lost its score gets no rank at all, rather than a short one
    labels = [control.Label("2mjw_A_P_1", "P25815", "p", "bad", "LIG", "A", "1", 20, 0.9,
                            Path("x"), ("LIG", "A", "1"))]
    rows, why = control.self_ranks(labels, scores)
    assert rows == [] and "no self-rank" in why[0][1]
    assert stats.read_df_results is not None   # the reader this stage now shares


# --- P4: the denominator is stated, not implied --------------------------------------------

def test_the_source_sentence_does_not_say_all_of_a_shrunken_panel():
    """"redocked against all 382 panel targets" is the sentence a reviewer quotes. When
    100 of 482 failed, it has to say so."""
    args = _args()
    dargs = control.dock_args(args)
    layout = RunLayout(Path("/tmp/nonexistent-run"))
    whole = control._source(layout, args, dargs, 24, 482, 482)
    part = control._source(layout, args, dargs, 24, 382, 482)
    assert "against all 482 of the panel's 482 targets" in whole
    assert "against 382 of the panel's 482 targets" in part
    assert "all 382" not in part


def test_the_json_carries_both_the_scored_and_the_total_target_count():
    res = _result(n_panel=382, n_total=482)
    assert res["panel_targets_scored"] == 382 and res["panel_targets_total"] == 482
    assert res["panel_scored_fraction"] == pytest.approx(382 / 482)
    # the worst case if every unscored target had beaten the home pocket
    assert res["self_rank"]["median_percentile_adversarial"] == pytest.approx(
        (0.004 * 382 + 100) / 482)


def test_a_half_scored_panel_is_unknown_not_a_pass(tmp_path):
    p = tmp_path / "pc.json"
    p.write_text(json.dumps(_result(n_panel=382, n_total=482)))
    verdict, _ = stats.positive_control_verdict(p)
    assert verdict.startswith("unknown:") and "382" in verdict and "482" in verdict
    p.write_text(json.dumps(_result(n_panel=482, n_total=482)))
    assert stats.positive_control_verdict(p)[0] == "yes"


# --- P8: one physical ligand, one readout --------------------------------------------------

def test_a_ligand_filling_a_pocket_and_its_subpocket_is_counted_once(run_folder):
    """DoGSiteScorer reports nested subpockets, and a ligand in a parent pocket usually
    also fills its subpocket. Keyed on the TARGET, that one molecule produced two ligand
    IDs, two SDF records and two self-ranks -- n inflated with nothing added."""
    import pandas as pd
    from polyphar.stages import dock

    layout = run_folder
    sub = "04_pockets/P25815/P25815_2mjw_A.pdb_res_P_1_1.pdb"
    shutil.copyfile(layout.run_dir / "04_pockets/P25815/P25815_2mjw_A.pdb_res_P_1.pdb",
                    layout.run_dir / sub)
    df = pd.read_csv(layout.final_targets_csv)
    extra = df.iloc[[0]].copy()
    extra["Pocket"] = "P_1_1"
    extra["Pocket_Path"] = sub
    pd.concat([df, extra]).to_csv(layout.final_targets_csv, index=False)

    targets = dock.read_targets(layout)
    assert len(targets) == 3
    labels, unlabelled = control.label_targets(layout, targets)
    assert len(labels) == 1, [l.ligand_id for l in labels]
    assert {t for t, _ in unlabelled} == {"AF-P01350-F1_P_1", "2mjw_A_P_1_1"}
    assert any("counted once" in why for _t, why in unlabelled)


def test_an_unreadable_pose_counts_against_the_readout_not_out_of_it(tmp_path):
    """A V3000 record parses to 0 atoms. Dropping it from BOTH the numerator and the
    denominator is the same shape of error as excluding a failed target from everyone's
    denominator."""
    rec = ("EDRUG3D_L1\n t\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n"
           "M  V30 BEGIN CTAB\nM  END\n> <minimizedAffinity>\n-9.0\n\n$$$$\n")
    sdf = tmp_path / "poses.sdf"
    sdf.write_text(rec)
    best = control.best_poses(sdf, "minimizedAffinity")
    assert best["L1"][1] == [] and "V3000" in best["L1"][2]
    # and the JSON counts it: 1 pose, 0 in pocket, 0 redocking successes
    rows = [{"Target": "t", "Ligand_ID": "L1", "Score": -9.0, "In_Pocket": False,
             "Self_Pair": True, "RMSD": "", "Redock_Success": False, "Reason": "V3000"}]
    res = control.positive_control_json(_rank_rows(12, 0.01), rows, [], 10,
                                        "t", lib_rows=_lib_rows(12), n_panel_total=10)
    assert res["pose_in_pocket"]["n_poses"] == 1
    assert res["pose_in_pocket_fraction"] == 0.0
    assert res["n_redock_rmsd"] == 1 and res["redock_rmsd_success_fraction"] == 0.0
    assert res["redock_rmsd"]["n_unmeasurable"] == 1


# --- the RMSD itself -----------------------------------------------------------------------

def test_the_rmsd_is_symmetry_insensitive_and_needs_no_superposition():
    ref = [(0.0, 0.0, 0.0, "C"), (1.4, 0.0, 0.0, "C"), (0.0, 1.4, 0.0, "N")]
    assert control.symmetry_insensitive_rmsd(list(reversed(ref)), ref) == 0.0
    # two interchangeable carbons swapped: order-preserving RMSD would be 1.4, this is 0
    swapped = [ref[1], ref[0], ref[2]]
    assert control.symmetry_insensitive_rmsd(swapped, ref) == 0.0
    # ... but an element is never matched to a different element
    wrong = [(0.0, 1.4, 0.0, "C"), (1.4, 0.0, 0.0, "C"), (0.0, 0.0, 0.0, "N")]
    assert control.symmetry_insensitive_rmsd(wrong, ref) > 1.0
    # a rigid 3 A translation is 3 A: nothing is superposed away
    moved = [(x + 3.0, y, z, e) for x, y, z, e in ref]
    assert control.symmetry_insensitive_rmsd(moved, ref) == pytest.approx(3.0)
    with pytest.raises(ValueError, match="composition differs"):
        control.symmetry_insensitive_rmsd(ref[:2], ref)


# --- P7: wiring ------------------------------------------------------------------------------

def test_the_positive_control_out_path_cannot_overwrite_another_stages_file():
    """It is a path the stage overwrites and, on a failing control, UNLINKS. Unvalidated,
    `--positive-control-out 07_docking/manifest.json` deleted the docking manifest."""
    import argparse

    for bad in ("07_docking/manifest.json", "../escape.json", "/etc/passwd",
                "inputs/../../x.json", "inputs/notes.txt", "08_control/manifest.json",
                "02_structures/P1/x.json"):
        with pytest.raises(argparse.ArgumentTypeError):
            control.out_path(bad)
    for good in ("inputs/positive_control.json", "08_control/positive_control.json",
                 "inputs/sub/pc.json"):
        assert control.out_path(good) == good


def test_run_refuses_a_bad_positive_control_out_before_touching_anything(run_folder):
    import argparse

    layout = run_folder
    victim = layout.run_dir / "07_docking" / "manifest.json"
    victim.parent.mkdir(parents=True, exist_ok=True)
    victim.write_text("{}")
    with pytest.raises(argparse.ArgumentTypeError):
        control.run(layout, _args(positive_control_out="07_docking/manifest.json"), Prog())
    assert victim.is_file()


def test_the_cli_rejects_a_bad_positive_control_out():
    from polyphar.cli import build_parser

    with pytest.raises(SystemExit):
        build_parser().parse_args(["control", "--run-dir", "/tmp/x",
                                   "--positive-control-out", "07_docking/manifest.json"])


def test_the_installed_path_is_the_one_stats_reads_not_a_second_literal():
    assert control.INPUTS_POSITIVE_CONTROL == stats.POSITIVE_CONTROL


@pytest.mark.openbabel
def test_the_control_reports_its_own_docking_failures(run_folder):
    """Folded in only as `dock_*`, a failed control docking was invisible to
    cli.flagged(): the stage exited 0 while `dock` exits 3 for the same condition."""
    from polyphar import cli

    layout = run_folder
    # a receptor that cannot be prepared: the target's own file is gone
    (layout.run_dir / "03_pdbfixer/P25815/2mjw_A.pdb").write_text("not a pdb\n")
    counts = control.run(layout, _args(), Prog())
    # dock reports the broken receptor as an incomplete target; whatever it calls it,
    # the count has to reach a FLAG_KEY at the top level of the control's own counts.
    assert counts["dock_incomplete_targets"] > 0
    assert counts["incomplete_targets"] >= counts["dock_incomplete_targets"]
    assert counts["failed"] >= counts.get("dock_failed", 0)
    assert cli.flagged(counts), counts
    assert any(k in cli.FLAG_KEYS for k in counts if counts.get(k))


def test_the_control_estimate_uses_docks_engine_and_docks_slot_count(monkeypatch, tmp_path):
    """control IS dock: same engine, same units. Left out of the per-engine branch it fell
    through to `calib.cost("dock", None)` -- which accepts a dock entry measured with ANY
    engine -- and to the engine-blind parallel-efficiency table, so it printed
    "calibrated, 8 slots" for a Uni-Dock run the dock row correctly called uncalibrated,
    and halved it again at --nprocs 16."""
    from polyphar import benchmark
    from polyphar.eta import Calibration

    monkeypatch.setattr(benchmark, "bench_engine", lambda: "unidock")
    calib = Calibration(tmp_path / "calibration.json",
                        {"dock": {"seconds_per_unit": 1.0, "tool": "smina",
                                  "n_observations": 1}})
    units = {"dock": 100, "control": 100}
    seen = {}

    def fake_estimate(stage, u, nprocs, c, pe, *, threads_per_job=1, tool=None, slots=None):
        seen[stage] = (tool, slots)
        from polyphar.eta import Estimate

        return Estimate(None, "uncalibrated")

    monkeypatch.setattr("polyphar.eta.estimate_stage", fake_estimate)
    from polyphar.common import WorkUnits

    monkeypatch.setattr(benchmark, "projected_units",
                        lambda stage, s: WorkUnits(units.get(stage, 1), "ligands"))
    monkeypatch.setattr(benchmark, "completed", lambda layout: list(benchmark.ORDER[:-3]))
    rows, _obs = benchmark.estimate_rows(RunLayout(tmp_path), 16, calib)
    assert seen["control"] == seen["dock"], seen
    assert seen["control"][0] is not None       # the engine, not None


def test_the_recalibrate_hint_names_a_stage_benchmark_accepts():
    from polyphar.eta import CALIBRATION_ALIAS, RECALIBRATE_HINT
    from polyphar.benchmark import BENCH_STAGES

    for stage, alias in CALIBRATION_ALIAS.items():
        hint = RECALIBRATE_HINT.get(stage, "polyphar benchmark --stages {stage}")
        named = hint.format(tool="smina", stage=stage).rsplit(" ", 1)[-1]
        assert named in BENCH_STAGES, f"{stage}: `{hint}` names {named!r}"


def test_the_control_rate_is_written_where_the_estimate_reads_it(tmp_path, monkeypatch):
    """`estimate_stage` reads CALIBRATION_ALIAS[stage]; Progress.done used to write under
    the stage's own name, so control's measured rate was stored and never read."""
    from polyphar.eta import CALIBRATION_ALIAS, Calibration, Estimate, Progress

    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path))
    calib = Calibration(tmp_path / "calibration.json", {})
    prog = Progress("control", 10, Estimate(1.0, "x"), lambda *a, **k: None,
                    calib=calib, slots=1.0, log_path=tmp_path / "l.log")
    prog.advance(10)
    prog.done({})
    assert calib.cost(CALIBRATION_ALIAS["control"], None) is not None
    assert calib.cost("control", None) is None


def test_the_control_unit_basis_matches_docks():
    """control borrows dock's calibration entry, and the basis is printed next to it."""
    from polyphar.benchmark import PRIORS, projected_units

    assert (projected_units("control", PRIORS).cost_basis
            == projected_units("dock", PRIORS).cost_basis == "ligands")


def test_the_gpu_stop_note_does_not_advise_a_tool_that_does_not_exist():
    """`_finish` puts the note FIRST in the Pushover message, so a control-specific stop
    followed by "run dock-export / dock-worker / dock-merge" contradicted itself."""
    import inspect

    from polyphar.stages import run as run_stage

    src = inspect.getsource(run_stage._run_stages)
    assert 'if name == "dock" else' in src
    assert "dock-export / dock-worker / dock-merge" in src
    # the advice is dock's only on dock's branch; control gets its own sentence
    note_block = src[src.index("note = ("):src.index("_finish(layout, rows, notifier,\n"
                                                     "                f\"run stopped")]
    assert "--engine smina" in note_block
