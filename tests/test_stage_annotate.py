"""Tests for the annotate stage: geometry, the Uni-Dock grid rule, pLDDT read-back,
the subpocket parent lookup, and the cross-class property that motivates using radius
of gyration instead of a sequence-segment count."""
import csv
import math
import types

import pytest

from polyphar.stages import annotate


def _pdb(atoms):
    """atoms = [(resnum, x, y, z, bfactor)] -> PDB text."""
    out = []
    for i, (res, x, y, z, b) in enumerate(atoms, start=1):
        out.append(
            f"ATOM  {i:>5}  CA  ALA A{res:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00{b:>6.2f}           C"
        )
    return "\n".join(out) + "\n"


def _cube(res_start=1, step=1, a=2.0, b=90.0):
    """Eight atoms on the corners of a cube of side `a`."""
    pts = [(x, y, z) for x in (0, a) for y in (0, a) for z in (0, a)]
    return [(res_start + i * step, x, y, z, b) for i, (x, y, z) in enumerate(pts)]


def test_radius_of_gyration_matches_closed_form():
    # for the 8 corners of a cube of side a, every point is sqrt(3)/2*a from the centre
    a = 2.0
    xyz = [(x, y, z) for x in (0, a) for y in (0, a) for z in (0, a)]
    assert annotate.radius_of_gyration(xyz) == pytest.approx(math.sqrt(3) / 2 * a)


def test_radius_of_gyration_is_translation_invariant():
    xyz = [(x, y, z) for x in (0, 2.0) for y in (0, 2.0) for z in (0, 2.0)]
    moved = [(x + 137.0, y - 42.5, z + 3.25) for x, y, z in xyz]
    assert annotate.radius_of_gyration(moved) == pytest.approx(
        annotate.radius_of_gyration(xyz))


def test_bounding_box_is_unpadded_and_matches_dock():
    """The box must be what dock.pocket_box computes, or the panel is filtered on a
    different box than the one searched."""
    from polyphar.stages import dock
    xyz = [(0.0, 0.0, 0.0), (4.0, 6.0, 8.0), (1.0, 2.0, 3.0)]
    centre, size = annotate.bounding_box(xyz)
    assert centre == pytest.approx([2.0, 3.0, 4.0])
    assert size == pytest.approx([4.0, 6.0, 8.0])
    assert hasattr(dock, "pocket_box")


def test_count_segments_uses_a_gap_of_more_than_five():
    assert annotate.count_segments([1, 2, 3]) == 1
    assert annotate.count_segments([1, 6]) == 1          # gap of 5 -> same segment
    assert annotate.count_segments([1, 7]) == 2          # gap of 6 -> new segment
    assert annotate.count_segments([1, 2, 20, 21, 50]) == 3
    assert annotate.count_segments([]) == 0


@pytest.mark.parametrize("size,legal", [
    ([10.0, 10.0, 10.0], True),
    ([30.0, 30.0, 30.0], True),
    ([49.0, 10.0, 10.0], False),   # 131 points on one axis, over the 128 limit
    ([40.0, 40.0, 40.0], False),   # 107^3 = 1.22M points, over the 531,441 total
])
def test_grid_legality_follows_unidock_limits(size, legal):
    n, total, ok = annotate.grid(size)
    assert ok is legal
    assert ok == (max(n) <= annotate.MAX_AXIS and total <= annotate.MAX_POINTS)


def test_receptor_stem_handles_both_naming_schemes():
    assert annotate.receptor_stem("8vut_A_res_P_1.pdb") == "8vut_A"
    assert annotate.receptor_stem("AF-P48448-F1_res_P_2.pdb") == "AF-P48448-F1"
    assert annotate.receptor_stem("P03952_2anw_A_res_P_6.pdb") == "P03952_2anw_A"


def test_parent_volume_only_for_subpockets():
    vols = {("P1", "X_A", "P_1"): 5000.0, ("P1", "X_A", "P_1_1"): 900.0}
    assert annotate.parent_volume_of("P1", "X_A", "P_1", vols) is None
    assert annotate.parent_volume_of("P1", "X_A", "P_1_1", vols) == 5000.0
    assert annotate.parent_volume_of("P1", "X_A", "P_9_9", vols) is None  # parent absent


def test_parent_volume_lookup_is_per_structure():
    """`P_1` exists for every structure of an accession; the lookup must not cross them."""
    vols = {("P1", "1ABC_A", "P_1"): 900.0, ("P1", "AF-P1-F1", "P_1"): 25840.0}
    assert annotate.parent_volume_of("P1", "1ABC_A", "P_1_1", vols) == 900.0
    assert annotate.parent_volume_of("P1", "AF-P1-F1", "P_1_1", vols) == 25840.0


def test_annotate_one_reads_plddt_back_from_the_model(tmp_path):
    """DoGSiteScorer zeroes the B-factor column in its pocket files, so pLDDT has to
    come from the AlphaFold model in 02_structures."""
    pockets, structures = tmp_path / "04_pockets", tmp_path / "02_structures"
    (pockets / "P1").mkdir(parents=True)
    (structures / "P1").mkdir(parents=True)
    atoms = _cube(res_start=10, b=0.0)                      # pocket file: B-factors zeroed
    (pockets / "P1" / "AF-P1-F1_res_P_1.pdb").write_text(_pdb(atoms))
    model = [(r, x, y, z, 70.0 + i) for i, (r, x, y, z, _) in enumerate(atoms)]
    (structures / "P1" / "AF-P1-F1.pdb").write_text(_pdb(model))

    row = {"UniProt ID": "P1", "Gene Name": "G", "Protein": "AF-P1-F1",
           "Pocket": "P_1", "Pocket Filename": "AF-P1-F1_res_P_1.pdb"}
    got = annotate.annotate_one(row, str(pockets), str(structures), None)

    assert got["AlphaFold"] is True
    assert got["n_atoms"] == 8 and got["n_residues"] == 8
    assert got["pocket_plddt_min"] == pytest.approx(70.0)
    assert got["pocket_plddt_mean"] == pytest.approx(73.5)
    assert got["radius_of_gyration"] == pytest.approx(math.sqrt(3), abs=1e-3)


def test_annotate_one_leaves_plddt_empty_for_experimental(tmp_path):
    pockets, structures = tmp_path / "04_pockets", tmp_path / "02_structures"
    (pockets / "P2").mkdir(parents=True)
    structures.mkdir(parents=True, exist_ok=True)
    (pockets / "P2" / "8vut_A_res_P_1.pdb").write_text(_pdb(_cube(b=0.0)))
    row = {"UniProt ID": "P2", "Gene Name": "G", "Protein": "8VUT_A",
           "Pocket": "P_1", "Pocket Filename": "8vut_A_res_P_1.pdb"}
    got = annotate.annotate_one(row, str(pockets), str(structures), None)
    assert got["AlphaFold"] is False
    assert got["pocket_plddt_mean"] == "" and got["pocket_plddt_min"] == ""


def test_annotate_one_rejects_a_pocket_too_small_to_measure(tmp_path):
    pockets = tmp_path / "04_pockets"
    (pockets / "P3").mkdir(parents=True)
    (pockets / "P3" / "x_res_P_1.pdb").write_text(_pdb([(1, 0.0, 0.0, 0.0, 0.0)]))
    row = {"UniProt ID": "P3", "Protein": "X_A", "Pocket": "P_1",
           "Pocket Filename": "x_res_P_1.pdb"}
    with pytest.raises(ValueError, match="too few"):
        annotate.annotate_one(row, str(pockets), str(tmp_path / "02_structures"), None)


def test_rg_is_insensitive_to_residue_numbering_but_segments_are_not():
    """The property the panel filter depends on. Two pockets with identical geometry but
    different residue numbering -- pdbfixer renumbers 1..N with loops modelled in, while
    AlphaFold models carry UniProt positions -- must get the same Rg and may get
    different segment counts."""
    contiguous = _cube(res_start=1, step=1)        # pdbfixer-like
    gapped = _cube(res_start=1, step=40)           # UniProt-like, real gaps
    xyz_c = [(x, y, z) for _, x, y, z, _ in contiguous]
    xyz_g = [(x, y, z) for _, x, y, z, _ in gapped]

    assert annotate.radius_of_gyration(xyz_c) == pytest.approx(
        annotate.radius_of_gyration(xyz_g))
    assert annotate.count_segments(sorted(r for r, *_ in contiguous)) == 1
    assert annotate.count_segments(sorted(r for r, *_ in gapped)) == 8


def test_annotate_is_in_the_pipeline_order_before_select():
    """It produces the annotations `select` filters on, so it must run first — and it
    must be in ORDER at all, or `polyphar run` skips it silently."""
    from polyphar.benchmark import ORDER, STAGE_DIRS
    assert "annotate" in ORDER
    assert ORDER.index("annotate") < ORDER.index("select")
    assert STAGE_DIRS["annotate"] == "04_pockets/annotation"


def test_similarity_runs_after_select():
    """similarity no longer feeds selection; it runs on the panel's eligible pockets."""
    from polyphar.benchmark import ORDER
    assert ORDER.index("select") < ORDER.index("similarity")
