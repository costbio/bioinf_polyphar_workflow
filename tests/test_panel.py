"""Tests for the panel filters: each one rejects what it should, the funnel is a chain
that sums, the TM policies differ as measured, and the one-row-per-protein rule is
order-independent."""
import pytest

from polyphar import panel
from polyphar.stages.annotate import encode_residues


def _row(acc="P1", pocket="P_1", fname=None, ds=0.9, vol=800.0):
    return {"UniProt ID": acc, "Gene Name": "G_HUMAN", "Protein": "AF-P1-F1",
            "Pocket": pocket, "Pocket Filename": fname or f"AF-P1-F1_res_{pocket}.pdb",
            "drugScore": str(ds), "volume": str(vol)}


def _ann(af=True, rg=6.0, plddt=90.0, size=(18.0, 18.0, 18.0), residues=(10, 11, 12)):
    return {"AlphaFold": str(af), "radius_of_gyration": str(rg),
            "pocket_plddt_mean": "" if plddt is None else str(plddt),
            "box_size_x": str(size[0]), "box_size_y": str(size[1]),
            "box_size_z": str(size[2]),
            "residues": encode_residues(list(residues)),
            "_residues": list(residues)}


OPTS = panel.Options()


@pytest.mark.parametrize("kw,ann_kw,expect", [
    ({}, {}, None),
    ({"ds": 0.4}, {}, "drugscore"),
    ({"vol": 250.0}, {}, "volume"),
    ({"vol": 2500.0}, {}, "volume"),
    ({}, {"plddt": 60.0}, "plddt"),
    ({}, {"plddt": None}, "plddt"),
    ({}, {"rg": 9.0}, "rg"),
    ({}, {"size": (50.0, 10.0, 10.0)}, "grid"),      # 134 points on one axis
])
def test_each_filter_rejects_at_its_own_step(kw, ann_kw, expect):
    assert panel.classify(_row(**kw), _ann(**ann_kw), [], None, OPTS) == expect


def test_plddt_floor_applies_only_to_alphafold():
    """An experimental pocket has no pLDDT; requiring one would reject every one of them."""
    assert panel.classify(_row(), _ann(af=False, plddt=None), [], None, OPTS) is None


def test_subpocket_of_an_oversized_parent_is_rejected():
    r = _row(pocket="P_1_1")
    assert panel.classify(r, _ann(), [], 900.0, OPTS) is None
    assert panel.classify(r, _ann(), [], 25840.0, OPTS) == "subpocket_parent"


def test_rg_is_checked_before_the_grid_guard():
    """Rg bounds the box, so a pocket failing both must be attributed to Rg — otherwise
    the funnel would credit the grid guard with work Rg already did."""
    assert panel.classify(_row(), _ann(rg=20.0, size=(60.0, 60.0, 60.0)),
                          [], None, OPTS) == "rg"


class TestTMPolicies:
    SOLUBLE, SINGLE, MULTI = [], [(5, 25)], [(5, 25), (40, 60)]

    def test_policy_a_drops_only_multi_pass(self):
        o = panel.Options(tm_policy="A")
        assert not panel.tm_reject(True, self.SINGLE, [10, 11], o)
        assert panel.tm_reject(True, self.MULTI, [10, 11], o)

    def test_policy_b_drops_any_tm_protein(self):
        o = panel.Options(tm_policy="B")
        assert panel.tm_reject(True, self.SINGLE, [100, 101], o)
        assert not panel.tm_reject(True, self.SOLUBLE, [100, 101], o)

    def test_policy_c_drops_a_pocket_sitting_in_the_membrane(self):
        o = panel.Options(tm_policy="C")
        # single-pass protein, pocket outside the TM helix -> kept
        assert not panel.tm_reject(True, self.SINGLE, [100, 101, 102], o)
        # same protein, pocket inside it -> dropped
        assert panel.tm_reject(True, self.SINGLE, [10, 11, 12], o)

    def test_policy_c_overlap_check_is_alphafold_only(self):
        """Experimental chains carry pdbfixer numbering, so comparing them against
        UniProt TM positions would be meaningless."""
        o = panel.Options(tm_policy="C")
        assert panel.tm_reject(True, self.SINGLE, [10, 11, 12], o)
        assert not panel.tm_reject(False, self.SINGLE, [10, 11, 12], o)

    def test_policy_c_threshold_is_a_fraction_of_the_pocket(self):
        o = panel.Options(tm_policy="C", tm_overlap_max=0.10)
        residues = list(range(100, 120))              # 20 residues
        assert not panel.tm_reject(True, [(100, 101)], residues, o)   # 2/20 = 10%, not >
        assert panel.tm_reject(True, [(100, 103)], residues, o)       # 4/20 = 20%


def test_funnel_is_a_chain_that_sums():
    rows = [_row(acc="A", ds=0.4), _row(acc="B", vol=99.0), _row(acc="C")]
    anns = {r["Pocket Filename"]: _ann() for r in rows}
    for r in rows:
        anns[r["Pocket Filename"]] = _ann()
    pan, funnel, rej = panel.select_panel(rows, anns, {}, {}, OPTS)
    assert funnel[0] == ("selected", 3)
    assert [s for s, _ in funnel][1:] == panel.STEPS
    counts = [n for _, n in funnel]
    assert counts == sorted(counts, reverse=True)      # monotone
    assert counts[-1] == 1                             # only C survives
    assert len(rej) == 2


def test_one_row_per_protein_takes_the_highest_drugscore():
    rows = [_row(acc="P1", pocket="P_1", ds=0.6), _row(acc="P1", pocket="P_2", ds=0.95)]
    anns = {r["Pocket Filename"]: _ann() for r in rows}
    pan, _, _ = panel.select_panel(rows, anns, {}, {}, OPTS)
    assert len(pan) == 1
    assert pan["P1"][0]["Pocket"] == "P_2"


def test_tie_on_drugscore_breaks_on_pocket_name_not_row_order():
    a = _row(acc="P1", pocket="P_1", ds=0.8)
    b = _row(acc="P1", pocket="P_2", ds=0.8)
    anns = {r["Pocket Filename"]: _ann() for r in (a, b)}
    first, _, _ = panel.select_panel([a, b], anns, {}, {}, OPTS)
    second, _, _ = panel.select_panel([b, a], anns, {}, {}, OPTS)
    assert first["P1"][0]["Pocket"] == second["P1"][0]["Pocket"] == "P_1"


def test_a_pocket_without_an_annotation_is_reported_not_silently_dropped():
    rows = [_row(acc="P1")]
    pan, funnel, rej = panel.select_panel(rows, {}, {}, {}, OPTS)
    assert pan == {}
    assert rej == [("P1/AF-P1-F1_res_P_1.pdb", "missing annotation")]


def test_options_rejects_an_unknown_tm_policy():
    with pytest.raises(ValueError, match="tm_policy"):
        panel.Options(tm_policy="Z")


def test_parent_key_includes_the_structure():
    """`P_1` exists for every structure of an accession, so keying the parent lookup on
    (accession, pocket) alone silently returns another structure's pocket volume."""
    assert panel.parent_key({"UniProt ID": "P1", "Protein": "AF-P1-F1",
                             "Pocket": "P_1"}) is None
    assert panel.parent_key({"UniProt ID": "P1", "Protein": "AF-P1-F1",
                             "Pocket": "P_3_2"}) == ("P1", "AF-P1-F1", "P_3")


def test_subpocket_parent_lookup_does_not_cross_structures():
    """Two structures of one accession both have a P_1. A subpocket of the small one
    must not be judged against the large one's parent."""
    rows = [{"UniProt ID": "P1", "Gene Name": "G_HUMAN", "Protein": "1ABC_A",
             "Pocket": "P_1_1", "Pocket Filename": "1abc_A_res_P_1_1.pdb",
             "drugScore": "0.9", "volume": "800"}]
    anns = {"1abc_A_res_P_1_1.pdb": _ann(af=False, plddt=None)}
    volumes = {("P1", "1ABC_A", "P_1"): 900.0,       # its real parent: fine
               ("P1", "AF-P1-F1", "P_1"): 25840.0}   # another structure's: irrelevant
    pan, _, rej = panel.select_panel(rows, anns, {}, volumes, OPTS)
    assert len(pan) == 1 and not rej
