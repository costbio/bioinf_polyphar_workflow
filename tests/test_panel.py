"""Tests for the panel filters: each one rejects what it should, the funnel is a chain
that sums, the TM policies differ as measured, and the one-row-per-protein rule is
order-independent."""
import pytest

from polyphar import panel
from polyphar.stages.annotate import encode_residues


def _row(acc="P1", pocket="P_1", fname=None, ds=0.9, vol=800.0, sv=0.40):
    return {"UniProt ID": acc, "Gene Name": "G_HUMAN", "Protein": "AF-P1-F1",
            "Pocket": pocket, "Pocket Filename": fname or f"AF-P1-F1_res_{pocket}.pdb",
            "drugScore": str(ds), "volume": str(vol), "surf/vol": str(sv)}


def _ann(af=True, plddt=90.0, size=(18.0, 18.0, 18.0), residues=(10, 11, 12)):
    return {"AlphaFold": str(af),
            "pocket_plddt_mean": "" if plddt is None else str(plddt),
            "box_size_x": str(size[0]), "box_size_y": str(size[1]),
            "box_size_z": str(size[2]),
            "residues": encode_residues(list(residues)),
            "_residues": list(residues)}


OPTS = panel.Options()
ONE = panel.Options(per_protein=1)


@pytest.mark.parametrize("kw,ann_kw,expect", [
    ({}, {}, None),
    ({"ds": 0.4}, {}, "drugscore"),
    ({"vol": 250.0}, {}, "volume"),
    ({"vol": 2500.0}, {}, "volume"),
    ({}, {"plddt": 60.0}, "plddt"),
    ({}, {"plddt": None}, "plddt"),
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


def test_volume_is_checked_before_the_tm_rule():
    """Order matters for the funnel: a pocket failing both must be attributed to the
    first step, or the counts stop forming a chain."""
    assert panel.classify(_row(vol=9000.0), _ann(), [(5, 25)], None, OPTS) == "volume" 


def test_volume_ceiling_is_what_rejects_oversized_cavities():
    """The oversized interdomain cavities run ~2,600 A^3 against ~270 for a typical
    pocket. surf/vol could not separate them (0.539 against 0.510) and was removed."""
    oversized = _row(vol=2600.0)
    assert panel.classify(oversized, _ann(), [], None, OPTS) == "volume"


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
    anns = {(r["UniProt ID"], r["Pocket Filename"]): _ann() for r in rows}
    pan, elig, funnel, rej = panel.select_panel(rows, anns, {}, {}, OPTS)
    assert funnel[0] == ("selected", 3)
    assert [s for s, _ in funnel][1:] == panel.STEPS
    counts = [n for _, n in funnel]
    assert counts == sorted(counts, reverse=True)      # monotone
    assert counts[-1] == 1                             # only C survives
    assert len(rej) == 2


def test_one_row_per_protein_takes_the_highest_drugscore():
    rows = [_row(acc="P1", pocket="P_1", ds=0.6), _row(acc="P1", pocket="P_2", ds=0.95)]
    anns = {(r["UniProt ID"], r["Pocket Filename"]): _ann() for r in rows}
    pan, _, _, _ = panel.select_panel(rows, anns, {}, {}, panel.Options(per_protein=1))
    assert len(pan) == 1
    assert pan["P1"][0][0]["Pocket"] == "P_2"


def test_tie_on_drugscore_breaks_on_pocket_name_not_row_order():
    a = _row(acc="P1", pocket="P_1", ds=0.8)
    b = _row(acc="P1", pocket="P_2", ds=0.8)
    anns = {(r["UniProt ID"], r["Pocket Filename"]): _ann() for r in (a, b)}
    first, _, _, _ = panel.select_panel([a, b], anns, {}, {}, ONE)
    second, _, _, _ = panel.select_panel([b, a], anns, {}, {}, ONE)
    assert first["P1"][0][0]["Pocket"] == second["P1"][0][0]["Pocket"] == "P_1"


def test_a_pocket_without_an_annotation_is_reported_not_silently_dropped():
    rows = [_row(acc="P1")]
    pan, elig, funnel, rej = panel.select_panel(rows, {}, {}, {}, OPTS)
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
             "drugScore": "0.9", "volume": "800", "surf/vol": "0.40"}]
    anns = {("P1", "1abc_A_res_P_1_1.pdb"): _ann(af=False, plddt=None)}
    volumes = {("P1", "1ABC_A", "P_1"): 1200.0,      # its real parent: under the ceiling
               ("P1", "AF-P1-F1", "P_1"): 25840.0}   # another structure's: irrelevant
    pan, _, _, rej = panel.select_panel(rows, anns, {}, volumes, OPTS)
    assert len(pan) == 1 and not rej


def test_eligible_returns_every_passing_pocket_not_one_per_protein():
    """`similarity` runs on the eligible set, so it must not be collapsed to the panel."""
    rows = [_row(acc="P1", pocket="P_1", ds=0.9), _row(acc="P1", pocket="P_2", ds=0.8),
            _row(acc="P2", pocket="P_1", ds=0.7)]
    anns = {(r["UniProt ID"], r["Pocket Filename"]): _ann() for r in rows}
    pan, elig, _, _ = panel.select_panel(rows, anns, {}, {}, OPTS)
    assert len(pan) == 2 and len(elig) == 3


def test_up_to_n_pockets_per_protein_best_drugscore_first():
    """Picking one pocket per protein was close to a coin flip against known ligand
    sites (17 of 29), so the best few are docked instead."""
    rows = [_row(acc="P1", pocket="P_1", ds=0.70), _row(acc="P1", pocket="P_2", ds=0.95),
            _row(acc="P1", pocket="P_3", ds=0.80)]
    anns = {(r["UniProt ID"], r["Pocket Filename"]): _ann() for r in rows}
    pan, _, _, _ = panel.select_panel(rows, anns, {}, {}, panel.Options(per_protein=2))
    assert [p[0]["Pocket"] for p in pan["P1"]] == ["P_2", "P_3"]


def test_a_protein_with_one_eligible_pocket_contributes_one():
    """`per_protein` is a cap, not a quota: nothing is padded in that fails the filters."""
    rows = [_row(acc="P1", pocket="P_1"), _row(acc="P1", pocket="P_2", vol=9000.0)]
    anns = {(r["UniProt ID"], r["Pocket Filename"]): _ann() for r in rows}
    pan, _, _, _ = panel.select_panel(rows, anns, {}, {}, panel.Options(per_protein=2))
    assert len(pan["P1"]) == 1


def test_funnel_still_sums_when_an_annotation_is_missing():
    """A missing annotation is a rejection; if it is not counted the funnel reports more
    eligible pockets than passed, and the manifest contradicts its own table."""
    rows = [_row(acc="A"), _row(acc="B"), _row(acc="C")]
    anns = {("C", rows[2]["Pocket Filename"]): _ann()}
    pan, elig, funnel, rej = panel.select_panel(rows, anns, {}, {}, OPTS)
    assert funnel[-1][1] == len(elig) == 1
    assert len(rej) == 2


def test_a_missing_box_dimension_fails_the_grid_guard():
    """Defaulting a missing size to 0 made a malformed row a 1-point box, which is
    'legal' and silently becomes a dockable target."""
    ann = _ann(); ann.pop("box_size_y")
    assert panel.classify(_row(), ann, [], None, OPTS) == "grid"


def test_defaults_come_from_the_module_not_from_test_literals():
    """A test that hardcodes thresholds stays green when the production default moves —
    the failure mode that has already bitten this repo twice."""
    o = panel.Options()
    assert (o.volume_min, o.volume_max, o.plddt_min, o.per_protein) == (
        panel.VOLUME_MIN, panel.VOLUME_MAX, panel.PLDDT_MIN, panel.PER_PROTEIN)


# --- the workable grid ceiling ----------------------------------------------------------

def test_the_grid_guard_stops_below_what_unidock_formally_accepts():
    """Uni-Dock accepts a box just under its own 531,441-point limit and then crawls.

    On net3 a unit at ~530k points took 2 h 17 min against the ~10 min a normal one takes,
    so one such unit consumed an entire Colab session and the rest of the bundle was never
    drawn -- `incomplete_targets` could not reach zero and "re-run until nothing is
    outstanding" never terminated. Measured over 212 net3 targets: of the 197 that
    completed, the largest grid was 508,518 points; of the 15 that never finished, 7 were
    above it.
    """
    from polyphar import panel

    assert panel.WORKABLE_POINTS < panel.MAX_POINTS, "the point is that it stops EARLIER"
    assert panel.WORKABLE_POINTS == 508518

    # a box in the band between the two: legal to the engine, refused by us
    spacing = panel.SPACING
    side = 508518 ** (1 / 3)
    big = tuple(((side + 2) - 1) * spacing for _ in range(3))     # ~520k points
    n = [int(b / spacing) + 1 for b in big]
    assert panel.MAX_POINTS >= n[0] * n[1] * n[2] > panel.WORKABLE_POINTS, n
    assert not panel.grid_legal(big), "a box in the unusable band must be refused"
    # ...and the engine's own question can still be asked explicitly
    assert panel.grid_legal(big, limit=panel.MAX_POINTS)

    # an ordinary drug-site box is unaffected
    assert panel.grid_legal((20.0, 20.0, 20.0))


def test_the_grid_rejection_says_what_it_means():
    """It is no longer "exceeds Uni-Dock's grid limits": the box is inside those limits
    and still undockable, and a reader chasing the reason must not be sent to the wrong
    constant."""
    from polyphar import panel

    assert "Uni-Dock" not in panel.STEP_REASON["grid"]
    assert "workable" in panel.STEP_REASON["grid"]
