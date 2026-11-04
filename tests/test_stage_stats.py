"""Tests for the stats stage: the guarantees the design rests on.

Each of these is a property the spec (docs/superpowers/specs/2026-11-03-statistics-design.md)
says must hold, not an incidental implementation detail:
percentiles uniform per target BY CONSTRUCTION; `M` weighting proteins and not targets;
the permutation test finding a planted effect and finding nothing in null data; the strata
edges stable, recorded and compared against the pre-registered ones; Holm applied to
exactly the two primary tests; the descriptive surface carrying no p-value; `--nprocs`
changing nothing; a finished stage reporting 0 units.
"""
import csv
import json
import random
import types
from pathlib import Path

import pytest

from polyphar.common import RunLayout
from polyphar.stages import stats

FIX = Path(__file__).parent / "fixtures"


# --- step 1: percentiles ---------------------------------------------------------------

def test_percentiles_are_uniform_on_every_target_by_construction():
    """A generous pocket scores everything well; its percentiles still have mean 0.5 and
    the same spread as every other target's. That is how pocket generosity is REMOVED
    rather than adjusted for."""
    generous = [-12.0, -11.5, -11.0, -10.5, -10.0]      # everything binds well
    mean_pocket = [-7.0, -6.0, -5.0, -4.0, -3.0]
    for scores in (generous, mean_pocket):
        p = stats.percentiles(scores)
        assert sum(p) / len(p) == pytest.approx(0.5)
        assert sorted(p) == pytest.approx([(i + 0.5) / 5 for i in range(5)])
    assert stats.percentiles(generous) == pytest.approx(stats.percentiles(mean_pocket))


def test_percentile_direction_is_more_negative_is_better():
    p = stats.percentiles([-12.0, -3.0])
    assert p[0] > p[1] and p[0] == pytest.approx(0.75)


def test_percentile_ties_take_the_midrank():
    p = stats.percentiles([-5.0, -5.0, -1.0])
    assert p[0] == p[1] == pytest.approx(1.0 - (1.5 - 0.5) / 3)
    assert sum(p) / 3 == pytest.approx(0.5)


def test_percentiles_of_a_random_target_are_uniform():
    rng = random.Random(7)
    p = stats.percentiles([rng.gauss(-8, 2) for _ in range(200)])
    assert sum(p) / len(p) == pytest.approx(0.5)
    assert min(p) == pytest.approx(0.5 / 200) and max(p) == pytest.approx(1 - 0.5 / 200)


# --- step 2: M -------------------------------------------------------------------------

def test_m_weights_proteins_not_targets():
    """A protein contributing two pockets must not count twice. P1 has two pockets, P2
    one: M is the mean of (best of P1's two) and (P2's one), i.e. (0.9 + 0.1)/2."""
    per_target = {"a_P_1": {"L": 0.9}, "a_P_2": {"L": 0.5}, "b_P_1": {"L": 0.1}}
    protein_of = {"a_P_1": "P1", "a_P_2": "P1", "b_P_1": "P2"}
    m, n_proteins, n_targets = stats.m_statistic(per_target, protein_of)["L"]
    assert m == pytest.approx(0.5) and n_proteins == 2 and n_targets == 3
    # the mean over TARGETS would be (0.9 + 0.5 + 0.1)/3 = 0.5 only by coincidence here,
    # so check a case where they differ
    per_target["a_P_2"] = {"L": 0.9}
    m2, _, _ = stats.m_statistic(per_target, protein_of)["L"]
    assert m2 == pytest.approx(0.5)            # unchanged: P1 still contributes its best
    assert (0.9 + 0.9 + 0.1) / 3 == pytest.approx(0.6333333)   # the wrong statistic


def test_m_ignores_targets_not_on_the_panel():
    per_target = {"a_P_1": {"L": 1.0}, "ghost_P_1": {"L": 0.0}}
    out = stats.m_statistic(per_target, {"a_P_1": "P1"})
    assert out["L"][0] == pytest.approx(1.0) and out["L"][2] == 1


def test_percentile_matrix_ranks_within_each_target_independently():
    m = stats.percentile_matrix({"t1": {"a": -10.0, "b": -1.0},
                                 "t2": {"a": -1.0, "b": -10.0}})
    assert m["t1"]["a"] > m["t1"]["b"] and m["t2"]["b"] > m["t2"]["a"]


# --- step 3: strata ---------------------------------------------------------------------

def test_strata_edges_are_the_preregistered_ones_for_the_real_library_profile():
    """The edges are a pre-registered constant; the code must produce them from quintiles
    and tertiles, not from anything tuned later."""
    assert stats.PREREGISTERED_EDGES == {"heavy_atoms": [17.0, 22.0, 26.0, 32.0],
                                         "rotatable_bonds": [4.0, 7.0]}
    values = list(range(1, 101))
    assert stats.quantile_edges(values, (0.2, 0.4, 0.6, 0.8)) == [20.8, 40.6, 60.4, 80.2]


def test_bin_of_uses_half_open_intervals_from_the_edges():
    edges = [17.0, 22.0, 26.0, 32.0]
    assert [stats.bin_of(v, edges) for v in (1, 16, 17, 21, 22, 31, 32, 99)] == \
        [0, 0, 1, 1, 2, 3, 4, 4]


def test_build_strata_is_stable_and_records_everything(tmp_path):
    library = [{"ID": str(i), "name": f"L{i}", "heavy_atoms": i,
                "rotatable_bonds": i % 9} for i in range(1, 61)]
    a = stats.build_strata(library, None)
    b = stats.build_strata(list(library), None)
    assert a["edges"] == b["edges"] and a["of"] == b["of"]
    assert sum(a["sizes"].values()) == 60
    assert set(a) >= {"edges", "sizes", "quantiles", "bin_rule", "preregistered",
                      "matches_preregistration", "rotatable_bond_definition"}
    assert a["matches_preregistration"] is False       # not the real library


def test_strata_use_the_ligands_stage_rotatable_bond_definition(tmp_path):
    """RDKit's Lipinski count differs for 1,316 of 2,083 e-Drug3D records, so the strata
    must come from ligands.csv and never be recomputed."""
    path = tmp_path / "ligands.csv"
    path.write_text("record_index,ID,name,heavy_atoms,mw,rotatable_bonds,kept,exclusion\n"
                    "0,1,A,20,300,5,True,\n"
                    "1,2,B,40,700,20,False,rotatable_bonds=20>15\n")
    lib = stats.read_library(path)
    assert [r["ID"] for r in lib] == ["1"]             # excluded records are not in it
    assert lib[0]["rotatable_bonds"] == 5
    assert "ligands.csv" in stats.build_strata(lib, None)["rotatable_bond_definition"] \
        or "ligsize" in stats.build_strata(lib, None)["rotatable_bond_definition"]


# --- step 4: the test --------------------------------------------------------------------

def _synthetic(n=400, effect=0.0, seed=1):
    """Two strata; positives are the small-M-by-size half unless `effect` lifts them."""
    rng = random.Random(seed)
    m, pos, strat = [], [], []
    for i in range(n):
        s = "A" if i % 2 else "B"
        base = 0.3 if s == "A" else 0.7          # the size effect: strata differ a lot
        is_pos = i % 20 in (0, 1)        # positives in BOTH strata
        m.append(base + rng.gauss(0, 0.05) + (effect if is_pos else 0.0))
        pos.append(is_pos)
        strat.append(s)
    return m, pos, strat


def test_permutation_recovers_a_planted_effect():
    m, pos, strat = _synthetic(effect=0.25)
    res = stats.permutation_test(m, pos, strat, n_perm=2000, seed=3)
    assert res["statistic"] > 0.2
    assert res["p_value"] < 0.01
    assert res["n_positive"] == sum(pos) and res["n_strata"] == 2


def test_permutation_finds_nothing_in_null_data():
    """Null data over 30 independent draws: the statistic stays near zero and the p-values
    behave like a uniform, i.e. the test does not manufacture a result. One dataset would
    only show that this dataset's p was large."""
    ps = []
    for seed in range(30):
        m, pos, strat = _synthetic(effect=0.0, seed=seed)
        res = stats.permutation_test(m, pos, strat, n_perm=400, seed=100 + seed)
        assert abs(res["statistic"]) < 0.08
        ps.append(res["p_value"])
    assert sum(1 for p in ps if p < 0.05) <= 4         # expect ~1.5 of 30
    assert 0.3 < sum(ps) / len(ps) < 0.7


def test_permutation_is_seeded_and_reproducible():
    m, pos, strat = _synthetic(effect=0.1)
    a = stats.permutation_test(m, pos, strat, n_perm=500, seed=11)
    b = stats.permutation_test(m, pos, strat, n_perm=500, seed=11)
    c = stats.permutation_test(m, pos, strat, n_perm=500, seed=12)
    assert a["p_value"] == b["p_value"]
    assert a["statistic"] == pytest.approx(b["statistic"]) == pytest.approx(c["statistic"])


def test_permutation_p_is_never_zero():
    m, pos, strat = _synthetic(effect=5.0)
    res = stats.permutation_test(m, pos, strat, n_perm=200, seed=5)
    assert res["p_value"] == pytest.approx(1 / 201)


def test_the_test_is_one_sided():
    """A strongly NEGATIVE effect must give a large p, not a small one."""
    m, pos, strat = _synthetic(effect=-0.3)
    res = stats.permutation_test(m, pos, strat, n_perm=1000, seed=4)
    assert res["statistic"] < 0 and res["p_value"] > 0.9


def test_size_matching_neutralises_a_pure_size_effect():
    """Positives that are all large molecules, where every large molecule ranks high:
    unstratified this is a huge difference, stratified it is nothing."""
    m, pos, strat = [], [], []
    rng = random.Random(2)
    for i in range(300):
        big = i % 3 == 0
        m.append((0.8 if big else 0.3) + rng.gauss(0, 0.02))
        pos.append(big and i % 9 == 0)          # positives only among the big molecules
        strat.append("big" if big else "small")
    stratified = stats.permutation_test(m, pos, strat, n_perm=2000, seed=6)
    pooled = stats.permutation_test(m, pos, ["one"] * len(m), n_perm=2000, seed=6)
    assert pooled["statistic"] > 0.3 and pooled["p_value"] < 0.01
    # the size effect is gone: what is left is noise, orders of magnitude smaller
    assert abs(stratified["statistic"]) < pooled["statistic"] / 20


def test_stratified_difference_weights_by_the_positives_composition():
    m = [1.0, 0.0, 0.0, 10.0, 0.0, 0.0]
    pos = [True, False, False, True, False, False]
    strat = ["A", "A", "A", "B", "B", "B"]
    stat, detail = stats.stratified_difference(m, pos, strat)
    assert detail["A"]["weight"] == pytest.approx(0.5)
    assert stat == pytest.approx(0.5 * 1.0 + 0.5 * 10.0)


def test_strata_without_a_control_or_a_positive_are_dropped():
    m = [1.0, 0.5, 0.9]
    pos = [True, False, True]
    strat = ["A", "A", "LONELY"]
    stat, detail = stats.stratified_difference(m, pos, strat)
    assert list(detail) == ["A"] and stat == pytest.approx(0.5)


# --- multiple testing ---------------------------------------------------------------------

def test_holm_matches_the_textbook():
    assert stats.holm([0.01, 0.04]) == pytest.approx([0.02, 0.04])
    assert stats.holm([0.5, 0.001]) == pytest.approx([0.5, 0.002])
    assert stats.holm([0.9, 0.9]) == pytest.approx([1.0, 1.0])


def test_holm_is_monotone():
    out = stats.holm([0.03, 0.04, 0.9])
    assert out == sorted(out)


def test_benjamini_hochberg_matches_the_textbook():
    q = stats.benjamini_hochberg([0.01, 0.02, 0.03, 0.9])
    assert q == pytest.approx([0.04, 0.04, 0.04, 0.9])


# --- step 7: the descriptive surface --------------------------------------------------------

def test_the_surface_refuses_a_p_value_column(tmp_path):
    rows = [{"panel": "x", "p_value": 0.01}]
    with pytest.raises(ValueError, match="no p-values"):
        stats.write_no_pvalues(tmp_path / "s.csv", ["panel", "p_value"], rows)
    assert not (tmp_path / "s.csv").exists()


@pytest.mark.parametrize("column", [
    "p", "pval", "p_adj", "q_value", "fdr", "holm", "significance", "signif", "stars",
    # every one of these walked straight through the `[_\s]`-bounded denylist that this
    # allowlist replaced
    "p-value", "p.value", "adj.P.Val", "Pr(>|t|)", "pFDR", "bonferroni", "Benjamini",
    "log10p", "neglog10p", "prob", "probability", "confidence", "CI95", "ci_low",
    "P(perm)", "p<0.05", "p:value", "p/value", "p|value", "pvals", "p1", "evidence",
    "verdict", "two_sided_p", "exact_p", "Prob(>F)", "PValue", "q.val", "FDR.adj",
    "р_value",            # Cyrillic er
    "p​value",            # zero-width space
    "\U0001d45d",              # mathematical italic small p
    # and anything else at all: the allowlist refuses unexpected names, p-value or not
    "scratch", "note", "comment",
])
def test_the_surface_refuses_every_spelling_of_a_p_value(tmp_path, column):
    with pytest.raises(ValueError):
        stats.write_no_pvalues(tmp_path / "s.csv", ["panel", column], [])
    assert not (tmp_path / "s.csv").exists()
    # and smuggled through a row rather than the column list
    with pytest.raises(ValueError):
        stats.write_no_pvalues(tmp_path / "s.csv", list(stats.SURFACE_COLUMNS),
                               [{"panel": "x", column: 0.01}])
    assert not (tmp_path / "s.csv").exists()


def test_the_surface_writer_is_an_allowlist_not_a_denylist(tmp_path):
    """The point of the change: a name does not have to LOOK like a p-value to be refused.
    Only the expected fields may be written."""
    stats.write_no_pvalues(tmp_path / "ok.csv", list(stats.SURFACE_COLUMNS), [])
    assert (tmp_path / "ok.csv").is_file()
    with pytest.raises(ValueError, match="may carry ONLY"):
        stats.write_no_pvalues(tmp_path / "s.csv",
                               list(stats.SURFACE_COLUMNS) + ["harmless_extra"], [])


def test_the_surface_markdown_is_guarded_too(tmp_path):
    """surface.md is what gets pasted into a draft; it had no guard at all."""
    rows = stats.surface_rows("net1", "hcc", {"a": [-10.0]}, {"b": [-5.0]},
                              thresholds=(-9.0,), min_targets=(1,))
    assert "No p-values" in stats.surface_markdown(rows)
    for r in rows:
        r["p-value"] = 0.001
    with pytest.raises(ValueError, match="no p-values"):
        stats.surface_markdown(rows)


def test_the_surface_refuses_a_p_value_smuggled_in_a_row(tmp_path):
    """Columns are not the only way in: a row carrying an extra key is refused too."""
    with pytest.raises(ValueError):
        stats.write_no_pvalues(tmp_path / "s.csv", ["panel"],
                               [{"panel": "x", "p_value": 0.01}])


def test_surface_rows_count_what_the_field_reports():
    pos = {"a": [-10.0, -9.0, -5.0], "b": [-5.0, -5.0, -5.0]}
    ctl = {"c": [-10.0, -10.0, -10.0]}
    rows = stats.surface_rows("panel", "hcc", pos, ctl, thresholds=(-9.0,),
                              min_targets=(1, 2))
    by_k = {r["min_targets"]: r for r in rows}
    assert by_k[1]["n_positive_meeting"] == 1 and by_k[2]["n_positive_meeting"] == 1
    assert by_k[2]["frac_control_meeting"] == 1.0


# --- the positive-control gate --------------------------------------------------------

def _pc(**over):
    """A complete, passing positive_control.json. The gate now needs the three readouts,
    their n, and the panel-coverage pair: a partial file is `unknown`, which is the point
    of the minimum-n and coverage rules, so every test that wants to exercise ONE
    threshold starts from a passing whole and overrides it."""
    base = {"median_self_rank_percentile": 0.01,
            "redock_rmsd_success_fraction": 0.80,
            "library_top_decile_fraction": 0.80,
            "pose_in_pocket_fraction": 0.9,
            "n_redocked": 24, "n_redock_rmsd": 24, "n_library_percentile": 24,
            "panel_targets_scored": 482, "panel_targets_total": 482}
    base.update(over)
    return base


def test_a_missing_positive_control_is_not_a_pass(tmp_path):
    verdict, detail = stats.positive_control_verdict(tmp_path / "absent.json")
    assert verdict.startswith("unknown") and "NOT interpretable" in verdict and detail == {}


def test_a_failing_positive_control_makes_the_result_uninterpretable(tmp_path):
    p = tmp_path / "pc.json"
    p.write_text(json.dumps(_pc(median_self_rank_percentile=0.4,
                                redock_rmsd_success_fraction=0.2,
                                library_top_decile_fraction=0.1)))
    verdict, _ = stats.positive_control_verdict(p)
    assert verdict.startswith("no:") and "UNINTERPRETABLE" in verdict
    assert "self-rank" in verdict and "redocking success" in verdict
    assert "top decile" in verdict


def test_a_passing_positive_control_says_yes(tmp_path):
    p = tmp_path / "pc.json"
    p.write_text(json.dumps(_pc()))
    assert stats.positive_control_verdict(p)[0] == "yes"


def test_the_preregistered_thresholds_are_the_spec_ones():
    """Every constant the pre-registration pins, pinned here too. Each of these could be
    changed, departing from the pre-registration, with an otherwise green suite."""
    assert stats.SELF_RANK_MAX == 0.05
    assert stats.REDOCK_RMSD_MAX == 2.0 and stats.REDOCK_SUCCESS_MIN == 0.50
    assert stats.LIBRARY_TOP_DECILE_MIN == 0.50
    assert (stats.MIN_SELF_RANK_N, stats.MIN_REDOCK_N, stats.MIN_LIBRARY_N) == (10, 10, 10)
    assert stats.MIN_PANEL_SCORED_FRACTION == 0.98
    assert not hasattr(stats, "POSE_IN_POCKET_MIN"), (
        "pose-in-pocket is descriptive now: its null is above the old 0.70 gate")
    assert stats.N_PERMUTATIONS == 10000
    assert stats.SEED == 42
    assert stats.MIN_COVERAGE == 0.98
    assert stats.HA_QUANTILES == (0.2, 0.4, 0.6, 0.8)
    assert stats.RB_QUANTILES == (1 / 3, 2 / 3)
    assert stats.HA_QUANTILES_FINE == (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    assert stats.BALANCE_FLAG_SD == 0.6
    assert stats.PREREGISTERED_EDGES == {"heavy_atoms": [17.0, 22.0, 26.0, 32.0],
                                         "rotatable_bonds": [4.0, 7.0]}


# --- the positive-control gate validates its input (an invalid control is never a pass) --

@pytest.mark.parametrize("data", [
    {"median_self_rank_percentile": -1.0, "pose_in_pocket_fraction": 99.0},
    {"median_self_rank_percentile": "0.0", "redock_rmsd_success_fraction": "1e9"},
    {"median_self_rank_percentile": False, "redock_rmsd_success_fraction": True},
    {"median_self_rank_percentile": "-Infinity", "library_top_decile_fraction": "Infinity"},
    {"median_self_rank_percentile": "NaN", "redock_rmsd_success_fraction": "NaN"},
    {"median_self_rank_percentile": 1.5},
    {"library_top_decile_fraction": 1.5},
    {"redock_rmsd_success_fraction": -0.5},
    # a readout that is simply absent is not a pass either
    {"redock_rmsd_success_fraction": None},
    {"library_top_decile_fraction": None},
    # ... and neither is one over too few observations, or a half-scored panel
    {"n_redocked": 9}, {"n_redock_rmsd": 1}, {"n_library_percentile": 0},
    {"n_redocked": None}, {"panel_targets_scored": 382, "panel_targets_total": 482},
    {"panel_targets_total": None},
])
def test_an_invalid_positive_control_is_never_a_pass(tmp_path, data):
    """`yes` means the redocking was run and passed. A negative rank, a pose fraction of
    99, a boolean or a non-finite value is not a near-miss: it must land in `unknown`."""
    p = tmp_path / "pc.json"
    full = _pc(**data)
    for k, v in data.items():
        if v is None:
            full.pop(k)
    p.write_text(json.dumps(full))
    verdict, _ = stats.positive_control_verdict(p)
    assert verdict.startswith("unknown"), verdict
    assert not verdict.startswith("yes")


def test_the_positive_control_failure_message_does_not_round_the_number_away(tmp_path):
    """0.0500001 fails, and the message must not print it as "0.05 is not within 0.05"."""
    p = tmp_path / "pc.json"
    p.write_text(json.dumps(_pc(median_self_rank_percentile=0.0500001)))
    verdict, _ = stats.positive_control_verdict(p)
    assert verdict.startswith("no:") and "0.0500001" in verdict


# --- the positive set -----------------------------------------------------------------

def _library(n=60):
    return [{"ID": str(i), "name": f"DRUG{i}", "heavy_atoms": 10 + (i % 30),
             "rotatable_bonds": i % 9} for i in range(1, n + 1)]


def test_a_missing_positive_set_fails_clearly(tmp_path):
    with pytest.raises(FileNotFoundError) as exc:
        stats.read_positive_sets(tmp_path / "nope.csv", _library())
    assert "guideline reading, not code" in str(exc.value)
    assert "positive_sets.csv" in str(exc.value)


def test_the_shipped_fixture_parses_and_names_what_it_drops():
    lib = [{"ID": "1", "name": "SORAFENIB", "heavy_atoms": 29, "rotatable_bonds": 9},
           {"ID": "2", "name": "LENVATINIB", "heavy_atoms": 30, "rotatable_bonds": 7},
           {"ID": "3", "name": "REGORAFENIB", "heavy_atoms": 30, "rotatable_bonds": 9},
           {"ID": "4", "name": "DOXORUBICIN", "heavy_atoms": 39, "rotatable_bonds": 5},
           {"ID": "5", "name": "DONEPEZIL", "heavy_atoms": 28, "rotatable_bonds": 6},
           {"ID": "6", "name": "MEMANTINE", "heavy_atoms": 12, "rotatable_bonds": 0}]
    sets, unresolved = stats.read_positive_sets(FIX / "positive_sets.csv", lib)
    assert set(sets) == {"hcc", "ad"}
    assert set(sets["hcc"]) == {"1", "2", "3", "4"}
    assert sets["hcc"]["2"]["name"] == "LENVATINIB"       # resolved by name, no ID given
    assert {v["subgroup"] for v in sets["hcc"].values()} == {"targeted", "cytotoxic"}
    assert len(unresolved) == 1 and "RAMUCIRUMAB" in unresolved[0][0]


def test_a_name_that_disagrees_with_the_id_is_reported_not_swallowed():
    lib = [{"ID": "1", "name": "SORAFENIB", "heavy_atoms": 29, "rotatable_bonds": 9}]
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as fh:
        fh.write("set,ID,name,subgroup,source\nhcc,1,LENVATINIB,targeted,x\n")
        path = fh.name
    sets, unresolved = stats.read_positive_sets(path, lib)
    assert sets["hcc"]["1"]["name"] == "SORAFENIB"
    assert unresolved and "not 'LENVATINIB'" in unresolved[0][1]


# --- the stage end to end -------------------------------------------------------------

def _write_csv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns)
        w.writeheader()
        for r in rows:
            w.writerow(r)


HCC_POSITIVES = ("3", "11", "19", "27", "35", "43")
AD_POSITIVES = ("8",)


def _run_folder(tmp_path, name="net1", n_ligands=120, n_proteins=8, effect=2.0,
                positives=HCC_POSITIVES, replicates=2, two_pockets_for=("P001",)):
    """A tiny but structurally real run folder: ligands.csv, final_targets.csv,
    df_results_<label>.csv and raw replicate SDFs."""
    rng = random.Random(42)
    run = tmp_path / name
    library = _library(n_ligands)
    _write_csv(run / "07_docking" / "library" / "ligands.csv",
               ["record_index", "ID", "name", "heavy_atoms", "rotatable_bonds", "kept"],
               [{"record_index": i, "kept": "True", **{k: r[k] for k in
                 ("ID", "name", "heavy_atoms", "rotatable_bonds")}}
                for i, r in enumerate(library)])

    targets = []
    for p in range(n_proteins):
        acc = f"P{p:03d}"
        pockets = ["P_1", "P_2"] if acc in two_pockets_for else ["P_1"]
        for pocket in pockets:
            targets.append({"Node": f"{acc}_x", "UniProt_ID": acc, "Gene_Name": "G",
                            "Organism": "Human", "Structure_ID": f"{acc}s",
                            "Similarity": "", "Protein": f"{acc}_1abc_A",
                            "Chain": "A", "Pocket": pocket,
                            "Pocket_Path": "x", "Receptor_Path": "y", "Is_Hub": False})
    _write_csv(run / "06_targets" / "final_targets.csv", list(targets[0]), targets)

    docking = run / "07_docking"
    for t in targets:
        label = f'{t["Protein"]}_{t["Pocket"]}'
        rows, per_rep = [], {i: {} for i in range(replicates)}
        for r in library:
            base = -5.0 - 0.1 * r["heavy_atoms"] + rng.gauss(0, 0.5)
            if r["ID"] in positives:
                base -= effect
            draws = [base + rng.gauss(0, 0.05) for _ in range(replicates)]
            mean = sum(draws) / len(draws)
            rows.append({"Ligand": r["name"], "ID": r["ID"], "Pocket": label,
                         "n": replicates, "score_mean": f"{mean:.4f}",
                         "score_sd": "0.1", "score_se": "0.05",
                         "score_min": f"{min(draws):.4f}", "score_max": f"{max(draws):.4f}"})
            for i, d in enumerate(draws):
                per_rep[i][r["ID"]] = d
        _write_csv(docking / f"df_results_{label}.csv",
                   ["Ligand", "ID", "Pocket", "n", "score_mean", "score_sd", "score_se",
                    "score_min", "score_max"], rows)
        d = docking / label / "replicates"
        d.mkdir(parents=True, exist_ok=True)
        for i in range(replicates):
            (d / f"r{i}.sdf").write_text("".join(_sdf_record(lid, s)
                                                 for lid, s in per_rep[i].items()))
            (d / f"r{i}.done").write_text("done\n")
    return run, library


def _sdf_record(lid, score):
    return (f"EDRUG3D_{lid}\n  polyphar\n\n"
            "  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
            f"> <ID>\n{lid}\n\n> <minimizedAffinity>\n{score:.4f}\n\n$$$$\n")


def _args(run, **kw):
    base = dict(run_dir=str(run), nprocs=1, dry_run=False, no_notify=True,
                positive_set="inputs/positive_sets.csv",
                positive_control="inputs/positive_control.json",
                control_panel=None, permutations=500, permutation_seed=7,
                min_coverage=stats.MIN_COVERAGE)
    base.update(kw)
    return types.SimpleNamespace(**base)


def _positive_set(run, hcc=HCC_POSITIVES, ad=AD_POSITIVES):
    rows = [{"set": "hcc", "ID": i, "name": f"DRUG{i}",
             "subgroup": "targeted" if n % 2 == 0 else "cytotoxic", "source": "NCCN"}
            for n, i in enumerate(hcc)]
    rows += [{"set": "ad", "ID": i, "name": f"DRUG{i}", "subgroup": "", "source": "AAN"}
             for i in ad]
    _write_csv(run / "inputs" / "positive_sets.csv",
               ["set", "ID", "name", "subgroup", "source"], rows)


class _Progress:
    def __init__(self):
        self.n = 0

    def advance(self, n, size=0.0):
        self.n += n


def _rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def test_stage_end_to_end_writes_every_output(tmp_path):
    run, _lib = _run_folder(tmp_path)
    _positive_set(run)
    layout = RunLayout(run)
    counts = stats.run(layout, _args(run), _Progress())

    d = run / "08_stats"
    for name in (stats.STRATA, stats.RESULTS, stats.PER_DRUG, stats.SURFACE,
                 stats.SURFACE_MD, stats.SUMMARY_MD, "manifest.json", "failed.tsv"):
        assert (d / name).is_file(), name
    assert (d / "m_statistic_net1.csv").is_file()
    assert counts["proteins"] == 8 and counts["targets"] == 9     # one protein, 2 pockets
    assert counts["positives"] == len(HCC_POSITIVES) + len(AD_POSITIVES)
    assert counts["failed"] == 0
    assert counts["replicates"] == 2

    results = _rows(d / stats.RESULTS)
    primary = [r for r in results if r["family"] == "primary"]
    assert len(primary) == 1 and primary[0]["drug_set"] == "hcc"
    assert float(primary[0]["statistic"]) > 0
    assert float(primary[0]["p_value"]) < 0.05
    # the planted effect survives per replicate (step 5), and no analytic SE is reported
    assert float(primary[0]["stat_r0"]) > 0 and float(primary[0]["stat_r1"]) > 0
    assert primary[0]["replicate_sign_jitter_check"] == "True"
    assert "score_se" not in results[0] and "se" not in results[0]


def test_the_control_set_is_described_in_the_output(tmp_path):
    """"Other approved drugs of the same size" is load-bearing and easy to miss."""
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    d = run / "08_stats"
    row = _rows(d / stats.RESULTS)[0]
    assert "approved" in row["control_set"] and "e-Drug3D" in row["control_set"]
    assert "approved" in (d / stats.SUMMARY_MD).read_text()
    assert "approved" in (d / stats.SURFACE_MD).read_text()


def test_the_surface_file_has_no_p_value_column(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    text = (run / "08_stats" / stats.SURFACE).read_text()
    header = text.splitlines()[0].split(",")
    assert not any(stats.PVALUE_LIKE.search(c) for c in header)
    md = (run / "08_stats" / stats.SURFACE_MD).read_text()
    assert "no p-values" in md.lower() or "No p-values" in md


def test_strata_json_is_written_before_any_outcome_is_read(tmp_path, monkeypatch):
    """The edges must exist before the scores do, or they could be tuned to the outcome."""
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    seen = {}
    real = stats._load

    def spy(jobs, meta, args, progress, failures, dropped=None):
        seen["strata_existed"] = (run / "08_stats" / stats.STRATA).is_file()
        return real(jobs, meta, args, progress, failures, dropped)

    monkeypatch.setattr(stats, "_load", spy)
    stats.run(RunLayout(run), _args(run), _Progress())
    assert seen["strata_existed"] is True
    data = json.loads((run / "08_stats" / stats.STRATA).read_text())
    assert data["edges"]["heavy_atoms"] and data["sizes"]
    assert data["preregistered"] == stats.PREREGISTERED_EDGES


def test_holm_is_applied_to_exactly_the_two_primary_tests(tmp_path):
    run, _ = _run_folder(tmp_path, name="net1")
    other, _ = _run_folder(tmp_path, name="cluster3", effect=0.0, n_proteins=6,
                           two_pockets_for=())
    _positive_set(run)
    stats.run(RunLayout(run), _args(run, control_panel=str(other)), _Progress())

    rows = _rows(run / "08_stats" / stats.RESULTS)
    primary = [r for r in rows if r["family"] == "primary"]
    assert len(primary) == 2
    assert {r["panel"] for r in primary} == {"net1", "cluster3"}
    for r in primary:
        assert r["adjustment"] == "holm(primary, m=2)"
        assert float(r["p_adjusted"]) >= float(r["p_value"])
    assert stats.holm([float(r["p_value"]) for r in primary]) == pytest.approx(
        [float(r["p_adjusted"]) for r in primary])
    # exploratory arms are NOT corrected, and say so, with their n on the row
    for r in rows:
        if r["family"] == "exploratory":
            assert r["p_adjusted"] == "" and "exploratory" in r["adjustment"]
            assert int(r["n_positive"]) >= 0
    subgroups = [r for r in rows if r["family"] == "subgroup"]
    assert subgroups and all(s["adjustment"].startswith("holm(subgroup") for s in subgroups)
    assert all(s["panel"] == "net1" for s in subgroups)


def test_the_negative_control_arm_is_absent_and_said_to_be_absent(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    rows = _rows(run / "08_stats" / stats.RESULTS)
    primary = [r for r in rows if r["family"] == "primary"]
    assert len(primary) == 1 and "--control-panel" in primary[0]["note"]
    assert "not supplied" in (run / "08_stats" / stats.SUMMARY_MD).read_text()


def test_the_negative_control_prediction_direction_holds_on_synthetic_panels(tmp_path):
    """HCC drugs are planted strong on net1 and neutral on the control panel, so
    HCC-on-HCC > HCC-on-AD must come out."""
    run, _ = _run_folder(tmp_path, name="net1", effect=0.8)
    other, _ = _run_folder(tmp_path, name="ad", effect=0.0, two_pockets_for=())
    _positive_set(run)
    stats.run(RunLayout(run), _args(run, control_panel=str(other)), _Progress())
    rows = {r["panel"]: r for r in _rows(run / "08_stats" / stats.RESULTS)
            if r["family"] == "primary"}
    assert float(rows["net1"]["statistic"]) > float(rows["ad"]["statistic"])
    assert float(rows["net1"]["p_value"]) < float(rows["ad"]["p_value"])
    assert "overlap" in (run / "08_stats" / stats.SUMMARY_MD).read_text()


def test_per_drug_results_are_labelled_exploratory_with_bh(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    rows = _rows(run / "08_stats" / stats.PER_DRUG)
    assert rows and all(r["status"] == "exploratory (BH FDR)" for r in rows)
    assert all(r["q_value_bh"] != "" for r in rows)
    assert all(float(r["q_value_bh"]) >= float(r["p_value"]) for r in rows)


def test_the_positive_control_verdict_reaches_every_primary_row(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    (run / "inputs").mkdir(parents=True, exist_ok=True)
    (run / "inputs" / "positive_control.json").write_text(json.dumps(
        _pc(median_self_rank_percentile=0.6, redock_rmsd_success_fraction=0.1,
            library_top_decile_fraction=0.1)))
    stats.run(RunLayout(run), _args(run), _Progress())
    rows = _rows(run / "08_stats" / stats.RESULTS)
    for r in rows:
        if r["family"] == "primary":
            assert r["interpretable"].startswith("no:")
            assert "UNINTERPRETABLE" in r["interpretable"]
    assert "UNINTERPRETABLE" in (run / "08_stats" / stats.SUMMARY_MD).read_text()


def test_unknown_positive_control_is_recorded_as_unknown(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    counts = stats.run(RunLayout(run), _args(run), _Progress())
    assert counts["positive_control"] == "unknown"


def test_nprocs_does_not_change_results(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    layout = RunLayout(run)
    stats.run(layout, _args(run, nprocs=1), _Progress())
    first = {p.name: p.read_bytes() for p in (run / "08_stats").glob("*")
             if p.name != "manifest.json"}
    for p in (run / "08_stats").glob("*"):
        p.unlink()
    stats.run(layout, _args(run, nprocs=4), _Progress())
    second = {p.name: p.read_bytes() for p in (run / "08_stats").glob("*")
              if p.name != "manifest.json"}
    assert first == second


def test_resume_reports_zero_units_and_touches_nothing(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    layout = RunLayout(run)
    args = _args(run)
    assert stats.units(layout, args).count > 0
    counts = stats.run(layout, args, _Progress())
    assert stats.units(layout, args).count == 0
    before = (run / "08_stats" / stats.RESULTS).stat().st_mtime_ns
    progress = _Progress()
    again = stats.run(layout, args, progress)
    assert again == counts and progress.n == 0
    assert (run / "08_stats" / stats.RESULTS).stat().st_mtime_ns == before


def test_a_changed_score_file_invalidates_the_resume(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    layout, args = RunLayout(run), _args(run)
    stats.run(layout, args, _Progress())
    assert stats.units(layout, args).count == 0
    path = sorted((run / "07_docking").glob("df_results_*.csv"))[0]
    path.write_text(path.read_text().replace("-5.", "-6."))
    assert stats.units(layout, args).count > 0


def test_units_fail_clearly_without_a_positive_set(tmp_path):
    run, _ = _run_folder(tmp_path)
    layout = RunLayout(run)
    with pytest.raises(FileNotFoundError, match="guideline reading"):
        stats.units(layout, _args(run))


def test_replicate_sdf_reader_takes_the_best_pose_per_ligand(tmp_path):
    path = tmp_path / "r0.sdf"
    path.write_text(_sdf_record("1", -5.0) + _sdf_record("1", -9.0)
                    + _sdf_record("2", -3.0))
    assert stats.read_replicate_best(path) == {"1": -9.0, "2": -3.0}


def test_replicate_sdf_reader_refuses_two_score_fields(tmp_path):
    path = tmp_path / "r0.sdf"
    path.write_text("EDRUG3D_1\n\n\n  0  0\nM  END\n> <ID>\n1\n\n"
                    "> <minimizedAffinity>\n-5.0\n\n> <docking_score>\n-6.0\n\n$$$$\n")
    with pytest.raises(ValueError, match="more than one score field"):
        stats.read_replicate_best(path)


def test_an_unfinished_replicate_is_not_read(tmp_path):
    run, _ = _run_folder(tmp_path)
    label = sorted(p.name for p in (run / "07_docking").glob("*") if p.is_dir()
                   and p.name != "library")[0]
    (run / "07_docking" / label / "replicates" / "r0.done").unlink()
    found = stats.replicate_paths(run / "07_docking", label)
    assert [i for i, _ in found] == [1]


# --- wiring ---------------------------------------------------------------------------

def test_stats_is_wired_as_a_stage():
    from polyphar import cli
    from polyphar.benchmark import ORDER, STAGE_DIRS
    assert cli.STAGES["stats"][0] == "polyphar.stages.stats"
    assert ORDER[-1] == "stats"
    assert ORDER.index("dock") < ORDER.index("stats")
    assert STAGE_DIRS["stats"] == "08_stats"


def test_stats_has_a_projected_estimate():
    from polyphar import benchmark
    state = dict(benchmark.PRIORS)
    assert benchmark.projected_units("stats", state).count > 0


def test_the_manifest_is_written_last_and_records_the_preregistration(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    manifest = json.loads((run / "08_stats" / "manifest.json").read_text())
    assert manifest["stage"] == "stats"
    assert manifest["versions"]["preregistration"].endswith("preregistration.md")
    assert manifest["versions"]["preregistered_heavy_atoms"] == [17.0, 22.0, 26.0, 32.0]
    assert manifest["counts"]["permutations"] == 500
    assert "signature" in manifest


def test_the_run_preflight_stops_without_a_positive_set(tmp_path):
    from polyphar.stages import run as run_stage
    run, _ = _run_folder(tmp_path)
    problems = run_stage.preflight_problems(RunLayout(run), _args(run), ["stats"])
    assert problems and "guideline reading" in problems[0]
    _positive_set(run)
    assert run_stage.preflight_problems(RunLayout(run), _args(run), ["stats"]) == []


def test_the_preregistration_exists_and_states_what_it_must():
    """The pre-registration is the point of the design; it must exist, must carry the
    actual bin edges, and must state the positive-control thresholds AND the consequence
    of failing them. This is the part most likely to be quietly dropped.

    A 10-line stub containing only the substrings the earlier version of this test looked
    for passed it, so the assertions below also cover the tie rule, the coverage rule, the
    seed, the statistic's definition, the word `tertile`, `Benjamini`, the 1,960-record
    basis and the `## Amendments` section itself. Even so, substring checks are a floor,
    not a guarantee: the sha256 pin in the next test is what actually fixes the text."""
    doc = stats.preregistration_path()
    text = doc.read_text()
    flat = " ".join(text.split())          # the document hard-wraps; phrases span lines
    assert "17, 22, 26, 32" in text and "4, 7" in text          # the computed edges
    assert "10,000" in text and "one-sided" in text
    assert "Holm" in text and "two primary tests" in text.lower()
    assert "HCC-on-HCC > HCC-on-AD" in text                     # direction, stated
    assert "top 5 %" in text or "top 5%" in text
    assert ">= 50 % within 2.0" in text          # the redocking-RMSD gate
    assert "top decile" in text                  # the library-percentile gate
    assert "minimum n" in text.lower()
    assert "UNINTERPRETABLE" in text
    assert "absent positive control is not a pass" in text.lower()
    # the code's thresholds are the document's
    assert f"{stats.SELF_RANK_MAX}" in text
    assert f"{stats.MIN_PANEL_SCORED_FRACTION:g}" in text
    assert f"{stats.REDOCK_RMSD_MAX:g}" in text
    # the statistic, exactly, and the pieces that make it what it is
    assert "mid-rank" in text.lower()                            # the tie rule
    assert "max over that protein's targets" in text
    assert "mean over PROTEINS" in text
    assert "1 - (rank - 0.5) / N" in text
    assert "tertile" in text.lower() and "quintile" in text.lower()
    assert "**Seed:** 42" in text or "Seed: 42" in text
    assert f"{stats.MIN_COVERAGE:g}" in text and "coverage" in text.lower()
    assert "Benjamini" in text
    assert "1,960" in text
    assert "## Amendments" in text
    # the two design flaws the reviews found must be stated, not only fixed in code
    assert "difference-of-significance" in flat
    assert "paired within-stratum permutation" in flat
    assert "within-stratum" in flat and "0.6 sd" in flat
    assert "Holm makes arm 2" in flat or "raw p is named as the honest figure" in flat
    assert "sensitivity" in text.lower()
    assert "allowlist" in text.lower()


def test_the_preregistration_is_the_one_this_code_was_written_against():
    """The document's sha256 is a module constant and is recorded in every manifest.

    Without this, nothing read the pre-registration: the manifest held its PATH as a
    string literal and the edges as code constants, so a document/code divergence passed
    silently and an amendment made after the numbers were known left no trace. Editing the
    document now fails this test until `stats.PREREGISTRATION_SHA256` and an entry under
    `## Amendments` are updated together."""
    doc = stats.preregistration_path()
    assert doc.is_file(), doc
    assert stats.PREREGISTRATION_SHA256, "PREREGISTRATION_SHA256 is not pinned"
    assert stats.preregistration_sha256() == stats.PREREGISTRATION_SHA256, (
        f"{doc} has changed. Update stats.PREREGISTRATION_SHA256 to "
        f"{stats.preregistration_sha256()} AND add a dated entry under '## Amendments' "
        f"saying what changed and why.")


def test_the_manifest_records_the_preregistration_hash(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    versions = json.loads((run / "08_stats" / "manifest.json").read_text())["versions"]
    assert versions["preregistration_sha256"] == stats.PREREGISTRATION_SHA256
    assert versions["preregistration_sha256_expected"] == stats.PREREGISTRATION_SHA256


def test_the_code_points_at_the_preregistration():
    source = Path(stats.__file__).read_text()
    assert "2026-11-04-preregistration.md" in source
    assert stats.PREREGISTERED_EDGES["heavy_atoms"] == [17.0, 22.0, 26.0, 32.0]


def test_nothing_here_ranks_on_pocketmatch_pmin():
    """Pmin == 1 is dominated by containment and ~49,600 artificial 1.0 rows sit at the
    top of any distribution built on it; a ranking added here must use Pmax."""
    source = Path(stats.__file__).read_text()
    assert "Pmin" in source and "Pmax, not Pmin" in source
    assert "Pocketsim" not in source


def test_no_mean_absolute_delta_standard_error_anywhere():
    """E|X - Y| = 1.128 sigma for a normal: mean |delta| / sqrt(N) is not a standard
    error and must not come back."""
    source = Path(stats.__file__).read_text()
    assert "1.128" in source
    assert "abs(" not in source.split("Step 5")[1].split("Step 6")[0]


# --- P1: the pre-registered prediction is TESTED, not read off two verdicts ------------

def test_the_prediction_has_its_own_contrast_arm(tmp_path):
    """HCC-on-HCC > HCC-on-AD must be a statistic with a p, not an inference from "one
    arm is significant and the other is not" (the difference-of-significance error)."""
    run, _ = _run_folder(tmp_path, name="net1", effect=0.8)
    other, _ = _run_folder(tmp_path, name="ad", effect=0.0, two_pockets_for=())
    _positive_set(run)
    counts = stats.run(RunLayout(run), _args(run, control_panel=str(other)), _Progress())
    rows = _rows(run / "08_stats" / stats.RESULTS)
    contrast = [r for r in rows if r["family"] == "contrast"]
    assert len(contrast) == 1 and counts["contrast_tests"] == 1
    c = contrast[0]
    assert "minus" in c["arm"]
    assert float(c["statistic"]) > 0 and float(c["p_value"]) < 0.05
    assert "paired" in c["direction"]
    # it is NOT in the Holm family of two
    assert c["p_adjusted"] == "" and "not in the Holm family" in c["adjustment"]
    assert len([r for r in rows if r["family"] == "primary"]) == 2
    text = (run / "08_stats" / stats.SUMMARY_MD).read_text()
    assert "difference-of-significance" in text
    assert "paired within-stratum permutation" in text
    assert "larger arm statistic" in text


def test_the_contrast_is_the_paired_difference_of_the_two_panels():
    a = {"1": 0.9, "2": 0.4, "3": 0.5}
    b = {"1": 0.2, "2": 0.4, "9": 0.1}
    assert stats.paired_difference(a, b) == pytest.approx({"1": 0.7, "2": 0.0})


def test_the_contrast_is_null_when_both_panels_behave_the_same(tmp_path):
    run, _ = _run_folder(tmp_path, name="net1", effect=0.8)
    other, _ = _run_folder(tmp_path, name="ad", effect=0.8, two_pockets_for=())
    _positive_set(run)
    stats.run(RunLayout(run), _args(run, control_panel=str(other)), _Progress())
    c = [r for r in _rows(run / "08_stats" / stats.RESULTS)
         if r["family"] == "contrast"][0]
    assert float(c["p_value"]) > 0.05


def test_the_negative_control_arm_names_the_raw_p_as_the_honest_one(tmp_path):
    """Holm makes the arm we WANT to fail harder to reject. Both p's are reported and the
    raw one is named."""
    run, _ = _run_folder(tmp_path, name="net1", effect=0.8)
    other, _ = _run_folder(tmp_path, name="ad", effect=0.0, two_pockets_for=())
    _positive_set(run)
    stats.run(RunLayout(run), _args(run, control_panel=str(other)), _Progress())
    rows = _rows(run / "08_stats" / stats.RESULTS)
    ctl = [r for r in rows if r["family"] == "primary" and r["panel"] == "ad"][0]
    assert ctl["p_value"] != "" and ctl["p_adjusted"] != ""
    assert "RAW p" in ctl["note"] and "backwards" in ctl["note"]
    text = (run / "08_stats" / stats.SUMMARY_MD).read_text()
    assert "honest figure" in text and "harder to\nreject" in text.replace("  ", " ") \
        or "harder to" in text


def test_a_control_panel_pointing_at_the_primary_run_is_refused(tmp_path):
    """It renamed the second panel <name>_control and Holm-corrected two identical
    p-values: a fabricated negative control that cannot fail."""
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    with pytest.raises(ValueError, match="this run folder"):
        stats.run(RunLayout(run), _args(run, control_panel=str(run)), _Progress())
    with pytest.raises(ValueError, match="cannot fail"):
        stats.units(RunLayout(run), _args(run, control_panel=str(run)))


# --- P2: the within-stratum size confound is diagnosed, not hidden --------------------

def test_balance_rows_fire_on_a_within_stratum_size_imbalance():
    """Stratification removes the BETWEEN-stratum size confound, not the within-stratum
    one; this is the diagnostic that says whether the objection applies here."""
    ids = [str(i) for i in range(40)]
    stratum = ["A"] * 40
    # positives are the 10 largest molecules of the only stratum
    ha = {str(i): 20 + i for i in range(40)}
    positive = [i >= 30 for i in range(40)]
    detail = {"A": {"weight": 1.0, "n_positive": 10, "n_control": 30}}
    rows, weighted, biggest = stats.balance_rows(detail, ids, positive, stratum, ha, {})
    assert len(rows) == 1
    assert rows[0]["ha_mean_positive"] > rows[0]["ha_mean_control"]
    assert weighted > stats.BALANCE_FLAG_SD and biggest == pytest.approx(weighted)
    assert rows[0]["flag"]
    assert "SIZE IMBALANCE" in stats.balance_flag(weighted, biggest)


def test_balance_stays_quiet_on_a_size_balanced_positive_set():
    ids = [str(i) for i in range(40)]
    stratum = ["A"] * 40
    ha = {str(i): 20 + i for i in range(40)}
    positive = [i % 4 == 0 for i in range(40)]      # evenly spread across the size range
    detail = {"A": {"weight": 1.0, "n_positive": 10, "n_control": 30}}
    _rows_, weighted, biggest = stats.balance_rows(detail, ids, positive, stratum, ha, {})
    assert abs(weighted) < stats.BALANCE_FLAG_SD
    assert stats.balance_flag(weighted, biggest) == ""


def test_the_stage_writes_a_balance_table_and_a_sensitivity_arm(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    counts = stats.run(RunLayout(run), _args(run), _Progress())
    balance = _rows(run / "08_stats" / stats.BALANCE)
    assert balance
    assert {"ha_mean_positive", "ha_mean_control", "ha_standardised_difference",
            "targets_mean_positive", "targets_mean_control"} <= set(balance[0])
    arms = {r["arm"] for r in balance}
    assert any("hcc on net1" == a for a in arms)
    sens = [r for r in _rows(run / "08_stats" / stats.RESULTS)
            if r["family"] == "sensitivity"]
    assert len(sens) == 1 and "fine HA bins" in sens[0]["arm"]
    assert sens[0]["p_adjusted"] == "" and "sensitivity" in sens[0]["adjustment"]
    assert "balance_flagged" in counts
    text = (run / "08_stats" / stats.SUMMARY_MD).read_text()
    assert "within-stratum" in text.lower() and "0.6 sd" in text
    assert "deciles" in text


def test_the_fine_strata_are_finer_and_recorded(tmp_path):
    library = [{"ID": str(i), "name": f"L{i}", "heavy_atoms": i,
                "rotatable_bonds": i % 9} for i in range(1, 101)]
    s = stats.build_strata(library, None)
    assert len(s["fine_sizes"]) > len(s["sizes"])
    assert len(s["fine_edges"]["heavy_atoms"]) == 9
    assert set(s["fine_of"]) == set(s["of"])


# --- P6: a file that reads but says nothing must raise, not report "no effect" ---------

def test_unparseable_scores_are_counted_and_not_absorbed(tmp_path):
    path = tmp_path / "df_results_x.csv"
    path.write_text("Ligand,ID,Pocket,n,score_mean\n"
                    "a,1,x,3,nan\nb,2,x,3,inf\nc,3,x,3,-7.0\nd,4,x,3,\ne,5,x,3,oops\n")
    scores, dropped = stats.read_df_results(path)
    assert scores == {"3": -7.0}
    assert len(dropped) == 4
    assert any("not finite" in why for _id, why in dropped)


def test_a_nan_score_would_corrupt_its_whole_target_so_it_is_rejected():
    """percentiles sorts with NaN present, and NaN compares false against everything, so
    the NaN is placed by position rather than by value: [-10, nan, -5, -8] comes out
    [.875, .625, .125, .375]. The NaN ligand ranks SECOND of four on this target, ahead of
    two real scores, and every ligand below it is displaced -- so the target's whole
    ranking is wrong, not just the NaN's own place in it. The mean is still exactly 0.5,
    so nothing downstream can notice."""
    p = stats.percentiles([-10.0, float("nan"), -5.0, -8.0])
    assert p == pytest.approx([0.875, 0.625, 0.125, 0.375])
    assert p[1] > p[2] and p[1] > p[3]     # the NaN outranks two real scores
    assert sum(p) / 4 == pytest.approx(0.5)      # and it looks perfectly healthy


def test_a_panel_whose_every_row_is_unparseable_raises(tmp_path):
    """It used to return panel.pooled = {label: {}} -- non-empty -- and write a manifest
    plus an arms table of blank statistics and blank p's that reads as "no effect"."""
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    for p in (run / "07_docking").glob("df_results_*.csv"):
        p.write_text("Ligand,ID,Pocket,n,score_mean\nA,1,x,2,NOT_A_NUMBER\n")
    with pytest.raises(RuntimeError, match="not one usable score"):
        stats.run(RunLayout(run), _args(run), _Progress())
    assert not (run / "08_stats" / "manifest.json").exists()
    failed = (run / "08_stats" / "failed.tsv").read_text()
    assert "unusable score_mean row" in failed


# --- P7: coverage ----------------------------------------------------------------------

def test_a_positive_below_coverage_is_named_not_silently_dropped(tmp_path):
    """Pre-registration section 9 promises every guideline drug dropped is named."""
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    for p in sorted((run / "07_docking").glob("df_results_*.csv")):
        rs = list(csv.DictReader(p.open()))
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rs[0]))
            w.writeheader()
            w.writerows([r for r in rs if r["ID"] != "3"])
    counts = stats.run(RunLayout(run), _args(run), _Progress())
    assert counts["positives_below_coverage"] == 1
    failed = (run / "08_stats" / "failed.tsv").read_text()
    assert "DRUG3" in failed and "dropped from the test" in failed
    summary = (run / "08_stats" / stats.SUMMARY_MD).read_text()
    assert "DRUG3" in summary and "did not clear the coverage rule" in summary


def test_coverage_is_measured_on_targets_as_well_as_proteins(tmp_path):
    """E[max of k uniforms] = k/(k+1): losing one pocket of a two-pocket protein moves M
    while protein coverage stays at 1.000, so a protein-only gate cannot see it."""
    run, _ = _run_folder(tmp_path, two_pockets_for=("P000", "P001", "P002", "P003"))
    _positive_set(run)
    # drop ligand "5" from every SECOND pocket: every protein is still covered
    for p in sorted((run / "07_docking").glob("df_results_*P_2.csv")):
        rs = list(csv.DictReader(p.open()))
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rs[0]))
            w.writeheader()
            w.writerows([r for r in rs if r["ID"] != "5"])
    stats.run(RunLayout(run), _args(run), _Progress())
    row = {r["ID"]: r for r in _rows(run / "08_stats" / "m_statistic_net1.csv")}["5"]
    assert float(row["coverage"]) == 1.0            # protein coverage sees nothing
    assert float(row["target_coverage"]) < 1.0      # target coverage does
    assert row["usable"] == "False"


# --- P9: the smaller ones --------------------------------------------------------------

def test_the_summary_arms_table_carries_the_interpretable_verdict(tmp_path):
    """The table is what gets pasted into a Methods draft, and it used to be
    byte-identical for pass, fail and absent positive control."""
    tables = {}
    for name, pc in (("absent", None),
                     ("pass", _pc()),
                     ("fail", _pc(median_self_rank_percentile=0.6,
                                  redock_rmsd_success_fraction=0.1,
                                  library_top_decile_fraction=0.1))):
        run, _ = _run_folder(tmp_path / name)
        _positive_set(run)
        if pc is not None:
            (run / "inputs").mkdir(parents=True, exist_ok=True)
            (run / "inputs" / "positive_control.json").write_text(json.dumps(pc))
        stats.run(RunLayout(run), _args(run), _Progress())
        text = (run / "08_stats" / stats.SUMMARY_MD).read_text()
        tables[name] = text[text.index("## Arms"):text.index("## The pre-registered")]
    assert tables["absent"] != tables["pass"] != tables["fail"]
    assert "UNINTERPRETABLE" in tables["absent"] and "UNINTERPRETABLE" in tables["fail"]
    assert "UNINTERPRETABLE" not in tables["pass"]
    assert "interpretable" in tables["pass"]


def test_bh_runs_over_the_resolvable_subset(tmp_path):
    """One drug with no control in its stratum used to wipe out every q-value."""
    ps = ["", 0.01, 0.02]
    rows = [{"p_value": p, "q_value_bh": "", "stratum": "A", "status": ""} for p in ps]
    usable = [r for r in rows if r["p_value"] != ""]
    q = stats.benjamini_hochberg([r["p_value"] for r in usable])
    assert len(q) == 2
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    drug = _rows(run / "08_stats" / stats.PER_DRUG)
    assert drug and all(r["q_value_bh"] != "" for r in drug if r["p_value"] != "")


def test_an_empty_positive_set_says_the_rows_are_missing_not_the_header(tmp_path):
    path = tmp_path / "positive_sets.csv"
    path.write_text("set,ID,name,subgroup,source\n")
    with pytest.raises(ValueError, match="header is there but the file has no rows"):
        stats.read_positive_sets(path, _library())
    path.write_text("wrong,columns\na,b\n")
    with pytest.raises(ValueError, match="needs a header"):
        stats.read_positive_sets(path, _library())


def test_one_positive_is_caveated_on_the_row(tmp_path):
    run, _ = _run_folder(tmp_path)
    _positive_set(run, hcc=("3",), ad=())
    stats.run(RunLayout(run), _args(run), _Progress())
    row = [r for r in _rows(run / "08_stats" / stats.RESULTS)
           if r["family"] == "primary"][0]
    assert row["n_positive"] == "1"
    assert "n_positive = 1" in row["note"] and "floor" in row["note"]


def test_the_replicate_sign_check_is_labelled_as_a_jitter_check(tmp_path):
    """corr(M_r0, M_r1) is about 0.9996, so this is True on essentially every arm. It must
    not sit beside a p-value reading as corroboration."""
    run, _ = _run_folder(tmp_path)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    header = _rows(run / "08_stats" / stats.RESULTS)[0]
    assert "replicate_sign_jitter_check" in header
    assert "replicate_sign_agrees" not in header
    text = (run / "08_stats" / stats.SUMMARY_MD).read_text()
    assert "not corroboration" in text.lower() and "0.9996" in text


def test_extra_replicates_are_reported_even_though_they_have_no_column(tmp_path):
    run, _ = _run_folder(tmp_path, replicates=5)
    _positive_set(run)
    stats.run(RunLayout(run), _args(run), _Progress())
    row = [r for r in _rows(run / "08_stats" / stats.RESULTS)
           if r["family"] == "primary"][0]
    assert "stat_r3+" in row["note"] and "r4=" in row["note"]


def test_run_stops_at_dock_by_default_but_stats_stays_a_stage():
    from polyphar import benchmark, cli
    assert benchmark.ORDER[-1] == "stats"
    args = cli.build_parser().parse_args(["run"])
    assert args.to_stage == "dock"
    assert cli.build_parser().parse_args(["run", "--to", "stats"]).to_stage == "stats"
    assert cli.build_parser().parse_args(["run", "--from", "stats"]).from_stage == "stats"


def test_the_runbook_documents_the_stats_stage():
    text = (Path(__file__).resolve().parents[1] / "docs" / "RUNBOOK.md").read_text()
    assert "stats" in text
    assert "inputs/positive_sets.csv" in text
    assert "inputs/positive_control.json" in text
