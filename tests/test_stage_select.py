import json
import os
import random
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from polyphar import cli
from polyphar.common import RunLayout
from polyphar.stages import select as sel

FIX = Path(__file__).parent / "fixtures"
SMALL = FIX / "pocketsim_small.csv"

HUB1 = "P00001_HUB1_HUMAN_1AAA_A_P_1"
HUB2 = "P00008_HUB2_HUMAN_1HHH_C_P_3"
P73 = "O15350_P73_HUMAN_3VD2_A_P_1"
N2 = "P00002_NEI2_HUMAN_1BBB_A_P_1"
N3 = "P00003_2AAA_HUMAN_3CCC_B_P_1"
N4 = "P00004_NEI4_HUMAN_AF-P00004-F1_P_12_1"
N9 = "P00009_NEI9_HUMAN_1III_A_P_1"


def small():
    return sel.read_pocketsim(SMALL)


def _graph(df=None, pmax=0.7):
    pairs = sel.unique_protein_pairs(small() if df is None else df, pmax)
    return pairs, sel.build_graph(pairs)


def _protein_info(extra_missing=()):
    """Protein_info.csv rows for every fixture protein (PDB rows and one AlphaFold row)."""
    rows = []
    for acc, entry, pdb, chain in [
            ("P00001", "HUB1_HUMAN", "1AAA", "A"), ("P00002", "NEI2_HUMAN", "1BBB", "A"),
            ("P00003", "2AAA_HUMAN", "3CCC", "B"), ("O15350", "P73_HUMAN", "3VD2", "A"),
            ("P00005", "NEI5_HUMAN", "1EEE", "A"), ("P00006", "NEI6_HUMAN", "1FFF", "A"),
            ("P00007", "NEI7_HUMAN", "1GGG", "A"), ("P00008", "HUB2_HUMAN", "1HHH", "C"),
            ("P00009", "NEI9_HUMAN", "1III", "A"), ("P00010", "LOW10_HUMAN", "1JJJ", "A")]:
        if acc in extra_missing:
            continue
        rows.append([acc, entry, pdb, chain, "False", "2.0", "X-ray", "0.9",
                     f"02_structures/{acc}/{pdb.lower()}_{chain}.pdb"])
    rows.append(["P00004", "NEI4_HUMAN", "-", "-", "True", "-", "-", "-",
                 "02_structures/P00004/AF-P00004-F1.pdb"])
    return pd.DataFrame(rows, columns=["UniProt ID", "Gene Name", "PDB", "Chain ID", "AlphaFold",
                                       "Resolution", "Structure method", "SIFTS coverage", "File"])


def _targets(tie_break, exclude=("O15350",), hub=None, info=True):
    _, g = _graph()
    ranking = sel.rank_hubs(g, tie_break)
    return sel.select_targets(g, ranking, set(exclude), tie_break, hub=hub,
                              protein_info=_protein_info() if info else None)


# --- pairs and graph -----------------------------------------------------------------

def test_reversed_pair_deduplicated():
    pairs, g = _graph()
    p1p2 = pairs[(pairs["UniProt ID"].isin(["P00001", "P00002"]))
                 & (pairs["UniProt ID_2"].isin(["P00001", "P00002"]))]
    assert len(p1p2) == 1                      # (A,B) and (B,A) give one row ...
    assert g.has_edge(HUB1, N2)                # ... and one edge
    assert g[HUB1][N2]["weight"] == pytest.approx(0.90)   # the max-Pmax pocket pair
    assert not g.has_node("P00002_NEI2_HUMAN_1BBB_A_P_2")
    keys = {frozenset(x) for x in zip(pairs["UniProt ID"], pairs["UniProt ID_2"])}
    assert len(keys) == len(pairs)             # one row per unordered protein pair


def test_pmax_filter_applies_after_the_dedup_and_keeps_threshold_rows():
    pairs, g = _graph()
    assert (pairs["Pmax"] >= 0.7).all()
    assert not any(n.startswith("P00010_") for n in g.nodes)    # 0.65 pair dropped
    assert sel.unique_protein_pairs(small(), 0.95)["Pmax"].tolist() == [0.95]


def test_pair_tie_goes_to_smallest_node_pair():
    # N4 - N9 (P_1) and N4 - N9b (P_3) tie at 0.72
    pairs, g = _graph()
    assert g.has_edge(N4, N9)
    assert not g.has_node("P00009_NEI9_HUMAN_1III_A_P_3")


def test_node_attributes_are_carried_not_parsed():
    _, g = _graph()
    a = g.nodes[N4]
    assert a["uniprot"] == "P00004" and a["entry_name"] == "NEI4_HUMAN"
    assert a["protein"] == "AF-P00004-F1" and a["pocket"] == "P_12_1"
    assert a["pocket_path"] == "04_pockets/P00004/P00004_af-p00004-f1.pdb_res_P_12_1.pdb"
    assert g.nodes[HUB2]["protein"] == "1HHH_C" and g.nodes[HUB2]["pocket"] == "P_3"


# --- ranking and hub choice ----------------------------------------------------------

def test_rank_hubs_columns_and_orders():
    _, g = _graph()
    w = sel.rank_hubs(g, "weighted")
    assert list(w.columns) == ["Node", "UniProt_ID", "Gene_Name", "Degree", "Weighted_Degree",
                               "Degree_Centrality", "Betweenness_Centrality",
                               "Closeness_Centrality"]
    assert w["Node"].iloc[0] == P73 and w["Degree"].iloc[0] == 5
    assert w["Node"].tolist()[1:3] == [HUB1, HUB2]       # tie at 4, weighted 3.45 > 3.20
    assert w["Weighted_Degree"].iloc[1] == pytest.approx(0.95 + 0.90 + 0.85 + 0.75)
    b = sel.rank_hubs(g, "betweenness")
    assert b["Degree"].is_monotonic_decreasing
    u = sel.rank_hubs(g, "union")
    assert u["Node"].tolist() == w["Node"].tolist()      # union ranks as weighted


def test_tied_top_skips_excluded():
    _, g = _graph()
    ranking = sel.rank_hubs(g, "weighted")
    assert sel.tied_top(ranking) == [P73]
    assert sel.tied_top(ranking, {"O15350"}) == [HUB1, HUB2]


def test_hub_skips_excluded_and_neighbours_exclude_it():
    t = _targets("weighted")
    hubs = t[t["Is_Hub"]]
    assert hubs["Node"].tolist() == [HUB1]               # O15350 has the top degree, skipped
    assert P73 not in set(t["Node"])                     # and it is a neighbour of the hub
    assert hubs["Similarity"].iloc[0] == 1.0
    sims = dict(zip(t["Node"], t["Similarity"]))
    assert sims == {HUB1: 1.0, N2: pytest.approx(0.90), N3: pytest.approx(0.85),
                    N4: pytest.approx(0.75)}


def test_without_exclusion_the_top_node_is_the_hub():
    t = _targets("weighted", exclude=())
    assert t[t["Is_Hub"]]["Node"].tolist() == [P73]
    assert HUB1 in set(t["Node"]) and len(t) == 6


def test_union_of_tied_hubs():
    t = _targets("union")
    assert t["Node"].is_unique
    assert set(t[t["Is_Hub"]]["Node"]) == {HUB1, HUB2}
    assert (t[t["Is_Hub"]]["Similarity"] == 1.0).all()
    sims = dict(zip(t["Node"], t["Similarity"]))
    assert set(sims) == {HUB1, HUB2, N2, N3, N4, N9}     # P73 excluded
    assert sims[N2] == pytest.approx(0.90)   # neighbour of both hubs: max(0.90, 0.78)
    assert sims[N3] == pytest.approx(0.85)   # max(0.85, 0.82)
    assert sims[N4] == pytest.approx(0.75)
    assert sims[N9] == pytest.approx(0.71)
    w = _targets("weighted")
    assert w["Is_Hub"].sum() == 1


def test_adjacent_tied_hubs_are_listed_once_as_hubs():
    import networkx as nx
    g = nx.Graph()
    for n in "ABCD":
        g.add_node(n, uniprot=n, entry_name=f"{n}_HUMAN", protein="1AAA_A", pocket="P_1",
                   pocket_path="x", structure_id="1AAA", chain="A")
    g.add_edge("A", "B", weight=0.8)
    g.add_edge("B", "C", weight=0.9)
    g.add_edge("A", "D", weight=0.75)
    g.add_edge("C", "D", weight=0.71)
    # A,B,C,D all degree 2 - make A and B the (adjacent) tied top by adding E
    g.add_node("E", uniprot="E", entry_name="E_HUMAN", protein="1AAA_A", pocket="P_1",
               pocket_path="x", structure_id="1AAA", chain="A")
    g.add_edge("A", "E", weight=0.7)
    g.add_edge("B", "E", weight=0.7)
    ranking = sel.rank_hubs(g, "union")
    assert sorted(sel.tied_top(ranking)) == ["A", "B"]
    t = sel.select_targets(g, ranking, set(), "union")
    assert t["Node"].is_unique
    assert t.set_index("Node").loc["B", "Is_Hub"]
    assert t.set_index("Node").loc["B", "Similarity"] == 1.0


def test_tie_break_deterministic():
    raw = small()
    expect = None
    for seed in range(10):
        order = list(raw.index)
        random.Random(seed).shuffle(order)
        shuffled = raw.loc[order]
        pairs, g = _graph(shuffled)
        for tb in ("weighted", "betweenness", "union"):
            ranking = sel.rank_hubs(g, tb)
            t = sel.select_targets(g, ranking, {"O15350"}, tb)
            got = (tb, tuple(sorted(t[t["Is_Hub"]]["Node"])), tuple(t["Node"]),
                   tuple(ranking["Node"]), tuple(pairs["Node"]), tuple(pairs["Node_2"]))
            if seed == 0:
                expect = expect or {}
                expect[tb] = got
            assert got == expect[tb], (seed, tb)
    assert expect["weighted"][1] == (HUB1,)


def test_no_non_excluded_node_raises():
    _, g = _graph()
    ranking = sel.rank_hubs(g, "weighted")
    every = set(ranking["UniProt_ID"])
    with pytest.raises(ValueError, match="no node"):
        sel.select_targets(g, ranking, every, "weighted")


# --- target table --------------------------------------------------------------------

def test_entry_name_starting_with_digits():
    t = _targets("union").set_index("Node")
    row = t.loc[N3]
    assert row["Gene_Name"] == "2AAA" and row["Organism"] == "HUMAN"
    assert row["Structure_ID"] == "3CCC"       # the real PDB ID, not the digits of the entry name
    assert row["Chain"] == "B" and row["Pocket"] == "P_1"
    af = t.loc[N4]
    assert af["Structure_ID"] == "AF-P00004-F1" and af["Pocket"] == "P_12_1"


def test_final_targets_columns_and_receptor_paths():
    t = _targets("union")
    assert list(t.columns) == ["Node", "UniProt_ID", "Gene_Name", "Organism", "Structure_ID",
                               "Similarity", "Protein", "Chain", "Pocket", "Pocket_Path",
                               "Receptor_Path", "Is_Hub"]
    t = t.set_index("Node")
    assert t.loc[HUB1, "Receptor_Path"] == "03_pdbfixer/P00001/1aaa_A.pdb"
    assert t.loc[HUB2, "Receptor_Path"] == "03_pdbfixer/P00008/1hhh_C.pdb"
    assert t.loc[N3, "Receptor_Path"] == "03_pdbfixer/P00003/3ccc_B.pdb"
    assert t.loc[N4, "Receptor_Path"] == "02_structures/P00004/AF-P00004-F1.pdb"
    assert t.loc[HUB1, "Pocket_Path"] == "04_pockets/P00001/P00001_1aaa_a.pdb_res_P_1.pdb"
    assert t["Similarity"].is_monotonic_decreasing


def test_receptor_path_comes_from_protein_info_not_the_protein_name():
    info = _protein_info()
    info.loc[info["UniProt ID"] == "P00001", "File"] = "02_structures/P00001/odd_name.pdb"
    _, g = _graph()
    ranking = sel.rank_hubs(g, "weighted")
    t = sel.select_targets(g, ranking, {"O15350"}, "weighted", protein_info=info)
    assert t.set_index("Node").loc[HUB1, "Receptor_Path"] == "03_pdbfixer/P00001/odd_name.pdb"


def test_missing_protein_info_row_is_a_clear_error():
    _, g = _graph()
    ranking = sel.rank_hubs(g, "weighted")
    with pytest.raises(ValueError, match="P00003.*3CCC.*B"):
        sel.select_targets(g, ranking, {"O15350"}, "weighted",
                           protein_info=_protein_info(extra_missing=("P00003",)))


# --- stage ---------------------------------------------------------------------------

def _run_folder(tmp_path, info=None):
    layout = RunLayout(tmp_path)
    layout.similarity.mkdir(parents=True)
    layout.pocketsim_csv.write_bytes(SMALL.read_bytes())
    layout.structures.mkdir(parents=True)
    (_protein_info() if info is None else info).to_csv(layout.protein_info, index=False)
    return layout


def _args(**kw):
    base = dict(pmax=0.7, exclude=["O15350"], hub_tie_break="union", nprocs=1, dry_run=False)
    base.update(kw)
    return SimpleNamespace(**base)


class Prog:
    def __init__(self):
        self.messages = []

    def advance(self, n=1, size=None):
        pass

    def notifier(self, title, message, priority=0):
        self.messages.append((title, message))


def test_run_writes_outputs_manifest_and_names_the_tie(tmp_path, capsys):
    layout = _run_folder(tmp_path)
    prog = Prog()
    counts = sel.run(layout, _args(), prog)
    d = layout.targets
    for name in ("final_targets.csv", "centrality.csv", "unique_pairs.csv", "manifest.json",
                 "failed.tsv"):
        assert (d / name).exists(), name
    assert (d / "failed.tsv").read_text() == "item\tstage\treason\n"
    final = pd.read_csv(layout.final_targets_csv)
    assert set(final["Node"]) == {HUB1, HUB2, N2, N3, N4, N9}
    assert final["Is_Hub"].sum() == 2
    assert pd.read_csv(d / "centrality.csv")["Node"].iloc[0] == P73
    assert len(pd.read_csv(d / "unique_pairs.csv")) == 12
    assert counts["failed"] == 0 and counts["targets"] == 6 and counts["hubs"] == 2
    m = json.loads((d / "manifest.json").read_text())
    assert m["stage"] == "select"
    assert m["args"]["exclude"] == {"O15350": "p73: problematic for docking (legacy decision)"}
    assert m["args"]["hub_tie_break"] == "union" and m["args"]["pmax"] == 0.7
    assert m["top_tie"] == [HUB1, HUB2] and m["hubs"] == [HUB1, HUB2]
    out = capsys.readouterr().out
    assert "tie" in out and HUB1 in out and HUB2 in out
    assert prog.messages and HUB2 in prog.messages[0][1]


def test_run_paths_are_relative_to_run_folder(tmp_path):
    layout = _run_folder(tmp_path)
    sel.run(layout, _args(), Prog())
    final = pd.read_csv(layout.final_targets_csv)
    for col in ("Pocket_Path", "Receptor_Path"):
        assert not final[col].str.startswith("/").any()


def test_rerun_on_same_inputs_is_skipped(tmp_path):
    layout = _run_folder(tmp_path)
    sel.run(layout, _args(), Prog())
    first = layout.targets.joinpath("final_targets.csv").stat().st_mtime_ns
    stamp = (layout.targets / "manifest.json").read_text()
    assert sel.units(layout, _args()).count == 0
    counts = sel.run(layout, _args(), Prog())
    assert counts["targets"] == 6 and counts["failed"] == 0
    assert layout.targets.joinpath("final_targets.csv").stat().st_mtime_ns == first
    assert (layout.targets / "manifest.json").read_text() == stamp
    # a different parameter is not skipped
    assert sel.units(layout, _args(hub_tie_break="weighted")).count == 1
    sel.run(layout, _args(hub_tie_break="weighted"), Prog())
    assert pd.read_csv(layout.final_targets_csv)["Is_Hub"].sum() == 1


def test_units_trivial_and_missing_input_clear_error(tmp_path):
    layout = RunLayout(tmp_path)
    with pytest.raises(FileNotFoundError, match="polyphar similarity"):
        sel.units(layout, _args())
    layout = _run_folder(tmp_path)
    wu = sel.units(layout, _args())
    assert wu.count == 1


def test_run_with_inconsistent_run_folder_raises_and_writes_no_outputs(tmp_path):
    layout = _run_folder(tmp_path, info=_protein_info(extra_missing=("P00003",)))
    with pytest.raises(ValueError, match="P00003"):
        sel.run(layout, _args(), Prog())
    assert not layout.final_targets_csv.exists()
    assert not (layout.targets / "manifest.json").exists()


def test_cli_dry_run_and_defaults(tmp_path, capsys):
    layout = _run_folder(tmp_path)
    args = cli.build_parser().parse_args(["select", "--run-dir", str(tmp_path)])
    assert args.pmax == 0.7 and args.exclude == ["O15350"] and args.hub_tie_break == "union"
    rc = cli.main(["select", "--run-dir", str(tmp_path), "--dry-run"])
    assert rc == 0 and "select: 1 units" in capsys.readouterr().out
    assert not layout.targets.exists()


def test_cli_end_to_end(tmp_path, monkeypatch):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    layout = _run_folder(tmp_path / "run")
    rc = cli.main(["select", "--run-dir", str(tmp_path / "run"), "--no-notify",
                   "--hub-tie-break", "weighted", "--exclude", "O15350", "P00002"])
    assert rc == 0
    final = pd.read_csv(layout.final_targets_csv)
    assert "P00002" not in set(final["UniProt_ID"])
    assert final["Is_Hub"].sum() == 1


# --- regression: 2023 target selection -----------------------------------------------

REFERENCE = Path("/home/onur/repos/Gene_Disease_Association/target_selection_results")
SUBTYPES = {   # reference stem -> path under $ADEP_DATA/disease_subtypes
    "liver_cancer_net2": "Liver_cancer/net2/PM_net2/Pocketsim_net2.csv",
    "liver_cancer_net3": "Liver_cancer/net3/PM_net3/Pocketsim_net3.csv",
    "alzheimer_Cluster1": "Alzhemier/Cluster1/PM_Cluster1/Pocketsim.csv",
    "alzheimer_Cluster3": "Alzhemier/Cluster3/PM_Cluster3/Pocketsim.csv",
    "alzheimer_Cluster4": "Alzhemier/Cluster4/PM_Cluster4/Pocketsim.csv",
}


def _legacy_hub(stem):
    """The 2023 hub: first row of its centrality.csv, the second if that is O15350."""
    c = pd.read_csv(REFERENCE / f"{stem}_centrality.csv")
    return c["Node"].iloc[1 if c["UniProt_ID"].iloc[0] == "O15350" else 0]


def _select(df, hub, ordered=False):
    pairs = sel.unique_protein_pairs(df, 0.7, ordered=ordered)
    g = sel.build_graph(pairs)
    ranking = sel.rank_hubs(g, "weighted")
    return sel.select_targets(g, ranking, {"O15350"}, "weighted", hub=hub), g, ranking


@pytest.mark.regression
@pytest.mark.parametrize("stem", list(SUBTYPES))
def test_reproduces_2023_targets(stem):
    path = Path(os.environ["ADEP_DATA"]) / "disease_subtypes" / SUBTYPES[stem]
    df = sel.read_pocketsim(path)
    ref = pd.read_csv(REFERENCE / f"{stem}_final_targets.csv")
    hub = _legacy_hub(stem)
    got, g, ranking = _select(df, hub)
    want = dict(zip(ref["Node"], ref["Similarity"].round(6)))
    have = dict(zip(got["Node"], got["Similarity"].round(6)))
    if have != want:
        only_ref, only_new = sorted(set(want) - set(have)), sorted(set(have) - set(want))
        diff_sim = sorted(n for n in set(want) & set(have) if want[n] != have[n])
        legacy, _, _ = _select(df, hub, ordered=True)
        legacy_have = dict(zip(legacy["Node"], legacy["Similarity"].round(6)))
        print(f"\n{stem}: hub {hub}\n  only in 2023 reference: {only_ref}\n"
              f"  only in new result: {only_new}\n  similarity differs: {diff_sim}\n"
              f"  ordered-pair dedup reproduces the reference: {legacy_have == want}\n"
              f"  top degree tie: {sel.tied_top(ranking, {'O15350'})}")
    assert have == want
