"""Tests for the panel-based select stage.

The hub-and-neighbours tests this file used to hold are gone with the mechanism. The
filters themselves are tested in test_panel.py; here the concern is the stage: schema
compatibility with `dock`, receptor path resolution, exclusions, the funnel, and the
outputs.
"""
import csv
import json
import types
from pathlib import Path

import pytest

from polyphar.common import RunLayout
from polyphar.stages import select
from polyphar.stages.annotate import encode_residues


def _args(**kw):
    base = dict(rg_max=8.5, volume_min=300.0, volume_max=2000.0, drugscore_min=0.5,
                plddt_min=70.0, tm_policy="C", tm_overlap=0.10, exclude=["O15350"],
                no_fetch_tm=True)
    base.update(kw)
    return types.SimpleNamespace(**base)


class _Progress:
    def advance(self, n, size=0.0):
        pass


def _make_run(tmp_path, drugscore_rows, annotation_rows, protein_info_rows, tm=None):
    """A realistic run folder. tm defaults to "no TM segments" for every accession the
    DrugScore rows mention, as a real tm.tsv would have, so no test touches the network."""
    if tm is None:
        tm = "".join(f"{a}\t0\t\n"
                     for a in sorted({r["UniProt ID"] for r in drugscore_rows}))
    run = tmp_path / "run"
    (run / "04_pockets" / "annotation").mkdir(parents=True)
    (run / "02_structures").mkdir(parents=True)
    (run / "01_uniprot").mkdir(parents=True)

    def write(path, rows, cols):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)

    write(run / "04_pockets" / "DrugScore.csv", drugscore_rows,
          ["UniProt ID", "Gene Name", "Protein", "Pocket", "Pocket Filename",
           "volume", "drugScore", "Selected"])
    write(run / "04_pockets" / "annotation" / "pocket_annotation.csv", annotation_rows,
          ["Pocket Filename", "AlphaFold", "radius_of_gyration", "pocket_plddt_mean",
           "box_size_x", "box_size_y", "box_size_z", "residues"])
    write(run / "02_structures" / "Protein_info.csv", protein_info_rows,
          ["UniProt ID", "PDB", "Chain ID", "AlphaFold", "File"])
    # pre-created so tmfetch never reaches the network during tests
    (run / "01_uniprot" / "tm.tsv").write_text("accession\ttm_count\tsegments\n" + tm)
    return RunLayout(run)


def _ds(acc, protein, pocket, fname, vol=800.0, ds=0.9):
    return {"UniProt ID": acc, "Gene Name": f"{acc}G_HUMAN", "Protein": protein,
            "Pocket": pocket, "Pocket Filename": fname, "volume": str(vol),
            "drugScore": str(ds), "Selected": "True"}


def _an(fname, af=True, rg=6.0, plddt=90.0, residues=(10, 11, 12)):
    return {"Pocket Filename": fname, "AlphaFold": str(af),
            "radius_of_gyration": str(rg),
            "pocket_plddt_mean": "" if plddt is None else str(plddt),
            "box_size_x": "18", "box_size_y": "18", "box_size_z": "18",
            "residues": encode_residues(list(residues))}


# --- helpers that survived the rewrite ------------------------------------------------

def test_split_protein_handles_alphafold_and_pdb_chains():
    assert select.split_protein("2MJW_B") == ("2MJW", "B")
    assert select.split_protein("AF-P06493-F1") == ("AF-P06493-F1", "A")
    assert select.split_protein("1ABC") == ("1ABC", select.NO_CHAIN)


def test_split_entry_name_splits_at_the_last_underscore():
    assert select.split_entry_name("S100P_HUMAN") == ("S100P", "HUMAN")
    assert select.split_entry_name("NOUNDERSCORE") == ("NOUNDERSCORE", "")


# --- receptor paths -------------------------------------------------------------------

def test_receptor_path_for_alphafold_uses_the_model_file():
    panel = {"P1": (_ds("P1", "AF-P1-F1", "P_1", "f.pdb"), _an("f.pdb"))}
    info = [{"UniProt ID": "P1", "PDB": "", "Chain ID": "A", "AlphaFold": "True",
             "File": "02_structures/P1/AF-P1-F1.pdb"}]
    paths, chains = select._receptor_paths(panel, info)
    assert paths["P1"] == "02_structures/P1/AF-P1-F1.pdb"
    assert chains["P1"] == "A"


def test_receptor_path_for_a_pdb_chain_points_at_pdbfixer_and_keeps_the_real_chain():
    """The Protein field upper-cases the chain, so the lookup must be case-insensitive
    while Chain reports the real author chain."""
    panel = {"P1": (_ds("P1", "2MJW_B", "P_1", "f.pdb"), _an("f.pdb", af=False))}
    info = [{"UniProt ID": "P1", "PDB": "2MJW", "Chain ID": "b", "AlphaFold": "False",
             "File": "02_structures/P1/2mjw_b.pdb"}]
    paths, chains = select._receptor_paths(panel, info)
    assert paths["P1"] == "03_pdbfixer/P1/2mjw_b.pdb"
    assert chains["P1"] == "b"


def test_chains_differing_only_by_case_are_ambiguous_not_guessed():
    panel = {"P1": (_ds("P1", "2MJW_B", "P_1", "f.pdb"), _an("f.pdb", af=False))}
    info = [{"UniProt ID": "P1", "PDB": "2MJW", "Chain ID": "b", "AlphaFold": "False",
             "File": "02_structures/P1/2mjw_b.pdb"},
            {"UniProt ID": "P1", "PDB": "2MJW", "Chain ID": "B", "AlphaFold": "False",
             "File": "02_structures/P1/2mjw_B.pdb"}]
    with pytest.raises(ValueError, match="ambiguous"):
        select._receptor_paths(panel, info)


def test_a_target_with_no_protein_info_row_is_an_error():
    panel = {"P1": (_ds("P1", "9XYZ_A", "P_1", "f.pdb"), _an("f.pdb", af=False))}
    with pytest.raises(ValueError, match="no row in Protein_info"):
        select._receptor_paths(panel, [])


# --- output schema ----------------------------------------------------------------------

def test_target_rows_keep_the_columns_dock_consumes():
    """dock reads Protein, Pocket, Receptor_Path and Pocket_Path; the schema must not
    drift even though Similarity and Is_Hub are now meaningless."""
    panel = {"P1": (_ds("P1", "AF-P1-F1", "P_2", "AF-P1-F1_res_P_2.pdb"), _an("x"))}
    rows = select.target_rows(panel, {"P1": "02_structures/P1/AF-P1-F1.pdb"}, {"P1": "A"})
    assert list(rows[0]) == select.TARGET_COLUMNS
    r = rows[0]
    assert r["Protein"] == "AF-P1-F1" and r["Pocket"] == "P_2"
    assert r["Pocket_Path"] == "04_pockets/P1/AF-P1-F1_res_P_2.pdb"
    assert r["Receptor_Path"] == "02_structures/P1/AF-P1-F1.pdb"
    assert r["Gene_Name"] == "P1G" and r["Organism"] == "HUMAN"
    assert r["Similarity"] == "" and r["Is_Hub"] is False


def test_target_rows_are_sorted_by_accession_not_input_order():
    panel = {
        "P2": (_ds("P2", "AF-P2-F1", "P_1", "b.pdb"), _an("b.pdb")),
        "P1": (_ds("P1", "AF-P1-F1", "P_1", "a.pdb"), _an("a.pdb")),
    }
    rows = select.target_rows(panel, {"P1": "x", "P2": "y"}, {"P1": "A", "P2": "A"})
    assert [r["UniProt_ID"] for r in rows] == ["P1", "P2"]


# --- the stage -----------------------------------------------------------------------

def test_run_writes_the_panel_the_funnel_and_a_reason_per_rejection(tmp_path):
    ds_rows = [
        _ds("P1", "AF-P1-F1", "P_1", "a.pdb"),                 # kept
        _ds("P2", "AF-P2-F1", "P_1", "b.pdb", vol=5000.0),     # rejected: volume
        _ds("P3", "AF-P3-F1", "P_1", "c.pdb"),                 # rejected: rg
    ]
    anns = [_an("a.pdb"), _an("b.pdb"), _an("c.pdb", rg=20.0)]
    info = [{"UniProt ID": a, "PDB": "", "Chain ID": "A", "AlphaFold": "True",
             "File": f"02_structures/{a}/AF-{a}-F1.pdb"} for a in ("P1", "P2", "P3")]
    layout = _make_run(tmp_path, ds_rows, anns, info)

    counts = select.run(layout, _args(), _Progress())

    assert counts["panel"] == 1 and counts["selected_pockets"] == 3
    targets = list(csv.DictReader(open(layout.targets / "final_targets.csv")))
    assert [t["UniProt_ID"] for t in targets] == ["P1"]

    funnel = list(csv.DictReader(open(layout.targets / "panel_funnel.tsv")))
    counts_seq = [int(f["pockets_remaining"]) for f in funnel]
    assert counts_seq == sorted(counts_seq, reverse=True)      # monotone chain
    assert counts_seq[0] == 3 and counts_seq[-1] == 1

    failed = (layout.targets / "failed.tsv").read_text().strip().splitlines()
    assert len(failed) == 3                                     # header + 2 rejections
    assert any("volume" in l for l in failed) and any("gyration" in l for l in failed)

    manifest = json.loads((layout.targets / "manifest.json").read_text())
    assert manifest["counts"]["panel"] == 1
    assert manifest["args"]["rg_max"] == 8.5


def test_excluded_accessions_are_dropped_and_recorded(tmp_path):
    ds_rows = [_ds("O15350", "AF-O15350-F1", "P_1", "a.pdb"),
               _ds("P1", "AF-P1-F1", "P_1", "b.pdb")]
    anns = [_an("a.pdb"), _an("b.pdb")]
    info = [{"UniProt ID": a, "PDB": "", "Chain ID": "A", "AlphaFold": "True",
             "File": f"02_structures/{a}/AF-{a}-F1.pdb"} for a in ("O15350", "P1")]
    layout = _make_run(tmp_path, ds_rows, anns, info)

    counts = select.run(layout, _args(), _Progress())

    assert counts["panel"] == 1 and counts["excluded"] == 1
    targets = list(csv.DictReader(open(layout.targets / "final_targets.csv")))
    assert [t["UniProt_ID"] for t in targets] == ["P1"]
    assert "p73" in (layout.targets / "failed.tsv").read_text()


def test_missing_annotation_file_is_a_clear_error(tmp_path):
    layout = _make_run(tmp_path, [_ds("P1", "AF-P1-F1", "P_1", "a.pdb")], [_an("a.pdb")],
                       [{"UniProt ID": "P1", "PDB": "", "Chain ID": "A",
                         "AlphaFold": "True", "File": "x.pdb"}])
    (layout.pockets / "annotation" / "pocket_annotation.csv").unlink()
    with pytest.raises(FileNotFoundError, match="annotate"):
        select.run(layout, _args(), _Progress())


def test_tm_policy_is_recorded_in_the_manifest(tmp_path):
    """The policy changes which pockets are eligible, so a run is not interpretable
    without it."""
    layout = _make_run(tmp_path, [_ds("P1", "AF-P1-F1", "P_1", "a.pdb")], [_an("a.pdb")],
                       [{"UniProt ID": "P1", "PDB": "", "Chain ID": "A",
                         "AlphaFold": "True", "File": "02_structures/P1/AF-P1-F1.pdb"}])
    select.run(layout, _args(tm_policy="B"), _Progress())
    manifest = json.loads((layout.targets / "manifest.json").read_text())
    assert manifest["args"]["tm_policy"] == "B"


def test_offline_mode_fails_clearly_when_tm_data_is_missing(tmp_path):
    """--no-fetch-tm must not silently skip the transmembrane filter."""
    layout = _make_run(tmp_path, [_ds("P1", "AF-P1-F1", "P_1", "a.pdb")], [_an("a.pdb")],
                       [{"UniProt ID": "P1", "PDB": "", "Chain ID": "A",
                         "AlphaFold": "True", "File": "02_structures/P1/AF-P1-F1.pdb"}],
                       tm="")        # tm.tsv exists but is empty
    with pytest.raises(FileNotFoundError, match="no-fetch-tm"):
        select.run(layout, _args(no_fetch_tm=True), _Progress())
