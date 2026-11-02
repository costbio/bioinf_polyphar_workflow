"""Item 1 of the final review: `run` must not skip a 0-unit stage whose inputs changed.

Real structures / fix / pockets stages (offline: canned SIFTS / RCSB / AFDB responses, the
fix outputs and DoGSiteScorer outputs already on disk, so every stage has 0 units on the
second run). After uniprot.csv shrinks, the removed accession must disappear from
Protein_info.csv and DrugScore.csv. (similarity and select are fingerprinted on DrugScore /
Pocketsim, see their own resume tests.)
"""
import pandas as pd

from polyphar import cli
from polyphar.common import RunLayout
from test_stage_structures import AF_JSON, FakeSession, _cif_2mjw

ROWS = [("P25815", "S100P_HUMAN", "2MJW"), ("Q15109", "RAGE_HUMAN", "2MJW")]


def _uniprot(layout, rows):
    layout.uniprot_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([(a, n, n.split("_")[0], "reviewed", "9606", p) for a, n, p in rows],
                 columns=["Entry", "Entry Name", "Gene Names (primary)", "Reviewed",
                          "Organism (ID)", "PDB"]).to_csv(layout.uniprot_csv, index=False)


def _outputs(layout, acc, name):
    """A fixed chain and DoGSite3 outputs for it (two pockets, 400 and 250 A^3)."""
    fixed = layout.pdbfixer / acc / f"{name}.pdb"
    fixed.parent.mkdir(parents=True, exist_ok=True)
    fixed.write_text("ATOM fixed\n")
    d = layout.pockets / acc
    d.mkdir(parents=True, exist_ok=True)
    for out in (name, f"AF-{acc}-F1"):
        (d / f"{out}_desc.txt").write_text("name\tvolume\nP_1\t400\nP_2\t250\n")
        (d / f"{out}_P_1_res.pdb").write_text("ATOM\n")
        (d / f"{out}_P_2_res.pdb").write_text("ATOM\n")


def test_shrunken_input_reaches_protein_info_and_drugscore(tmp_path, monkeypatch):
    from polyphar.stages import run as run_mod
    from polyphar.stages import structures

    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    monkeypatch.setattr(run_mod, "_needs_benchmark", lambda stage: False)
    af = AF_JSON.read_text()
    session = FakeSession(cifs={"2MJW": _cif_2mjw(tmp_path).read_text()},
                          af={"P25815": af, "Q15109": af.replace("P25815", "Q15109")})
    monkeypatch.setattr(structures, "http_session", lambda: session)
    layout = RunLayout(tmp_path / "run")
    _uniprot(layout, ROWS)
    _outputs(layout, "P25815", "2mjw_B")
    _outputs(layout, "Q15109", "2mjw_A")
    argv = ["run", "--run-dir", str(layout.run_dir), "--from", "structures", "--to", "pockets",
            "--no-notify", "--nprocs", "1"]

    assert cli.main(argv) == 0
    info = pd.read_csv(layout.protein_info, dtype=str, keep_default_na=False)
    assert set(info["UniProt ID"]) == {"P25815", "Q15109"}
    assert set(pd.read_csv(layout.drugscore_csv, index_col=0)["UniProt ID"]) == {"P25815",
                                                                                 "Q15109"}
    _uniprot(layout, ROWS[:1])   # a corrected, smaller input; nothing left to download
    for stage in (structures,):
        assert stage.units(layout, None).count == 0
    assert cli.main(argv) == 0
    info = pd.read_csv(layout.protein_info, dtype=str, keep_default_na=False)
    assert set(info["UniProt ID"]) == {"P25815"}
    assert set(pd.read_csv(layout.drugscore_csv, index_col=0)["UniProt ID"]) == {"P25815"}
    assert "Q15109" not in (layout.structures / "alphafold.tsv").read_text()

    # a changed pocket selection reruns pockets although no model is left to run
    selected = pd.read_csv(layout.drugscore_csv, index_col=0)["Selected"]
    assert list(selected) == [True, False, True, False]
    assert cli.main([*argv, "--from", "pockets", "--min-volume", "200"]) == 0
    assert pd.read_csv(layout.drugscore_csv, index_col=0)["Selected"].all()

