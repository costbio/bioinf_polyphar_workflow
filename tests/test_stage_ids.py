import json
import os
from pathlib import Path

import pandas as pd
import pytest

from polyphar import cli
from polyphar.common import RunLayout, sha256_file
from polyphar.stages.ids import build_id_list

ROWS = [("P1", "reviewed"), ("P2", "unreviewed"), (" P1 ", "Reviewed"), ("P3;P4", "reviewed")]


def _xlsx(path, rows=ROWS, sheet="net1"):
    df = pd.DataFrame(rows, columns=["uniprot Entry", "Reviewed/unreviewed"])
    with pd.ExcelWriter(path) as xw:
        df.to_excel(xw, sheet_name=sheet, index=False)
    return path


def test_dedup_reviewed_and_order(tmp_path):
    x = _xlsx(tmp_path / "w.xlsx")
    assert build_id_list(x, "net1", "uniprot Entry", "Reviewed/unreviewed") == ["P1", "P3", "P4"]


def test_no_reviewed_column_keeps_all_rows(tmp_path):
    x = _xlsx(tmp_path / "w.xlsx")
    assert build_id_list(x, "net1", "uniprot Entry", None) == ["P1", "P2", "P3", "P4"]


def test_split_on_comma_and_whitespace_and_blank_cells(tmp_path):
    x = _xlsx(tmp_path / "w.xlsx", rows=[("A1, A2  A3", "reviewed"), (None, "reviewed"), ("A2;;A4", "reviewed")])
    assert build_id_list(x, "net1", "uniprot Entry", "Reviewed/unreviewed") == ["A1", "A2", "A3", "A4"]


def test_missing_sheet_names_available_sheets(tmp_path):
    x = _xlsx(tmp_path / "w.xlsx")
    with pytest.raises(ValueError, match=r"nope.*net1"):
        build_id_list(x, "nope", "uniprot Entry", None)


def test_missing_column_names_available_columns(tmp_path):
    x = _xlsx(tmp_path / "w.xlsx")
    with pytest.raises(ValueError, match=r"'bogus'.*uniprot Entry"):
        build_id_list(x, "net1", "bogus", None)
    with pytest.raises(ValueError, match=r"'bogus'.*Reviewed/unreviewed"):
        build_id_list(x, "net1", "uniprot Entry", "bogus")


def test_cli_xlsx_writes_ids_and_manifest(tmp_path, monkeypatch):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    run = tmp_path / "run"
    (run / "inputs").mkdir(parents=True)
    x = _xlsx(run / "inputs" / "w.xlsx")
    rc = cli.main(["ids", "--run-dir", str(run), "--no-notify", "--xlsx", "inputs/w.xlsx",
                   "--sheet", "net1", "--id-column", "uniprot Entry",
                   "--reviewed-column", "Reviewed/unreviewed"])
    assert rc == 0
    layout = RunLayout(run)
    assert layout.ids_txt.read_text() == "P1\nP3\nP4\n"
    m = json.loads((layout.stage_dir("00_ids") / "manifest.json").read_text())
    assert m["stage"] == "ids"
    assert m["args"]["sheet"] == "net1"
    assert m["args"]["id_column"] == "uniprot Entry"
    assert m["args"]["reviewed_column"] == "Reviewed/unreviewed"
    assert m["inputs"][0]["sha256"] == sha256_file(x)
    assert m["counts"] == {"rows_read": 4, "ids_kept": 3, "duplicates_dropped": 1,
                           "unreviewed_dropped": 1, "failed": 0}


def test_cli_txt_mode_same_normalisation_no_reviewed_filter(tmp_path, monkeypatch):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    run = tmp_path / "run"
    (run / "inputs").mkdir(parents=True)
    t = run / "inputs" / "ids.txt"
    t.write_bytes(b"P1\r\n P2;P3 \r\nP1\r\n\r\n")
    assert cli.main(["ids", "--run-dir", str(run), "--no-notify", "--txt", "inputs/ids.txt"]) == 0
    layout = RunLayout(run)
    assert layout.ids_txt.read_text() == "P1\nP2\nP3\n"
    m = json.loads((layout.stage_dir("00_ids") / "manifest.json").read_text())
    assert m["args"]["txt"] == "inputs/ids.txt"
    assert m["inputs"][0]["sha256"] == sha256_file(t)
    assert m["counts"]["ids_kept"] == 3 and m["counts"]["duplicates_dropped"] == 1
    assert m["counts"]["failed"] == 0


def test_cli_xlsx_and_txt_are_mutually_exclusive(tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["ids", "--run-dir", str(tmp_path), "--xlsx", "a.xlsx", "--txt", "a.txt"])


def test_cli_requires_a_source(tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["ids", "--run-dir", str(tmp_path)])


def test_cli_xlsx_needs_sheet_and_id_column(tmp_path, monkeypatch):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    with pytest.raises(SystemExit, match="--sheet"):
        cli.main(["ids", "--run-dir", str(tmp_path), "--no-notify", "--xlsx", "a.xlsx"])


def test_stage_add_arguments_hook(monkeypatch):
    """A stage module's optional add_arguments(parser) extends its own subcommand only."""
    import sys
    import types
    mod = types.ModuleType("fake_args_stage")
    mod.add_arguments = lambda p: p.add_argument("--frobnicate", default="x")
    monkeypatch.setitem(sys.modules, "fake_args_stage", mod)
    monkeypatch.setitem(cli.STAGES, "fix", ("fake_args_stage", "fake"))
    parser = cli.build_parser()
    assert parser.parse_args(["fix", "--frobnicate", "y"]).frobnicate == "y"
    with pytest.raises(SystemExit):
        parser.parse_args(["uniprot", "--frobnicate", "y"])


# ---------------------------------------------------------------- regression
HCC = Path(os.environ.get("ADEP_DATA", "")) / "disease_subtypes" / "Liver_cancer"


@pytest.fixture(scope="module")
def hcc_ids(tmp_path_factory):
    """Read the (slow, /mnt/c) workbook once; build all three tabs from a local copy."""
    import shutil
    local = tmp_path_factory.mktemp("hcc") / "w.xlsx"
    shutil.copy(HCC / "clusters_features_ppi.xlsx", local)
    return {t: build_id_list(local, t, "uniprot Entry", "Reviewed/unreviewed")
            for t in ("net1", "net2", "net3")}


def _ref(rel):
    return {l.strip() for l in (HCC / rel).read_text().splitlines() if l.strip()}


@pytest.mark.regression
def test_hcc_tabs_match_known_counts(hcc_ids):
    assert len(hcc_ids["net1"]) == 459
    assert set(hcc_ids["net1"]) != _ref("net1/net1UniprotID.txt")
    assert set(hcc_ids["net2"]) == _ref("net2/net2.txt") and len(hcc_ids["net2"]) == 156
    assert set(hcc_ids["net3"]) == _ref("net3/net3.txt") and len(hcc_ids["net3"]) == 197
