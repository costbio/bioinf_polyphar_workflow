import json
import math
import re
from pathlib import Path

import pandas as pd
import pytest
import requests

from polyphar import cli
from polyphar.common import RunLayout, sha256_file
from polyphar.stages import uniprot
from polyphar.stages.uniprot import classify, fetch_uniprot, parse_uniprot_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "uniprot_sample.tsv"


def _sample():
    return parse_uniprot_tsv(FIXTURE.read_text())


class FakeResponse:
    def __init__(self, text="", status=200, release="2026_03"):
        self.text = text
        self.status_code = status
        self.headers = {"x-uniprot-release": release} if release else {}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(f"{self.status_code} error", response=self)


class FakeSession:
    """Answers `accession:X` queries from a table {accession: tsv row (no newline)}.
    A secondary accession maps to the row of its primary."""

    HEADER = "Entry\tEntry Name\tGene Names (primary)\tReviewed\tOrganism (ID)\tPDB"

    def __init__(self, table, release="2026_03", fail_on=None):
        self.table = table
        self.release = release
        self.fail_on = fail_on
        self.calls = []

    def get(self, url, params=None, timeout=None):
        assert timeout, "every request must carry a timeout"
        self.calls.append(params["query"])
        accs = re.findall(r"accession:(\w+)", params["query"])
        if self.fail_on and self.fail_on in accs:
            raise requests.exceptions.RetryError("too many 503 responses")
        rows = []
        for a in accs:
            row = self.table.get(a)
            if row and row not in rows:
                rows.append(row)
        fields = params["fields"].split(",")
        text = self.HEADER + "\n" if len(fields) > 1 else "Entry\n"
        for r in rows:
            text += (r if len(fields) > 1 else r.split("\t")[0]) + "\n"
        return FakeResponse(text, release=self.release)


ROWS = {
    "P25815": "P25815\tS100P_HUMAN\tS100P\treviewed\t9606\t1J55;1OZO;2MJW;7NMI;",
    "Q9NPD8": "Q9NPD8\tUBE2T_HUMAN\tUBE2T\treviewed\t9606\t1YH2;",
    "UNREV1": "UNREV1\tX_HUMAN\tX\tunreviewed\t9606\t",
    "ECOLI1": "ECOLI1\tRECA_ECOLI\trecA\treviewed\t562\t1REA;",
}
ROWS["OLD1"] = ROWS["P25815"]   # secondary accession -> primary row


def test_parse_strips_trailing_semicolon():
    df = _sample()
    pdb = df.set_index("Entry")["PDB"]
    assert pdb["P25815"] == "1J55;1OZO;2MJW;7NMI"
    assert list(df.columns) == ["Entry", "Entry Name", "Gene Names (primary)", "Reviewed",
                                "Organism (ID)", "PDB"]
    assert pdb["A0A024RBG1"] == ""          # no structures: empty, not "nan"
    assert df["Organism (ID)"].tolist() == [9606, 9606, 9606, 562]


def test_parse_empty_response_gives_empty_frame_with_columns():
    for text in ("", "Entry\tEntry Name\tGene Names (primary)\tReviewed\tOrganism (ID)\tPDB\n"):
        df = parse_uniprot_tsv(text)
        assert df.empty
        assert "Entry" in df.columns and "PDB" in df.columns


def test_classify_splits_categories():
    df = _sample()
    requested = ["P25815", "Q9NPD8", "A0A024RBG1", "P0A7G6", "NOPE01"]
    res = classify(requested, df)
    assert res["kept"]["Entry"].tolist() == ["P25815", "Q9NPD8"]
    assert res["not_found"] == ["NOPE01"]
    assert res["unreviewed"] == ["A0A024RBG1"]
    assert res["non_human"] == ["P0A7G6"]
    assert res["renamed"] == {}


def test_classify_reports_renamed():
    df = _sample()
    seen = []

    def requery(acc):
        seen.append(acc)
        return {"OLD1": "P25815"}.get(acc)

    res = classify(["OLD1", "Q9NPD8", "GONE01"], df, requery=requery)
    assert seen == ["OLD1", "GONE01"]            # only the missing ones, one at a time
    assert res["renamed"] == {"OLD1": "P25815"}
    assert res["not_found"] == ["GONE01"]
    assert "P25815" in res["kept"]["Entry"].tolist()


def test_classify_without_requery_calls_missing_not_found():
    res = classify(["OLD1"], _sample())
    assert res["not_found"] == ["OLD1"] and res["renamed"] == {}


def test_fetch_batches_dedupes_and_returns_release():
    s = FakeSession(ROWS)
    df, release = fetch_uniprot(["P25815", "OLD1", "Q9NPD8", "UNREV1", "ECOLI1"], s, batch=2)
    assert len(s.calls) == 3
    assert release == "2026_03"
    assert df["Entry"].tolist() == ["P25815", "Q9NPD8", "UNREV1", "ECOLI1"]  # OLD1 dup dropped
    assert df.set_index("Entry")["PDB"]["P25815"] == "1J55;1OZO;2MJW;7NMI"


def test_fetch_batch_failure_is_loud_and_names_the_batch():
    s = FakeSession(ROWS, fail_on="Q9NPD8")
    with pytest.raises(RuntimeError, match=r"batch 2.*Q9NPD8"):
        fetch_uniprot(["P25815", "Q9NPD8"], s, batch=1)


def test_fetch_http_error_includes_server_message():
    class Bad(FakeSession):
        def get(self, url, params=None, timeout=None):
            return FakeResponse("Error messages\nThe 'accession' filter value 'ZZ' has invalid format",
                                status=400)
    with pytest.raises(RuntimeError, match="invalid format"):
        fetch_uniprot(["ZZ"], Bad({}), batch=100)


def test_fetch_release_changing_mid_run_is_an_error():
    class Drift(FakeSession):
        def get(self, url, params=None, timeout=None):
            r = super().get(url, params, timeout)
            r.headers["x-uniprot-release"] = f"2026_0{len(self.calls)}"
            return r
    with pytest.raises(RuntimeError, match="release"):
        fetch_uniprot(["P25815", "Q9NPD8"], Drift(ROWS), batch=1)


def test_units_is_ceil_of_ids_over_100(tmp_path):
    layout = RunLayout(tmp_path)
    layout.ids_txt.parent.mkdir(parents=True)
    layout.ids_txt.write_text("".join(f"P{i:05d}\n" for i in range(201)))
    wu = uniprot.units(layout, None)
    assert (wu.count, wu.cost_basis) == (math.ceil(201 / 100), "batch")


def test_units_missing_ids_is_a_clear_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="polyphar ids"):
        uniprot.units(RunLayout(tmp_path), None)


def _setup_run(tmp_path, monkeypatch, ids, session):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    monkeypatch.setattr(uniprot, "http_session", lambda *a, **k: session)
    run = tmp_path / "run"
    layout = RunLayout(run)
    layout.ids_txt.parent.mkdir(parents=True)
    layout.ids_txt.write_text("".join(f"{i}\n" for i in ids))
    return run, layout


def test_run_writes_all_outputs_and_manifest(tmp_path, monkeypatch):
    s = FakeSession(ROWS)
    ids = ["P25815", "OLD1", "UNREV1", "ECOLI1", "Q9NPD8", "GONE01"]
    run, layout = _setup_run(tmp_path, monkeypatch, ids, s)
    assert cli.main(["uniprot", "--run-dir", str(run), "--no-notify"]) == 0
    d = layout.stage_dir("01_uniprot")
    kept = pd.read_csv(layout.uniprot_csv, dtype=str)
    assert kept["Entry"].tolist() == ["P25815", "Q9NPD8"]
    assert kept["PDB"].tolist()[0] == "1J55;1OZO;2MJW;7NMI"
    assert (d / "not_found.txt").read_text() == "GONE01\n"
    assert (d / "renamed.tsv").read_text() == "requested\tprimary\nOLD1\tP25815\n"
    assert (d / "excluded.tsv").read_text() == (
        "accession\treason\nUNREV1\tunreviewed\nECOLI1\tnon_human\n")
    assert (d / "failed.tsv").read_text() == "item\tstage\treason\n"
    m = json.loads((d / "manifest.json").read_text())
    assert m["stage"] == "uniprot"
    assert m["versions"] == {"uniprot": "2026_03"}
    assert m["inputs"][0]["sha256"] == sha256_file(layout.ids_txt)
    assert m["counts"] == {"requested": 6, "kept": 2, "not_found": 1, "unreviewed": 1,
                           "non_human": 1, "renamed": 1, "failed": 0}


def test_run_is_resumable_and_refetches_when_ids_change(tmp_path, monkeypatch):
    s = FakeSession(ROWS)
    run, layout = _setup_run(tmp_path, monkeypatch, ["P25815", "Q9NPD8"], s)
    assert cli.main(["uniprot", "--run-dir", str(run), "--no-notify"]) == 0
    n = len(s.calls)
    assert uniprot.units(layout, None).count == 0
    assert cli.main(["uniprot", "--run-dir", str(run), "--no-notify"]) == 0
    assert len(s.calls) == n                                  # no network on the rerun
    m = json.loads((layout.stage_dir("01_uniprot") / "manifest.json").read_text())
    assert m["counts"]["kept"] == 2
    layout.ids_txt.write_text("P25815\n")                    # ids changed -> refetch
    assert uniprot.units(layout, None).count == 1
    assert cli.main(["uniprot", "--run-dir", str(run), "--no-notify"]) == 0
    assert len(s.calls) > n
    assert pd.read_csv(layout.uniprot_csv, dtype=str)["Entry"].tolist() == ["P25815"]


def test_run_batch_failure_raises_and_leaves_no_manifest(tmp_path, monkeypatch):
    s = FakeSession(ROWS, fail_on="Q9NPD8")
    run, layout = _setup_run(tmp_path, monkeypatch, ["P25815", "Q9NPD8"], s)
    with pytest.raises(RuntimeError):
        cli.main(["uniprot", "--run-dir", str(run), "--no-notify"])
    assert not (layout.stage_dir("01_uniprot") / "manifest.json").exists()
    assert not layout.uniprot_csv.exists()


@pytest.mark.network
def test_fetch_live():
    from polyphar.common import http_session
    df, release = fetch_uniprot(["P25815", "Q9NPD8"], http_session())
    assert re.fullmatch(r"\d{4}_\d{2}", release)
    assert "2MJW" in df.set_index("Entry").loc["P25815", "PDB"].split(";")
