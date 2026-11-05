import json
import os
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
    def __init__(self, text="", status=200, release="2026_03", next_url=None):
        self.text = text
        self.status_code = status
        self.headers = {"x-uniprot-release": release} if release else {}
        self.links = {"next": {"url": next_url}} if next_url else {}

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


def test_fetch_uses_the_search_endpoint_and_follows_next_pages():
    """One batch spread over pages: the `next` link is followed with no params of ours
    (it carries its own cursor), repeated headers are dropped, and every row is kept."""
    pages = [
        (FakeSession.HEADER + "\n" + ROWS["P25815"] + "\n", "https://rest.uniprot.org/next1"),
        (FakeSession.HEADER + "\n" + ROWS["Q9NPD8"] + "\n", None),
    ]

    class Paged(FakeSession):
        def get(self, url, params=None, timeout=None):
            assert timeout, "every request must carry a timeout"
            self.calls.append((url, params))
            text, nxt = pages[len(self.calls) - 1]
            return FakeResponse(text, release=self.release, next_url=nxt)

    s = Paged(ROWS)
    df, release = fetch_uniprot(["P25815", "Q9NPD8"], s, batch=100)
    assert [c[0] for c in s.calls] == [uniprot.URL, "https://rest.uniprot.org/next1"]
    assert uniprot.URL.endswith("/uniprotkb/search")     # not /stream: it 500s on OR-queries
    assert s.calls[0][1]["size"] == uniprot.PAGE_SIZE
    assert s.calls[1][1] is None
    assert df["Entry"].tolist() == ["P25815", "Q9NPD8"] and release == "2026_03"


def test_fetch_stops_if_a_next_link_repeats():
    """A server that keeps pointing at the same page must not loop forever."""
    class Looping(FakeSession):
        def get(self, url, params=None, timeout=None):
            self.calls.append(url)
            return FakeResponse(self.HEADER + "\n" + ROWS["P25815"] + "\n",
                                release=self.release, next_url=uniprot.URL)

    s = Looping(ROWS)
    df, _ = fetch_uniprot(["P25815"], s, batch=100)
    assert s.calls == [uniprot.URL] and df["Entry"].tolist() == ["P25815"]


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
    assert m["counts"] == {"requested": 6, "kept": 2, "not_found": 1, "inactive": 0,
                           "unreviewed": 1, "non_human": 1, "renamed": 1, "failed": 0}


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


# ------------------------------------------------- one protein per workbook row
ROW_TABLE = {
    "P0DP23": "P0DP23\tCALM1_HUMAN\tCALM1\treviewed\t9606\t",
    "P0DP24": "P0DP24\tCALM2_HUMAN\tCALM2\treviewed\t9606\t",
    "P0DP25": "P0DP25\tCALM3_HUMAN\tCALM3\treviewed\t9606\t",
    "Q96HY3": "Q96HY3\tX_HUMAN\tX\tunreviewed\t9606\t",
    "P11111": "P11111\tA_HUMAN\tAAA\treviewed\t9606\t",
}
CALM_IDS = ["Q96HY3", "P0DP23", "P0DP24", "P0DP25"]


def _with_rows(tmp_path, monkeypatch, ids, rows_text, table=ROW_TABLE):
    run, layout = _setup_run(tmp_path, monkeypatch, ids, FakeSession(table))
    if rows_text is not None:
        layout.rows_tsv.write_text("row\tgene\tids\tsolo\n" + rows_text)
    return run, layout


def _go(run):
    return cli.main(["uniprot", "--run-dir", str(run), "--no-notify"])


def test_rows_keep_gene_matching_accession(tmp_path, monkeypatch):
    run, layout = _with_rows(tmp_path, monkeypatch, CALM_IDS + ["P11111"],
                             "3\tcalm3 \tQ96HY3;P0DP23;P0DP24;P0DP25\t\n")
    assert _go(run) == 0
    assert pd.read_csv(layout.uniprot_csv, dtype=str)["Entry"].tolist() == ["P0DP25", "P11111"]
    d = layout.stage_dir("01_uniprot")
    ex = (d / "excluded.tsv").read_text()
    assert "P0DP23\tsame workbook row as P0DP25 (gene calm3)\n" in ex
    assert "P0DP24\tsame workbook row as P0DP25 (gene calm3)\n" in ex
    assert "Q96HY3\tunreviewed\n" in ex
    m = json.loads((d / "manifest.json").read_text())
    assert m["counts"]["row_collapsed"] == 2 and m["counts"]["kept"] == 2
    assert [Path(i["path"]).name for i in m["inputs"]] == ["ids.txt", "rows.tsv"]
    assert m["inputs"][1]["sha256"] == sha256_file(layout.rows_tsv)


def test_rows_ambiguous_or_unmatched_gene_stops_with_candidates(tmp_path, monkeypatch):
    run, layout = _with_rows(tmp_path, monkeypatch, CALM_IDS,
                             "7\tNOPE\tQ96HY3;P0DP23;P0DP24;P0DP25\t\n")
    with pytest.raises(ValueError, match=r"row 7.*NOPE.*P0DP23 \(CALM1\).*P0DP25 \(CALM3\)"):
        _go(run)
    assert not (layout.stage_dir("01_uniprot") / "manifest.json").exists()
    assert not layout.uniprot_csv.exists()


def test_rows_two_gene_matches_is_ambiguous(tmp_path, monkeypatch):
    table = dict(ROW_TABLE, P0DP24="P0DP24\tCALM2_HUMAN\tCALM3\treviewed\t9606\t")
    run, _ = _with_rows(tmp_path, monkeypatch, CALM_IDS,
                        "7\tCALM3\tP0DP23;P0DP24;P0DP25\t\n", table)
    with pytest.raises(ValueError, match=r"row 7"):
        _go(run)


def test_rows_accession_kept_via_another_row_is_not_dropped(tmp_path, monkeypatch):
    # P0DP23 also is the sole accession of another row (solo) -> stays alongside CALM3
    run, layout = _with_rows(tmp_path, monkeypatch, ["P0DP23", "P0DP25"],
                             "3\tCALM3\tP0DP23;P0DP25\tP0DP23\n")
    assert _go(run) == 0
    assert pd.read_csv(layout.uniprot_csv, dtype=str)["Entry"].tolist() == ["P0DP23", "P0DP25"]
    m = json.loads((layout.stage_dir("01_uniprot") / "manifest.json").read_text())
    assert m["counts"]["row_collapsed"] == 0


def test_rows_accession_winning_another_row_is_not_dropped(tmp_path, monkeypatch):
    run, layout = _with_rows(tmp_path, monkeypatch, ["P0DP23", "P0DP24", "P0DP25"],
                             "3\tCALM3\tP0DP23;P0DP25\t\n9\tCALM1\tP0DP23;P0DP24\t\n")
    assert _go(run) == 0
    assert pd.read_csv(layout.uniprot_csv, dtype=str)["Entry"].tolist() == ["P0DP23", "P0DP25"]


def test_rows_zero_or_one_kept_candidate_is_untouched(tmp_path, monkeypatch):
    run, layout = _with_rows(tmp_path, monkeypatch, ["Q96HY3", "P11111"],
                             "2\tWHATEVER\tQ96HY3;P11111\t\n5\tGONE\tQ96HY3\t\n")
    assert _go(run) == 0
    assert pd.read_csv(layout.uniprot_csv, dtype=str)["Entry"].tolist() == ["P11111"]


def test_no_rows_tsv_means_no_row_collapsed_count_and_all_kept(tmp_path, monkeypatch):
    run, layout = _with_rows(tmp_path, monkeypatch, CALM_IDS, None)
    assert _go(run) == 0
    assert pd.read_csv(layout.uniprot_csv, dtype=str)["Entry"].tolist() == [
        "P0DP23", "P0DP24", "P0DP25"]
    m = json.loads((layout.stage_dir("01_uniprot") / "manifest.json").read_text())
    assert "row_collapsed" not in m["counts"]


def test_rows_tsv_change_invalidates_resume(tmp_path, monkeypatch):
    run, layout = _with_rows(tmp_path, monkeypatch, ["P0DP23", "P0DP25"],
                             "3\tCALM3\tP0DP23;P0DP25\t\n")
    assert _go(run) == 0
    assert uniprot.units(layout, None).count == 0
    layout.rows_tsv.write_text("row\tgene\tids\tsolo\n3\tCALM1\tP0DP23;P0DP25\t\n")
    assert uniprot.units(layout, None).count == 1
    assert _go(run) == 0
    assert pd.read_csv(layout.uniprot_csv, dtype=str)["Entry"].tolist() == ["P0DP23"]
    layout.rows_tsv.unlink()                               # removal also invalidates
    assert uniprot.units(layout, None).count == 1


# ------------------------------------------------------- AD regression (2023 lists)
AD = Path(os.environ.get("ADEP_DATA", "")) / "disease_subtypes" / "Alzhemier"


@pytest.fixture(scope="module")
def ad_workbook(tmp_path_factory):
    import shutil
    local = tmp_path_factory.mktemp("ad") / "w.xlsx"
    shutil.copy(AD / "PPIN_ROSMAP_clusterFeatures_v7.xlsx", local)
    return local


def _ad_2023(n):
    if n == 1:
        df = pd.read_csv(AD / "Cluster1" / "Cluster10_UniProt.csv", dtype=str)
        return set(df["Entry"])
    return {l.strip() for l in (AD / f"Cluster{n}" / f"Cluster{n}.txt").read_text().splitlines()
            if l.strip()}


@pytest.mark.network
@pytest.mark.regression
@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_ad_cluster_kept_matches_2023(n, ad_workbook, tmp_path, monkeypatch):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    run = tmp_path / "run"
    (run / "inputs").mkdir(parents=True)
    import shutil
    shutil.copy(ad_workbook, run / "inputs" / "w.xlsx")
    assert cli.main(["ids", "--run-dir", str(run), "--no-notify", "--xlsx", "inputs/w.xlsx",
                     "--sheet", f"cluster{n}", "--id-column", "UNIPROTID",
                     "--gene-column", "GeneID"]) == 0
    assert cli.main(["uniprot", "--run-dir", str(run), "--no-notify"]) == 0
    layout = RunLayout(run)
    got = set(pd.read_csv(layout.uniprot_csv, dtype=str)["Entry"])
    ref = _ad_2023(n)
    d = layout.stage_dir("01_uniprot")
    print(f"cluster{n}: got {len(got)} ref {len(ref)} only_now={sorted(got - ref)} "
          f"only_2023={sorted(ref - got)}")
    print((d / "renamed.tsv").read_text())
    # Cluster3: the 2023 list holds 'P8529', a truncated accession (not a UniProt entry);
    # the v7 workbook has the full P85299. This is the one explained difference.
    known = {3: ({"P85299"}, {"P8529"})}.get(n, (set(), set()))
    assert (got - ref, ref - got) == known


# --- inactive (deleted or demerged) entries ---------------------------------

INACTIVE_TSV = (
    "Entry\tEntry Name\tGene Names (primary)\tReviewed\tOrganism (ID)\tPDB\n"
    "P25815\tS100P_HUMAN\tS100P\treviewed\t9606\t1J55;\n"
    "A0A024QYU9\tA0A024QYU9_HUMAN\t\t\t\t\n"
)


def test_parse_keeps_inactive_stub_row_as_na_organism():
    """UniProt answers for a deleted entry with accession and entry name only. int64 on
    the empty organism raised ValueError and took the whole stage down."""
    df = parse_uniprot_tsv(INACTIVE_TSV)
    assert df["Entry"].tolist() == ["P25815", "A0A024QYU9"]
    assert df["Organism (ID)"].tolist()[0] == 9606
    assert df["Organism (ID)"].isna().tolist() == [False, True]


def test_classify_reports_inactive_entries_and_does_not_requery_them():
    df = parse_uniprot_tsv(INACTIVE_TSV)
    calls = []

    def requery(acc):
        calls.append(acc)
        return None

    res = classify(["P25815", "A0A024QYU9", "NOPE01"], df, requery=requery)
    assert res["kept"]["Entry"].tolist() == ["P25815"]
    assert res["inactive"] == ["A0A024QYU9"]
    assert res["not_found"] == ["NOPE01"]
    assert res["unreviewed"] == [] and res["non_human"] == []
    assert calls == ["NOPE01"]      # the inactive entry was returned, so never re-queried
