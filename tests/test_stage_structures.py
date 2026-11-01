import gzip
import json
import os
import shutil
from pathlib import Path

import pandas as pd
import pytest
import requests

from polyphar import cli
from polyphar.common import RunLayout
from polyphar.stages import structures
from polyphar.stages.structures import (ChainRejected, extract_chain, fetch_alphafold,
                                        load_sifts, select_chains, verify_chain)

FIX = Path(__file__).parent / "fixtures"
SIFTS = FIX / "sifts_subset.tsv"
AF_JSON = FIX / "af_api_P25815.json"
AF_PDB = ("HEADER                                            01-AUG-25\n"
          "ATOM      1  N   MET A   1      -1.000   2.000   3.000  1.00 90.00           N\n"
          "END\n")
LAST_MODIFIED = "Sun, 25 Oct 2026 02:30:33 GMT"


def _pairs(*pairs):
    return pd.DataFrame(pairs, columns=["UniProt ID", "PDB"])


def _cif_2mjw(tmp_path):
    out = tmp_path / "2mjw.cif"
    with gzip.open(FIX / "2mjw.cif.gz", "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def _strip_loops(text, *prefixes):
    """Drop the mmCIF loops whose first tag starts with one of `prefixes`."""
    out, skipping, lines = [], False, text.splitlines()
    for i, line in enumerate(lines):
        if skipping:
            skipping = line.strip() != "#"
            if not skipping:
                out.append(line)
            continue
        if line.strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].startswith(prefixes):
            skipping = True
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def _cif_hetatm_only(text, chain="B"):
    """2MJW whose author chain `chain` is made of HETATM records only."""
    lines = []
    for line in text.splitlines():
        f = line.split()
        if f and f[0] == "ATOM" and f[-3] == chain:
            line = "HETATM" + line[4:]
        lines.append(line)
    return "\n".join(lines) + "\n"


# --- fake HTTP ----------------------------------------------------------------

class FakeResponse:
    def __init__(self, content=b"", status=200, headers=None):
        self.content = content if isinstance(content, bytes) else content.encode()
        self.status_code = status
        self.headers = headers or {}

    @property
    def text(self):
        return self.content.decode("utf-8", errors="replace")

    def json(self):
        return json.loads(self.text)

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(f"{self.status_code} error", response=self)


class FakeSession:
    """Routes SIFTS, RCSB and AFDB URLs to canned responses; records every GET."""

    def __init__(self, cifs=None, af=None, af_pdb=AF_PDB, sifts_status=200, sifts_text=None):
        self.cifs = cifs or {}            # PDB ID (upper) -> cif text
        self.af = af or {}                # accession -> API json text
        self.af_pdb = af_pdb
        self.sifts_status = sifts_status
        self.sifts_text = SIFTS.read_text() if sifts_text is None else sifts_text
        self.calls = []

    def get(self, url, timeout=None, **kw):
        assert timeout is not None, "every request needs a timeout"
        self.calls.append(url)
        if url == structures.SIFTS_URL:
            if self.sifts_status != 200:
                return FakeResponse(b"", self.sifts_status)
            return FakeResponse(gzip.compress(self.sifts_text.encode()),
                                headers={"Last-Modified": LAST_MODIFIED})
        if url.startswith("https://files.rcsb.org/download/"):
            pdb = url.rsplit("/", 1)[1][:-4]
            if pdb in self.cifs:
                return FakeResponse(self.cifs[pdb])
            return FakeResponse(b"not found", 404)
        if url.startswith(structures.AF_API):
            acc = url.rsplit("/", 1)[1]
            if acc in self.af:
                return FakeResponse(self.af[acc])
            return FakeResponse(b'{"detail":"not found"}', 404)
        if url.startswith("https://alphafold.ebi.ac.uk/files/"):
            return FakeResponse(self.af_pdb)
        raise AssertionError(f"unexpected URL {url}")

    def downloads(self):
        return [u for u in self.calls]


def _af_entries(*entries):
    base = json.loads(AF_JSON.read_text())[0]
    out = []
    for acc, entry_id, version in entries:
        e = dict(base, uniprotAccession=acc, entryId=entry_id, modelEntityId=entry_id,
                 latestVersion=version,
                 pdbUrl=f"https://alphafold.ebi.ac.uk/files/{entry_id}-model_v{version}.pdb")
        out.append(e)
    return json.dumps(out)


# --- chain selection ------------------------------------------------------------

def test_load_sifts_reads_comment_header_and_keeps_strings():
    s = load_sifts(SIFTS)
    assert list(s.columns[:3]) == ["PDB", "CHAIN", "SP_PRIMARY"]
    assert set(s["PDB"]) >= {"2mjw", "7kzr", "5cvd", "2hqw", "1q4k"}
    assert all(isinstance(v, str) for v in s["PDB_BEG"])   # empty stays "", not NaN
    assert "" in set(s["PDB_BEG"])


def test_selects_target_chain_not_first():
    sel, unmapped = select_chains(
        _pairs(("P25815", "2MJW"), ("Q9NPD8", "7KZR"), ("P49450", "5CVD")), load_sifts(SIFTS))
    got = {(r["UniProt ID"], r["PDB"]): r["Chain ID"] for _, r in sel.iterrows()}
    assert got == {("P25815", "2MJW"): "B", ("Q9NPD8", "7KZR"): "X", ("P49450", "5CVD"): "D"}
    assert unmapped.empty
    assert sel.set_index("PDB").loc["2MJW", "SP coverage"] == 94


def test_ortholog_entry_unmapped():
    sel, unmapped = select_chains(_pairs(("Q05586", "2HQW")), load_sifts(SIFTS))
    assert sel.empty
    assert unmapped.to_dict("records") == [
        {"UniProt ID": "Q05586", "PDB": "2HQW", "reason": "no SIFTS chain for accession"}]


def test_absent_pdb_unmapped():
    sel, unmapped = select_chains(_pairs(("P25815", "1ZZZ")), load_sifts(SIFTS))
    assert sel.empty
    assert unmapped.to_dict("records") == [
        {"UniProt ID": "P25815", "PDB": "1ZZZ", "reason": "PDB absent from SIFTS"}]


CHIMERA_SIFTS = ("# test\nPDB\tCHAIN\tSP_PRIMARY\tRES_BEG\tRES_END\tPDB_BEG\tPDB_END\t"
                 "SP_BEG\tSP_END\n"
                 "9xyz\tA\tP11111\t1\t100\t1\t100\t1\t100\n"
                 "9xyz\tA\tP22222\t106\t205\t105\t204\t1\t100\n"
                 "9xyz\tA\tP11111\t210\t219\t209\tNone\t101\t110\n"
                 "9xyz\tB\tP11111\t1\t50\t1\t50\t1\t50\n")


def test_chimeric_chain_records_segments_and_partners(tmp_path):
    """A chain SIFTS maps to several accessions (fusion construct) keeps this accession's
    author-numbered segments and names the partner accessions."""
    path = tmp_path / "sifts.tsv"
    path.write_text(CHIMERA_SIFTS)
    sel, _ = select_chains(_pairs(("P11111", "9XYZ"), ("P22222", "9XYZ")), load_sifts(path))
    rows = {r["UniProt ID"]: r for r in sel.to_dict("records")}
    assert rows["P11111"]["Chain ID"] == "A"          # coverage 110 > 50
    assert rows["P11111"]["SIFTS segments"] == "1:100;209:None"
    assert rows["P11111"]["SIFTS other accessions"] == "P22222"
    assert rows["P22222"]["SIFTS segments"] == "105:204"
    assert rows["P22222"]["SIFTS other accessions"] == "P11111"


def test_single_protein_chain_has_no_partners():
    sel, _ = select_chains(_pairs(("P25815", "2MJW"),), load_sifts(SIFTS))
    assert sel.loc[0, "SIFTS segments"] == "1:94"
    assert sel.loc[0, "SIFTS other accessions"] == ""


def test_tie_breaks_alphabetically():
    sel, _ = select_chains(_pairs(("P99999", "9ZZA"), ("P53350", "1Q4K")), load_sifts(SIFTS))
    assert sel["Chain ID"].tolist() == ["A", "A"]


# --- extraction -------------------------------------------------------------------

def _atoms(pdb_path):
    return [l for l in Path(pdb_path).read_text().splitlines()
            if l.startswith(("ATOM", "HETATM"))]


def test_extract_chain_writes_only_requested_chain(tmp_path):
    cif = _cif_2mjw(tmp_path)
    sifts = load_sifts(SIFTS).set_index(["PDB", "CHAIN"])
    for chain in ("B", "D"):      # D: author chain D is label_asym C (P25815, not RAGE)
        out = tmp_path / f"2mjw_{chain}.pdb"
        info = extract_chain(cif, chain, out)
        assert info == {"Resolution": "-", "Structure method": "SOLUTION NMR"}
        atoms = _atoms(out)
        assert {l[21] for l in atoms} == {chain}
        resnums = sorted({int(l[22:26]) for l in atoms})
        beg, end = (int(sifts.loc[("2mjw", chain), c]) for c in ("PDB_BEG", "PDB_END"))
        assert resnums == list(range(beg, end + 1))
        text = out.read_text()
        assert text.count("\nMODEL") + text.startswith("MODEL") == 10   # all NMR models kept
    assert not list(tmp_path.glob("*.tmp"))


def test_extract_multichar_chain_written_as_A(tmp_path):
    lines = []
    for line in _cif_2mjw(tmp_path).read_text().splitlines():
        f = line.split()
        if f and f[0] in ("ATOM", "HETATM") and f[-3] == "B":
            f[-3] = "BB"
            line = " ".join(f)
        lines.append(line)
    cif = tmp_path / "multi.cif"
    cif.write_text("\n".join(lines) + "\n")
    out = tmp_path / "multi_BB.pdb"
    extract_chain(cif, "BB", out)
    atoms = _atoms(out)
    assert atoms and {l[21] for l in atoms} == {"A"}


def test_extract_hetatm_only_chain_rejected(tmp_path):
    cif = tmp_path / "het.cif"
    cif.write_text(_cif_hetatm_only(_cif_2mjw(tmp_path).read_text()))
    out = tmp_path / "het_B.pdb"
    with pytest.raises(ChainRejected, match="chain has no ATOM records"):
        extract_chain(cif, "B", out)
    assert not out.exists() and not list(tmp_path.glob("*.tmp"))


def test_cif_check_accepts_target_chain(tmp_path):
    cif = _cif_2mjw(tmp_path)
    assert verify_chain(cif, "B", "P25815") == ("verified", ["P25815"])
    assert verify_chain(cif, "D", "P25815") == ("verified", ["P25815"])   # author D (label C)


def test_cif_check_reports_mismatch(tmp_path):
    # author chain A of 2MJW is RAGE: the cif's own annotation disagrees with P25815
    assert verify_chain(_cif_2mjw(tmp_path), "A", "P25815") == ("mismatch", ["Q15109"])


def test_cif_check_without_annotation_is_unverified(tmp_path):
    cif = tmp_path / "bare.cif"
    cif.write_text(_strip_loops(_cif_2mjw(tmp_path).read_text(), "_struct_ref."))
    assert "_struct_ref." not in cif.read_text()
    assert verify_chain(cif, "B", "P25815") == ("unverified", [])


def test_extract_missing_chain_raises_and_writes_nothing(tmp_path):
    out = tmp_path / "2mjw_Z.pdb"
    with pytest.raises(ValueError, match="chain Z"):
        extract_chain(_cif_2mjw(tmp_path), "Z", out)
    assert not out.exists()


# --- AlphaFold --------------------------------------------------------------------

def test_af_missing_model_writes_nothing(tmp_path):
    s = FakeSession()
    assert fetch_alphafold("Q8WZ42", tmp_path, s) == (None, None)
    assert list(tmp_path.iterdir()) == []
    s = FakeSession(af={"P11111": _af_entries(("P11111", "AF-P11111-2-F1", 6))})
    assert fetch_alphafold("P11111", tmp_path, s) == (None, None)    # no F1 entry
    assert list(tmp_path.iterdir()) == []


def test_af_error_body_rejected(tmp_path):
    s = FakeSession(af={"P25815": AF_JSON.read_text()}, af_pdb="<html>Forbidden</html>")
    with pytest.raises(ValueError, match="PDB"):
        fetch_alphafold("P25815", tmp_path, s)
    assert list(tmp_path.iterdir()) == []


def test_af_picks_F1(tmp_path):
    api = _af_entries(("Q05586", "AF-Q05586-2-F1", 6), ("Q05586", "AF-Q05586-F1", 5))
    s = FakeSession(af={"Q05586": api})
    path, version = fetch_alphafold("Q05586", tmp_path, s)
    assert path == tmp_path / "AF-Q05586-F1.pdb"
    assert version == 5
    assert s.calls[-1] == "https://alphafold.ebi.ac.uk/files/AF-Q05586-F1-model_v5.pdb"
    assert path.read_text() == AF_PDB


def test_af_real_api_fixture(tmp_path):
    s = FakeSession(af={"P25815": AF_JSON.read_text()})
    path, version = fetch_alphafold("P25815", tmp_path, s)
    assert (path.name, version) == ("AF-P25815-F1.pdb", 6)


# --- stage --------------------------------------------------------------------------

def _setup(tmp_path, monkeypatch, rows, session):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    monkeypatch.setattr(structures, "http_session", lambda *a, **k: session)
    run = tmp_path / "run"
    layout = RunLayout(run)
    layout.uniprot_csv.parent.mkdir(parents=True)
    df = pd.DataFrame(rows, columns=["Entry", "Entry Name", "PDB"])
    df.insert(2, "Gene Names (primary)", "")
    df.insert(3, "Reviewed", "reviewed")
    df.insert(4, "Organism (ID)", "9606")
    df.to_csv(layout.uniprot_csv, index=False)
    return run, layout


ROWS = [("P25815", "S100P_HUMAN", "2MJW;1ZZZ"), ("Q05586", "NMDZ1_HUMAN", "2HQW"),
        ("O00000", "NOPDB_HUMAN", "")]


def _session(tmp_path):
    return FakeSession(cifs={"2MJW": _cif_2mjw(tmp_path).read_text()},
                       af={"P25815": AF_JSON.read_text()})


def _main(run):
    return cli.main(["structures", "--run-dir", str(run), "--no-notify", "--nprocs", "2"])


def test_units_counts_pairs_plus_accessions(tmp_path, monkeypatch):
    run, layout = _setup(tmp_path, monkeypatch, ROWS, _session(tmp_path))
    wu = structures.units(layout, None)
    assert (wu.count, wu.cost_basis) == (3 + 3, "download")


def test_run_writes_outputs(tmp_path, monkeypatch):
    s = _session(tmp_path)
    run, layout = _setup(tmp_path, monkeypatch, ROWS, s)
    assert _main(run) == 0
    d = layout.structures
    info = pd.read_csv(layout.protein_info, dtype=str, keep_default_na=False)
    assert list(info.columns) == ["UniProt ID", "Gene Name", "PDB", "Chain ID", "AlphaFold",
                                  "Resolution", "Structure method", "SIFTS coverage", "File",
                                  "SIFTS segments", "SIFTS other accessions",
                                  "SIFTS segments source"]
    assert info.to_dict("records") == [
        {"UniProt ID": "P25815", "Gene Name": "S100P_HUMAN", "PDB": "2MJW", "Chain ID": "B",
         "AlphaFold": "False", "Resolution": "-", "Structure method": "SOLUTION NMR",
         "SIFTS coverage": "94", "File": "02_structures/P25815/2mjw_B.pdb",
         "SIFTS segments": "1:94", "SIFTS other accessions": "",
         "SIFTS segments source": "pdb"},
        {"UniProt ID": "P25815", "Gene Name": "S100P_HUMAN", "PDB": "-", "Chain ID": "A",
         "AlphaFold": "True", "Resolution": "-", "Structure method": "-",
         "SIFTS coverage": "-", "File": "02_structures/P25815/AF-P25815-F1.pdb",
         "SIFTS segments": "", "SIFTS other accessions": "", "SIFTS segments source": ""},
    ]
    assert (d / "P25815" / "2mjw.cif").read_text().startswith("data_")
    assert not (d / "P25815" / "2mjw_A.pdb").exists()
    assert not list(d.rglob("*.tmp"))
    unmapped = (d / "unmapped.tsv").read_text()
    assert unmapped == ("UniProt ID\tPDB\treason\n"
                        "P25815\t1ZZZ\tPDB absent from SIFTS\n"
                        "Q05586\t2HQW\tno SIFTS chain for accession\n")
    failed = (d / "failed.tsv").read_text()
    assert "P25815/1ZZZ\tstructures\tPDB absent from SIFTS" in failed
    assert "Q05586/2HQW\tstructures\tno SIFTS chain for accession" in failed
    assert (d / "sifts" / "pdb_chain_uniprot.tsv.gz").is_file()
    m = json.loads((d / "manifest.json").read_text())
    assert m["versions"]["sifts_last_modified"] == LAST_MODIFIED
    assert m["versions"]["alphafold_versions"] == [6]
    assert m["counts"] == {"accessions": 3, "pdb_entries_listed": 3, "chains_extracted": 1,
                           "unmapped": {"PDB absent from SIFTS": 1,
                                        "no SIFTS chain for accession": 1},
                           "chimeric_chains": 0, "segments_from_seqres": 0,
                           "segments_unverified": 0,
                           "af_found": 1, "af_missing": 2, "cif_unverified": 0, "annotation_mismatch": 0,
                           "failed": 0}


def test_resume_skips_existing(tmp_path, monkeypatch):
    s = _session(tmp_path)
    run, layout = _setup(tmp_path, monkeypatch, ROWS, s)
    assert _main(run) == 0
    first = (layout.protein_info).read_text()
    assert structures.units(layout, None).count == 0
    n = len(s.calls)
    assert _main(run) == 0
    assert len(s.calls) == n                       # second run: 0 downloads
    assert layout.protein_info.read_text() == first
    # a lost extraction is redone from the cif already on disk, still no download
    (layout.structures / "P25815" / "2mjw_B.pdb").unlink()
    assert structures.units(layout, None).count == 1
    assert _main(run) == 0
    assert len(s.calls) == n
    # a lost cif + extraction and a half-written AF model are fetched again
    (layout.structures / "P25815" / "2mjw_B.pdb").unlink()
    (layout.structures / "P25815" / "2mjw.cif").unlink()
    (layout.structures / "P25815" / "AF-P25815-F1.pdb").write_text("")
    assert structures.units(layout, None).count == 2
    assert _main(run) == 0
    assert len([u for u in s.calls[n:] if "rcsb" in u]) == 1
    assert layout.protein_info.read_text() == first


def test_entry_download_failure_is_recorded_and_stage_continues(tmp_path, monkeypatch):
    s = _session(tmp_path)
    rows = [("P25815", "S100P_HUMAN", "2MJW"), ("P99999", "TEST_HUMAN", "9ZZA")]
    run, layout = _setup(tmp_path, monkeypatch, rows, s)
    assert _main(run) == 3        # finished with a failed unit (item 13)
    failed = (layout.structures / "failed.tsv").read_text()
    assert "P99999/9ZZA\tstructures\t" in failed and "404" in failed
    m = json.loads((layout.structures / "manifest.json").read_text())
    assert m["counts"]["failed"] == 1
    assert m["counts"]["chains_extracted"] == 1
    assert structures.units(layout, None).count == 1     # retried on the next run


def test_annotation_mismatch_warns_and_extracts(tmp_path, monkeypatch):
    # SIFTS (the authority) says 2MJW chain A is P25815 with the larger coverage;
    # the cif's depositor annotation says A is RAGE (Q15109): extract, but list it.
    sifts = SIFTS.read_text().replace("2mjw\tA\tQ15109\t3\t101\t23\t121\t23\t121",
                                      "2mjw\tA\tP25815\t3\t101\t23\t121\t1\t99")
    assert "2mjw\tA\tP25815" in sifts
    s = FakeSession(cifs={"2MJW": _cif_2mjw(tmp_path).read_text()}, sifts_text=sifts)
    run, layout = _setup(tmp_path, monkeypatch, [("P25815", "S100P_HUMAN", "2MJW")], s)
    assert _main(run) == 0
    d = layout.structures
    assert (d / "P25815" / "2mjw_A.pdb").is_file()
    expected = "UniProt ID\tPDB\tChain ID\tcif_accessions\nP25815\t2MJW\tA\tQ15109\n"
    assert (d / "annotation_mismatch.tsv").read_text() == expected
    assert (d / "failed.tsv").read_text() == "item\tstage\treason\n"
    info = pd.read_csv(layout.protein_info, dtype=str, keep_default_na=False)
    assert info.loc[info["PDB"] == "2MJW", "Chain ID"].tolist() == ["A"]
    counts = json.loads((d / "manifest.json").read_text())["counts"]
    assert (counts["chains_extracted"], counts["annotation_mismatch"], counts["failed"]) == (1, 1, 0)
    assert _main(run) == 0                                     # rerun: same, from disk
    assert (d / "annotation_mismatch.tsv").read_text() == expected
    assert json.loads((d / "manifest.json").read_text())["counts"] == counts


def test_rejected_chain_recorded_and_not_retried(tmp_path, monkeypatch):
    reason = "chain has no ATOM records"
    s = FakeSession(cifs={"2MJW": _cif_hetatm_only(_cif_2mjw(tmp_path).read_text())})
    run, layout = _setup(tmp_path, monkeypatch, [("P25815", "S100P_HUMAN", "2MJW")], s)
    assert _main(run) == 3        # finished with a failed unit (item 13)
    d = layout.structures
    assert not (d / "P25815" / "2mjw_B.pdb").exists()
    assert not list(d.rglob("*.tmp"))
    assert f"P25815/2MJW\tstructures\tChainRejected: {reason}" in (d / "failed.tsv").read_text()
    info = pd.read_csv(layout.protein_info, dtype=str, keep_default_na=False)
    assert "2MJW" not in set(info["PDB"])
    counts = json.loads((d / "manifest.json").read_text())["counts"]
    assert (counts["chains_extracted"], counts["failed"]) == (0, 1)
    assert structures.units(layout, None).count == 0          # not pending forever
    n = len(s.calls)
    assert _main(run) == 3        # finished with a failed unit (item 13)
    assert len(s.calls) == n
    assert reason in (d / "failed.tsv").read_text()           # still visible on rerun
    assert json.loads((d / "manifest.json").read_text())["counts"] == counts


def test_cif_without_annotation_is_extracted_and_counted(tmp_path, monkeypatch):
    s = FakeSession(cifs={"2MJW": _strip_loops(_cif_2mjw(tmp_path).read_text(), "_struct_ref.")})
    run, layout = _setup(tmp_path, monkeypatch, [("P25815", "S100P_HUMAN", "2MJW")], s)
    assert _main(run) == 0
    assert (layout.structures / "P25815" / "2mjw_B.pdb").is_file()
    m = json.loads((layout.structures / "manifest.json").read_text())
    assert (m["counts"]["chains_extracted"], m["counts"]["cif_unverified"]) == (1, 1)
    assert m["counts"]["failed"] == 0
    assert _main(run) == 0                                     # rerun: same count from disk
    m = json.loads((layout.structures / "manifest.json").read_text())
    assert m["counts"]["cif_unverified"] == 1


def test_sifts_failure_aborts(tmp_path, monkeypatch):
    s = FakeSession(sifts_status=404)
    run, layout = _setup(tmp_path, monkeypatch, ROWS, s)
    with pytest.raises(RuntimeError, match="SIFTS"):
        _main(run)
    assert not (layout.structures / "manifest.json").exists()
    assert not (layout.structures / "sifts" / "pdb_chain_uniprot.tsv.gz").exists()


def test_uniprot_csv_is_not_modified(tmp_path, monkeypatch):
    s = _session(tmp_path)
    run, layout = _setup(tmp_path, monkeypatch, ROWS, s)
    before = layout.uniprot_csv.read_bytes()
    assert _main(run) == 0
    assert layout.uniprot_csv.read_bytes() == before


# --- regression / network --------------------------------------------------------------

@pytest.mark.regression
def test_extraction_matches_2023_where_chain_was_right(tmp_path):
    import prody
    legacy_dir = (Path(os.environ["ADEP_DATA"]) / "disease_subtypes" / "Liver_cancer" / "net2"
                  / "structures" / "P53350")
    out = tmp_path / "1q4k_B.pdb"
    extract_chain(legacy_dir / "1q4k.cif", "B", out)
    prody.LOGGER.verbosity = "none"
    new, old = prody.parsePDB(str(out)), prody.parsePDB(str(legacy_dir / "1q4k_B.pdb"))
    assert new.numAtoms() == old.numAtoms()
    for get in ("getNames", "getResnames", "getResnums", "getChids"):
        assert list(getattr(new, get)()) == list(getattr(old, get)())
    assert abs(new.getCoords() - old.getCoords()).max() == 0.0


@pytest.mark.network
def test_end_to_end_two_proteins(tmp_path, monkeypatch):
    from polyphar.common import http_session
    run, layout = _setup(tmp_path, monkeypatch,
                         [("P25815", "S100P_HUMAN", "2MJW"), ("Q9NPD8", "UBE2T_HUMAN", "7KZR")],
                         http_session())
    assert _main(run) == 0
    info = pd.read_csv(layout.protein_info, dtype=str, keep_default_na=False)
    pdb_rows = {(r["PDB"], r["Chain ID"]) for _, r in info.iterrows() if r["AlphaFold"] == "False"}
    assert pdb_rows == {("2MJW", "B"), ("7KZR", "X")}
    assert ("2MJW", "A") not in pdb_rows
    x = layout.structures / "Q9NPD8" / "7kzr_X.pdb"
    assert {l[21] for l in _atoms(x)} == {"X"}
    assert len(info[info["AlphaFold"] == "True"]) == 2


# --- SIFTS segment bounds from SEQRES when PDB_BEG/PDB_END are unusable ------------------

_SCHEME_TAGS = ("asym_id entity_id seq_id mon_id ndb_seq_num pdb_seq_num auth_seq_num "
                "pdb_mon_id auth_mon_id pdb_strand_id pdb_ins_code hetero").split()


def scheme_cif(path, chains):
    """A minimal mmCIF with only _pdbx_poly_seq_scheme. chains: author chain ->
    [(seq_id, author number or None if unobserved, insertion code)]."""
    rows = []
    for k, (chain, residues) in enumerate(chains.items()):
        for seq, num, icode in residues:
            auth = "?" if num is None else str(num)
            pdb_seq = str(num if num is not None else seq)
            mon = "?" if num is None else "ALA"
            rows.append(f"{chr(65 + k)} 1 {seq} ALA {seq} {pdb_seq} {auth} {mon} {mon} "
                        f"{chain} {icode or '.'} n")
    path.write_text("data_9XYZ\n#\nloop_\n" + "".join(f"_pdbx_poly_seq_scheme.{t}\n"
                                                      for t in _SCHEME_TAGS)
                    + "\n".join(rows) + "\n#\n")
    return path


def _scheme(tmp_path):
    # chain A: seq 1..219; author = seq for 1..105, seq-1 after (a deletion in the
    # numbering), seq 215..219 unobserved; seq 53 is author 52A (insertion code)
    res = []
    for s in range(1, 220):
        num = s if s <= 105 else s - 1
        icode = ""
        if s == 53:
            num, icode = 52, "A"
        res.append((s, None if s >= 215 else num, icode))
    return structures._cif_dict(scheme_cif(tmp_path / "9xyz.cif", {"A": res}))


def test_usable_pdb_bounds_are_kept(tmp_path):
    cd = _scheme(tmp_path)
    assert structures.author_segments("1:100;105:204", "1:100;106:205", cd, "A") == \
        ("1:100;105:204", "pdb")


def test_empty_pdb_bounds_come_from_seqres(tmp_path):
    cd = _scheme(tmp_path)
    # seq 210..219 -> author 209..213 (214..218 unobserved); segment 1 keeps its PDB bounds
    assert structures.author_segments("1:100;209:None", "1:100;210:219", cd, "A") == \
        ("1:100;209:213", "seqres")
    assert structures.author_segments(":", "50:60", cd, "A") == ("50:60", "seqres")
    assert structures.author_segments("None:None", "53:60", cd, "A") == ("52A:60", "seqres")


def test_unobserved_or_unmappable_segments(tmp_path):
    cd = _scheme(tmp_path)
    # a wholly unobserved segment is dropped; the rest stays
    assert structures.author_segments("1:10;:", "1:10;215:219", cd, "A") == ("1:10", "seqres")
    # nothing observed, unknown chain, or unusable SEQRES bounds: left as SIFTS gave it
    assert structures.author_segments(":", "215:219", cd, "A") == (":", "")
    assert structures.author_segments(":", "1:10", cd, "Z") == (":", "")
    assert structures.author_segments(":", "None:10", cd, "A") == (":", "")
    assert structures.author_segments("", "", cd, "A") == ("", "")


def test_run_maps_empty_pdb_bounds_through_the_cif_seqres(tmp_path, monkeypatch):
    """Real 2MJW: author chain D numbers P25815 95..188 for SEQRES 1..94. With SIFTS
    PDB_BEG/PDB_END empty the segment is recovered from _pdbx_poly_seq_scheme."""
    text = SIFTS.read_text().splitlines()
    out = []
    for line in text:
        f = line.split("\t")
        if f[:2] == ["2mjw", "B"]:
            continue                      # chain D is then the one selected for P25815
        if f[:2] == ["2mjw", "D"]:
            f[5] = f[6] = ""
        out.append("\t".join(f))
    s = FakeSession(cifs={"2MJW": _cif_2mjw(tmp_path).read_text()},
                    af={"P25815": AF_JSON.read_text()}, sifts_text="\n".join(out) + "\n")
    run, layout = _setup(tmp_path, monkeypatch, [("P25815", "S100P_HUMAN", "2MJW")], s)
    assert _main(run) == 0
    info = pd.read_csv(layout.protein_info, dtype=str, keep_default_na=False)
    row = info[info["PDB"] == "2MJW"].iloc[0]
    assert (row["Chain ID"], row["SIFTS segments"], row["SIFTS segments source"]) == \
        ("D", "95:188", "seqres")
    assert info[info["AlphaFold"] == "True"].iloc[0]["SIFTS segments source"] == ""
    m = json.loads((layout.structures / "manifest.json").read_text())
    assert m["counts"]["segments_from_seqres"] == 1 and m["counts"]["segments_unverified"] == 0
