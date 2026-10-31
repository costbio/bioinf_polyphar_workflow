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


def _cif_relettered(text):
    """2MJW whose own annotation says author chain B is RAGE (Q15109), as after a
    remediation that re-lettered chains while the reused SIFTS still says B = S100P."""
    old = "\n2 2 2MJW B 1 ? 94  ? P25815"
    assert old in text
    return text.replace(old, "\n2 1 2MJW B 1 ? 94  ? Q15109")


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

    def __init__(self, cifs=None, af=None, af_pdb=AF_PDB, sifts_status=200):
        self.cifs = cifs or {}            # PDB ID (upper) -> cif text
        self.af = af or {}                # accession -> API json text
        self.af_pdb = af_pdb
        self.sifts_status = sifts_status
        self.calls = []

    def get(self, url, timeout=None, **kw):
        assert timeout is not None, "every request needs a timeout"
        self.calls.append(url)
        if url == structures.SIFTS_URL:
            if self.sifts_status != 200:
                return FakeResponse(b"", self.sifts_status)
            return FakeResponse(gzip.compress(SIFTS.read_bytes()),
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
    assert verify_chain(cif, "B", "P25815") == "verified"
    assert verify_chain(cif, "D", "P25815") == "verified"   # author D (label C)


def test_cif_check_rejects_chain_of_other_protein(tmp_path):
    with pytest.raises(ChainRejected, match="does not map chain A to P25815"):
        verify_chain(_cif_2mjw(tmp_path), "A", "P25815")


def test_cif_check_without_annotation_is_unverified(tmp_path):
    cif = tmp_path / "bare.cif"
    cif.write_text(_strip_loops(_cif_2mjw(tmp_path).read_text(), "_struct_ref."))
    assert "_struct_ref." not in cif.read_text()
    assert verify_chain(cif, "B", "P25815") == "unverified"


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
                                  "Resolution", "Structure method", "SIFTS coverage", "File"]
    assert info.to_dict("records") == [
        {"UniProt ID": "P25815", "Gene Name": "S100P_HUMAN", "PDB": "2MJW", "Chain ID": "B",
         "AlphaFold": "False", "Resolution": "-", "Structure method": "SOLUTION NMR",
         "SIFTS coverage": "94", "File": "02_structures/P25815/2mjw_B.pdb"},
        {"UniProt ID": "P25815", "Gene Name": "S100P_HUMAN", "PDB": "-", "Chain ID": "A",
         "AlphaFold": "True", "Resolution": "-", "Structure method": "-",
         "SIFTS coverage": "-", "File": "02_structures/P25815/AF-P25815-F1.pdb"},
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
                           "af_found": 1, "af_missing": 2, "cif_unverified": 0, "failed": 0}


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
    assert _main(run) == 0
    failed = (layout.structures / "failed.tsv").read_text()
    assert "P99999/9ZZA\tstructures\t" in failed and "404" in failed
    m = json.loads((layout.structures / "manifest.json").read_text())
    assert m["counts"]["failed"] == 1
    assert m["counts"]["chains_extracted"] == 1
    assert structures.units(layout, None).count == 1     # retried on the next run


@pytest.mark.parametrize("broken, reason", [
    (_cif_relettered, "cif annotation does not map chain B to P25815"),
    (_cif_hetatm_only, "chain has no ATOM records"),
])
def test_rejected_chain_recorded_and_not_retried(tmp_path, monkeypatch, broken, reason):
    s = FakeSession(cifs={"2MJW": broken(_cif_2mjw(tmp_path).read_text())})
    run, layout = _setup(tmp_path, monkeypatch, [("P25815", "S100P_HUMAN", "2MJW")], s)
    assert _main(run) == 0
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
    assert _main(run) == 0
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
