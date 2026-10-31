import json
import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from polyphar import cli
from polyphar.common import RunLayout
from polyphar.stages import fix
from polyphar.stages.fix import fix_one, select_inputs

FIX = Path(__file__).parent / "fixtures"
COLUMNS = ["UniProt ID", "Gene Name", "PDB", "Chain ID", "AlphaFold", "Resolution",
           "Structure method", "SIFTS coverage", "File"]
HETS = ("HETATM 9001 CA    CA B 201      10.000  10.000  10.000  1.00  0.00          CA\n"
        "HETATM 9002  O   HOH B 301      12.000  12.000  12.000  1.00  0.00           O\n")


def _model1(text, hets=""):
    """First model of the 2mjw_B fixture (keeps tests fast), optionally with HETATMs."""
    out = []
    for line in text.splitlines(keepends=True):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("MODEL", "REMARK")):
            out.append(line)
    return "".join(out) + hets + "END\n"


@pytest.fixture(scope="module")
def small_pdb_text():
    return _model1((FIX / "2mjw_B.pdb").read_text())


def _setup(tmp_path, files):
    """files: {relative name: text}; every PDB-like entry is listed in Protein_info.csv,
    AlphaFold rows for names starting with AF-. Returns layout."""
    layout = RunLayout(tmp_path)
    for rel, text in files.items():
        p = layout.structures / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    return layout


def _write_info(layout, listed):
    rows = []
    for rel in listed:
        acc, name = rel.split("/")
        af = name.startswith("AF-")
        rows.append([acc, "G", "-" if af else "2MJW", "A" if af else "B", af, "-", "-", "-",
                     f"02_structures/{rel}"])
    pd.DataFrame(rows, columns=COLUMNS).to_csv(layout.protein_info, index=False)


def _args(nprocs=2):
    return SimpleNamespace(nprocs=nprocs, ph=7.0)


class Prog:
    def __init__(self):
        self.calls = []

    def advance(self, n=1, size=None):
        self.calls.append((n, size))


def _atoms(path):
    return [l for l in Path(path).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]


def test_select_inputs_excludes_af(tmp_path, small_pdb_text):
    layout = _setup(tmp_path, {"P25815/2mjw_B.pdb": small_pdb_text,
                               "P25815/AF-P25815-F1.pdb": small_pdb_text,
                               "P25815/stray_A.pdb": small_pdb_text})
    _write_info(layout, ["P25815/2mjw_B.pdb", "P25815/AF-P25815-F1.pdb"])
    assert select_inputs(layout) == [layout.structures / "P25815" / "2mjw_B.pdb"]


def test_fix_adds_hydrogens_and_drops_hets(tmp_path, small_pdb_text):
    src = tmp_path / "02_structures" / "P25815" / "2mjw_B.pdb"
    src.parent.mkdir(parents=True)
    heavy = "".join(l for l in small_pdb_text.splitlines(keepends=True)
                    if l[76:78].strip() != "H" and l != "END\n")   # NMR model: drop its H
    src.write_text(heavy + HETS + "END\n")
    out_dir = tmp_path / "03_pdbfixer"
    fix_one(src, tmp_path / "02_structures", out_dir)
    out = out_dir / "P25815" / "2mjw_B.pdb"
    atoms = _atoms(out)
    elements = [l[76:78].strip() for l in atoms]
    assert not any(l[76:78].strip() == "H" for l in _atoms(src))
    assert elements.count("H") > 100
    assert not any(" CA    CA" in l or "HOH" in l for l in atoms)
    assert not any(l.startswith("HETATM") for l in atoms)
    assert not out.with_name(out.name + ".tmp").exists()


def test_failure_keeps_input(tmp_path, small_pdb_text):
    bad = "this is not a pdb file\nATOM      1  N   garbage\n"
    layout = _setup(tmp_path, {"P25815/2mjw_B.pdb": small_pdb_text, "P25815/bad_B.pdb": bad})
    _write_info(layout, ["P25815/2mjw_B.pdb", "P25815/bad_B.pdb"])
    before = {p: p.read_bytes() for p in layout.structures.rglob("*.pdb")}
    counts = fix.run(layout, _args(), Prog())
    assert counts["failed"] == 1 and counts["fixed"] == 1 and counts["inputs"] == 2
    assert {p: p.read_bytes() for p in layout.structures.rglob("*.pdb")} == before
    rows = (layout.pdbfixer / "failed.tsv").read_text().splitlines()
    assert rows[0] == "item\tstage\treason"
    assert len(rows) == 2
    item, stage, reason = rows[1].split("\t")
    assert (item, stage) == ("P25815/bad_B.pdb", "fix") and reason
    assert not (layout.pdbfixer / "P25815" / "bad_B.pdb").exists()
    assert not list(layout.pdbfixer.rglob("*.tmp"))
    assert (layout.pdbfixer / "P25815" / "2mjw_B.pdb").is_file()


def test_resume_ignores_partial_tmp(tmp_path, small_pdb_text):
    layout = _setup(tmp_path, {"P25815/2mjw_B.pdb": small_pdb_text,
                               "P25815/2mjw_C.pdb": small_pdb_text})
    _write_info(layout, ["P25815/2mjw_B.pdb", "P25815/2mjw_C.pdb"])
    out_dir = layout.pdbfixer / "P25815"
    out_dir.mkdir(parents=True)
    done = out_dir / "2mjw_B.pdb"
    done.write_text("ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
    stale = out_dir / "2mjw_C.pdb.tmp"
    stale.write_text("ATOM      1  N   MET A   1 truncat")

    wu = fix.units(layout, _args())
    in_c = layout.structures / "P25815" / "2mjw_C.pdb"
    assert wu.count == 1 and wu.cost_basis == "bytes" and wu.sizes == [in_c.stat().st_size]

    prog = Prog()
    counts = fix.run(layout, _args(), prog)
    assert counts == {"inputs": 2, "fixed": 1, "skipped_existing": 1, "failed": 0}
    assert prog.calls == [(1, in_c.stat().st_size)]
    assert len(_atoms(done)) == 1                     # completed output untouched
    assert len(_atoms(out_dir / "2mjw_C.pdb")) > 100  # leftover .tmp was redone
    assert not stale.exists()


def test_empty_output_is_not_done(tmp_path, small_pdb_text):
    layout = _setup(tmp_path, {"P25815/2mjw_B.pdb": small_pdb_text})
    _write_info(layout, ["P25815/2mjw_B.pdb"])
    out = layout.pdbfixer / "P25815" / "2mjw_B.pdb"
    out.parent.mkdir(parents=True)
    out.write_text("")
    assert fix.units(layout, _args()).count == 1


def test_nprocs_gives_identical_output(tmp_path, small_pdb_text):
    texts = {f"P25815/2mjw_{c}.pdb": small_pdb_text for c in "BCD"}
    results = []
    for n in (1, 3):
        layout = _setup(tmp_path / f"n{n}", texts)
        _write_info(layout, list(texts))
        fix.run(layout, _args(n), Prog())
        results.append({p.name: p.read_bytes() for p in layout.pdbfixer.rglob("*.pdb")})
    assert results[0] == results[1] and len(results[0]) == 3


def test_manifest_written_last(tmp_path, small_pdb_text):
    layout = _setup(tmp_path, {"P25815/2mjw_B.pdb": small_pdb_text})
    _write_info(layout, ["P25815/2mjw_B.pdb"])
    fix.run(layout, _args(1), Prog())
    m = json.loads((layout.pdbfixer / "manifest.json").read_text())
    assert m["stage"] == "fix"
    assert m["args"]["ph"] == 7.0
    assert m["versions"]["pdbfixer"] and m["versions"]["openmm"]
    assert m["counts"] == {"inputs": 1, "fixed": 1, "skipped_existing": 0, "failed": 0}
    assert m["inputs"][0]["path"].endswith("Protein_info.csv")


def test_interrupt_propagates_and_leaves_no_partial(tmp_path, small_pdb_text):
    texts = {f"P25815/2mjw_{c}.pdb": small_pdb_text for c in "BCDE"}
    layout = _setup(tmp_path, texts)
    _write_info(layout, list(texts))

    class Boom(Prog):
        def advance(self, n=1, size=None):
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        fix.run(layout, _args(2), Boom())
    assert not (layout.pdbfixer / "manifest.json").exists()
    assert not list(layout.pdbfixer.rglob("*.tmp"))
    # every output present is complete, so a rerun skips it or redoes the rest
    counts = fix.run(layout, _args(2), Prog())
    assert counts["fixed"] + counts["skipped_existing"] == 4 and counts["failed"] == 0


def test_missing_protein_info_is_clear(tmp_path):
    with pytest.raises(FileNotFoundError, match="polyphar structures"):
        fix.units(RunLayout(tmp_path), _args())


def test_cli_dry_run(tmp_path, small_pdb_text, capsys):
    layout = _setup(tmp_path, {"P25815/2mjw_B.pdb": small_pdb_text})
    _write_info(layout, ["P25815/2mjw_B.pdb"])
    assert cli.main(["fix", "--run-dir", str(tmp_path), "--dry-run", "--no-notify",
                     "--ph", "7.4"]) == 0
    assert "fix: 1 units (bytes)" in capsys.readouterr().out
    assert not layout.pdbfixer.exists()


@pytest.mark.regression
def test_regression_1q4k_B_matches_2023(tmp_path):
    base = Path(os.environ["ADEP_DATA"]) / "disease_subtypes/Liver_cancer/net2"
    src = base / "structures/P53350/1q4k_B.pdb"
    ref = base / "pdbfixer/P53350/1q4k_B.pdb"
    if not (src.is_file() and ref.is_file()):
        pytest.skip("1q4k_B reference data not present")
    work = tmp_path / "structures" / "P53350" / "1q4k_B.pdb"
    work.parent.mkdir(parents=True)
    shutil.copy(src, work)
    fix_one(work, tmp_path / "structures", tmp_path / "out")

    def heavy(path):
        d = {}
        for l in _atoms(path):
            if l[76:78].strip() == "H":
                continue
            d[(l[21], l[22:27].strip(), l[17:20].strip(), l[12:16].strip())] = (
                round(float(l[30:38]), 3), round(float(l[38:46]), 3), round(float(l[46:54]), 3))
        return d

    new, old = heavy(tmp_path / "out" / "P53350" / "1q4k_B.pdb"), heavy(ref)
    assert len(new) == len(old) and new.keys() == old.keys()
    # Input atoms are identical to 3 decimals. The terminal OXT is built by pdbfixer, and its
    # position differs by ~0.2 A between the 2023 and the current OpenMM/pdbfixer version.
    moved = {k for k in new if new[k] != old[k]}
    assert all(k[3] == "OXT" for k in moved), moved
    for k in moved:
        assert sum((a - b) ** 2 for a, b in zip(new[k], old[k])) ** 0.5 < 0.5
