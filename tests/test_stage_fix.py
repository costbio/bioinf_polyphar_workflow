import json
import multiprocessing
import os
import shutil
import time
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
    assert stale.read_text() == "ATOM      1  N   MET A   1 truncat"   # harmless leftover, ignored


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
    assert m["args"]["h_seed"] == fix.H_SEED
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
    # Orphaned workers may still be running: wait (bounded, failing on expiry) for them to
    # exit. Each cleans up its unique temp file or finishes and renames it.
    deadline = time.monotonic() + 60
    while multiprocessing.active_children() and time.monotonic() < deadline:
        time.sleep(0.1)
    assert not multiprocessing.active_children(), "worker processes still alive after 60 s"
    assert not list(layout.pdbfixer.rglob("*.tmp"))
    for out in layout.pdbfixer.rglob("*.pdb"):
        assert len(_atoms(out)) > 100 and out.read_text().rstrip().endswith("END")
    counts = fix.run(layout, _args(2), Prog())
    assert counts["fixed"] + counts["skipped_existing"] == 4 and counts["failed"] == 0
    assert not list(layout.pdbfixer.rglob("*.tmp"))


def test_concurrent_fix_one_same_output(tmp_path, small_pdb_text):
    from concurrent.futures import ThreadPoolExecutor

    src = tmp_path / "02_structures" / "P25815" / "2mjw_B.pdb"
    src.parent.mkdir(parents=True)
    src.write_text(small_pdb_text)
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(fix_one, src, tmp_path / "02_structures", tmp_path / "out")
                for _ in range(4)]
        for f in futs:
            f.result()   # none raises
    outs = list((tmp_path / "out").rglob("*"))
    assert [p.name for p in outs if p.is_file()] == ["2mjw_B.pdb"]
    out = tmp_path / "out" / "P25815" / "2mjw_B.pdb"
    assert len(_atoms(out)) > 100 and out.read_text().rstrip().endswith("END")
    assert oct(out.stat().st_mode & 0o777) == "0o644"


def test_worker_crash_is_attributed_and_stage_continues(tmp_path, small_pdb_text, monkeypatch):
    texts = {f"P25815/2mjw_{c}.pdb": small_pdb_text for c in "BCDE"}
    layout = _setup(tmp_path, texts)
    _write_info(layout, list(texts))
    monkeypatch.setenv(fix._CRASH_ENV, "2mjw_C.pdb")   # inherited by spawned workers
    monkeypatch.setattr(fix, "_new_token", lambda: "tok")
    # temp a killed sibling of the crasher would leave (this run's token), and a foreign one
    out_dir = layout.pdbfixer / "P25815"
    out_dir.mkdir(parents=True)
    ours = out_dir / "2mjw_C.pdb.tok.killed.tmp"
    foreign = out_dir / "2mjw_C.pdb.othertok.killed.tmp"
    ours.write_text("ATOM trunc")
    foreign.write_text("ATOM trunc")
    prog = Prog()
    counts = fix.run(layout, _args(2), prog)
    assert counts == {"inputs": 4, "fixed": 3, "skipped_existing": 0, "failed": 1}
    assert len(prog.calls) == 4
    rows = (layout.pdbfixer / "failed.tsv").read_text().splitlines()
    assert rows == ["item\tstage\treason", "P25815/2mjw_C.pdb\tfix\t" + fix.REASON_DIED]
    for c in "BDE":
        assert len(_atoms(layout.pdbfixer / "P25815" / f"2mjw_{c}.pdb")) > 100
    assert not (layout.pdbfixer / "P25815" / "2mjw_C.pdb").exists()
    # the run's own leftovers of units redone in isolation are swept; other runs' are not
    assert list(layout.pdbfixer.rglob("*.tmp")) == [foreign]
    assert (layout.structures / "P25815" / "2mjw_C.pdb").read_text() == small_pdb_text
    assert (layout.pdbfixer / "manifest.json").is_file()
    assert fix.units(layout, _args()).count == 1   # a rerun only retries the crasher


def test_submit_race_routes_remaining_units_to_isolation(tmp_path, small_pdb_text, monkeypatch):
    import concurrent.futures as cf
    from concurrent.futures.process import BrokenProcessPool

    class BreakingPool(cf.ProcessPoolExecutor):
        """The main pool breaks while units are still being submitted (3rd submit on)."""
        n = 0

        def submit(self, *a, **kw):
            if self._max_workers > 1:
                BreakingPool.n += 1
                if BreakingPool.n >= 3:
                    raise BrokenProcessPool("simulated")
            return super().submit(*a, **kw)

    monkeypatch.setattr(cf, "ProcessPoolExecutor", BreakingPool)
    texts = {f"P25815/2mjw_{c}.pdb": small_pdb_text for c in "BCDE"}
    layout = _setup(tmp_path, texts)
    _write_info(layout, list(texts))
    prog = Prog()
    counts = fix.run(layout, _args(2), prog)
    assert counts == {"inputs": 4, "fixed": 4, "skipped_existing": 0, "failed": 0}
    assert len(prog.calls) == 4
    assert BreakingPool.n >= 3   # the race really was hit
    assert len(list(layout.pdbfixer.rglob("*.pdb"))) == 4
    assert not list(layout.pdbfixer.rglob("*.tmp"))


def test_worker_exception_crosses_process_boundary_as_message(tmp_path):
    with pytest.raises(fix.FixError, match="^ValueError|^Exception|^\\w+Error"):
        fix._worker(tmp_path / "missing.pdb", tmp_path, tmp_path / "out", 7.0)


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


# --- reproducibility -------------------------------------------------------------------

def _fix_bytes(tmp_dir, text):
    """Run fix_one on `text` in tmp_dir; return the output bytes without OpenMM's dated header."""
    src = Path(tmp_dir) / "structures" / "P25815" / "2mjw_B.pdb"
    src.parent.mkdir(parents=True)
    src.write_text(text)
    fix_one(src, Path(tmp_dir) / "structures", Path(tmp_dir) / "out")
    out = (Path(tmp_dir) / "out" / "P25815" / "2mjw_B.pdb").read_text()
    return "".join(l for l in out.splitlines(keepends=True)
                   if not l.startswith("REMARK   1 CREATED WITH OPENMM")).encode()


def _fix_bytes_in_child(tmp_dir, text):   # module level: picklable for the spawn pool
    return _fix_bytes(tmp_dir, text)


def _drop_residue_atoms(text):
    """Strip the side chain of one residue so pdbfixer also has missing atoms to add."""
    out, seen = [], 0
    for line in text.splitlines(keepends=True):
        if line.startswith("ATOM") and line[22:27].strip() == "5" and line[12:16].strip() not in ("N", "CA", "C", "O"):
            seen += 1
            continue
        out.append(line)
    assert seen, "fixture has no residue 5 side chain"
    return "".join(out)


def test_fix_one_is_reproducible_in_process(tmp_path, small_pdb_text):
    text = _drop_residue_atoms(small_pdb_text)
    a = _fix_bytes(tmp_path / "a", text)
    b = _fix_bytes(tmp_path / "b", text)
    assert any(l[76:78].strip() == "H" for l in a.decode().splitlines() if l.startswith("ATOM"))
    assert a == b


def test_fix_one_is_reproducible_across_processes(tmp_path, small_pdb_text):
    from concurrent.futures import ProcessPoolExecutor

    text = _drop_residue_atoms(small_pdb_text)
    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=1, mp_context=ctx) as ex1:
        a = ex1.submit(_fix_bytes_in_child, tmp_path / "a", text).result()
    with ProcessPoolExecutor(max_workers=1, mp_context=ctx) as ex2:
        b = ex2.submit(_fix_bytes_in_child, tmp_path / "b", text).result()
    assert a == b


def test_fix_one_leaves_global_random_state_alone(tmp_path, small_pdb_text):
    import random

    random.seed(123)
    expected = [random.random() for _ in range(3)]
    random.seed(123)
    _fix_bytes(tmp_path, small_pdb_text)
    assert [random.random() for _ in range(3)] == expected
