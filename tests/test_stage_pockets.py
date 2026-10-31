import json
import os
import shutil
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from polyphar.common import RunLayout, sha256_file
from polyphar.stages import pockets

FIX = Path(__file__).parent / "fixtures"
NET2 = FIX / "pockets_net2_P53350"
FAKE_DATA = FIX / "fake_dogsite_data"
INFO_COLUMNS = ["UniProt ID", "Gene Name", "PDB", "Chain ID", "AlphaFold", "Resolution",
                "Structure method", "SIFTS coverage", "File"]
LEGACY_HEADER = (
    ",UniProt ID,Gene Name,Protein,Pocket,lig_cov,poc_cov,lig_name,4A_crit,volume,hull,"
    "surface,lid,depth,surf/vol,lid/hull,ellVol,ell c/a,ell b/a,siteAtms,accept,donor,aromat,"
    "hydrophobicity,metal,Cs,Ns,Os,Ss,Xs,negAA,posAA,polarAA,apolarAA,sumAA,ALA,ARG,ASN,ASP,"
    "CSO,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,A,C,G,U,I,N,DA,DC,"
    "DG,DT,DN,UNK,simpleScore,drugScore,Pocket Filename")
WARNING = ("No DoGSiteScorer worker is consuming the queue — "
           "start bin/polyphar-dogsite-worker on the host")


def _row(acc, gene, pdb, chain, af):
    if af:
        return [acc, gene, "-", "A", True, "-", "-", "-", f"02_structures/{acc}/AF-{acc}-F1.pdb"]
    return [acc, gene, pdb.upper(), chain, False, "2.0", "X-RAY", "0.9",
            f"02_structures/{acc}/{pdb.lower()}_{chain}.pdb"]


def _write_info(layout, rows):
    layout.structures.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=INFO_COLUMNS).to_csv(layout.protein_info, index=False)


def _put(path, text="ATOM\n"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _model_run(tmp_path, *, missing_fixed=False):
    """Run folder with one fixed PDB chain, one AlphaFold model and (optionally) a
    PDB row whose pdbfixer output is missing."""
    layout = RunLayout(tmp_path)
    rows = [_row("P25815", "S100P_HUMAN", "2mjw", "B", False),
            _row("P25815", "S100P_HUMAN", "-", "A", True)]
    _put(layout.structures / "P25815" / "2mjw_B.pdb", "ATOM raw\n")
    _put(layout.pdbfixer / "P25815" / "2mjw_B.pdb", "ATOM fixed 12345\n")
    _put(layout.structures / "P25815" / "AF-P25815-F1.pdb", "ATOM alphafold model\n")
    if missing_fixed:
        rows.append(_row("P53350", "PLK1_HUMAN", "2ogq", "A", False))
        _put(layout.structures / "P53350" / "2ogq_A.pdb")
    _write_info(layout, rows)
    return layout


class Prog:
    def __init__(self):
        self.calls = []
        self.notes = []

        def notifier(title, message, priority=0):
            self.notes.append((title, message, priority))
            return True
        self.notifier = notifier

    def advance(self, n=1, size=None):
        self.calls.append((n, size))


def _net2_run(tmp_path):
    """The real net2 pocket outputs (desc + exact res names) in 04_pockets/."""
    layout = RunLayout(tmp_path)
    shutil.copytree(NET2 / "P53350", layout.pockets / "P53350")
    shutil.copytree(NET2 / "O14746", layout.pockets / "O14746")
    _write_info(layout, [_row("P53350", "PLK1_HUMAN", "2ogq", "A", False),
                         _row("P53350", "PLK1_HUMAN", "4o56", "A", False),
                         _row("P53350", "PLK1_HUMAN", "4h5x", "B", False),
                         _row("O14746", "TERT_HUMAN", "-", "A", True)])
    return layout


# --- naming / parsing -------------------------------------------------------------

def test_pocket_filename_preserves_case(tmp_path):
    layout = _net2_run(tmp_path)
    df = pockets.parse_pockets(layout.pockets, layout.protein_info)
    assert len(df) == 221
    on_disk = {p.name for p in layout.pockets.rglob("*_res_*.pdb")}
    names = set(df["Pocket Filename"])
    assert "P53350_2ogq_A.pdb_res_P_1.pdb" in names
    assert "P53350_4h5x_B.pdb_res_P_2_1.pdb" in names
    assert names == on_disk      # exact case, every file referenced and present
    assert "P53350_2ogq_a.pdb_res_P_1.pdb" not in names   # the legacy lowercasing

    # columns and values match the legacy DrugScore.csv (except the case-fixed filename)
    out = tmp_path / "DrugScore.csv"
    df.to_csv(out)
    assert out.read_text().splitlines()[0] == LEGACY_HEADER
    legacy = pd.read_csv(NET2 / "legacy_DrugScore_rows.csv", index_col=0)
    legacy["Protein"] = legacy["Protein"].replace("AF-O14746", "AF-O14746-F1")
    key = ["UniProt ID", "Protein", "Pocket"]
    merged = df.merge(legacy, on=key, suffixes=("", "_legacy"), validate="one_to_one")
    assert len(merged) == 221
    assert (merged["Gene Name"] == merged["Gene Name_legacy"]).all()
    for col in ("volume", "simpleScore", "drugScore", "siteAtms"):
        assert (merged[col] - merged[f"{col}_legacy"]).abs().max() < 1e-9, col
    assert (merged["Pocket Filename"].str.lower()
            == merged["Pocket Filename_legacy"].str.lower()).all()


def test_parse_pockets_missing_file_raises(tmp_path):
    layout = _net2_run(tmp_path)
    (layout.pockets / "P53350" / "P53350_4o56_A.pdb_res_P_3.pdb").unlink()
    with pytest.raises(FileNotFoundError, match="P53350_4o56_A.pdb_res_P_3.pdb"):
        pockets.parse_pockets(layout.pockets, layout.protein_info)


def test_parse_pockets_wrong_case_on_disk_raises(tmp_path):
    # Review Focus 2: a file that exists only under another case is not found,
    # even on a case-insensitive filesystem.
    layout = _net2_run(tmp_path)
    p = layout.pockets / "P53350" / "P53350_2ogq_A.pdb_res_P_1.pdb"
    p.rename(p.with_name("P53350_2ogq_a.pdb_res_P_1.pdb"))
    with pytest.raises(FileNotFoundError, match="P53350_2ogq_A.pdb_res_P_1.pdb"):
        pockets.parse_pockets(layout.pockets, layout.protein_info)


def test_af_rows_named_without_prefix(tmp_path):
    layout = _model_run(tmp_path)
    models, missing = pockets.select_models(layout)
    assert missing == []
    by_id = {m.job_id: m for m in models}
    assert set(by_id) == {"P25815_2mjw_B", "AF-P25815-F1"}
    af = by_id["AF-P25815-F1"]
    assert af.model == "02_structures/P25815/AF-P25815-F1.pdb"
    assert af.out == "04_pockets/P25815/AF-P25815-F1.pdb"     # no "P25815_" prefix
    pdb = by_id["P25815_2mjw_B"]
    assert pdb.model == "03_pdbfixer/P25815/2mjw_B.pdb"
    assert pdb.out == "04_pockets/P25815/P25815_2mjw_B.pdb"

    # parsed: Protein is the AF model id, the filename the actual stem
    layout2 = _net2_run(tmp_path / "net2")
    df = pockets.parse_pockets(layout2.pockets, layout2.protein_info)
    af_rows = df[df["UniProt ID"] == "O14746"]
    assert set(af_rows["Protein"]) == {"AF-O14746-F1"}
    assert set(af_rows["Gene Name"]) == {"TERT_HUMAN"}
    assert "AF-O14746-F1.pdb_res_P_12_1.pdb" in set(af_rows["Pocket Filename"])
    assert set(df.loc[df["UniProt ID"] == "P53350", "Protein"]) == {"2OGQ_A", "4O56_A", "4H5X_B"}


# --- queue (container side) ---------------------------------------------------------

def test_enqueue_idempotent_and_relative_paths(tmp_path):
    layout = _model_run(tmp_path, missing_fixed=True)
    q = layout.pockets / "queue"
    q.mkdir(parents=True)
    (q / "CLOSED").write_text("")                    # left over from an earlier run
    first = pockets.enqueue(layout)
    assert sorted(first) == ["AF-P25815-F1", "P25815_2mjw_B"]
    assert not (q / "CLOSED").exists()
    snapshot = {p.name: p.read_text() for p in q.iterdir()}
    second = pockets.enqueue(layout)
    assert sorted(second) == sorted(first)
    assert {p.name: p.read_text() for p in q.iterdir()} == snapshot
    assert sorted(snapshot) == ["AF-P25815-F1.job", "P25815_2mjw_B.job"]   # no temp files
    job = json.loads((q / "P25815_2mjw_B.job").read_text())
    assert job == {"model": "03_pdbfixer/P25815/2mjw_B.pdb",
                   "out": "04_pockets/P25815/P25815_2mjw_B.pdb",
                   "size": (layout.pdbfixer / "P25815" / "2mjw_B.pdb").stat().st_size}
    af = json.loads((q / "AF-P25815-F1.job").read_text())
    assert af["model"] == "02_structures/P25815/AF-P25815-F1.pdb"
    assert af["out"] == "04_pockets/P25815/AF-P25815-F1.pdb"
    for spec in (job, af):
        assert not os.path.isabs(spec["model"]) and not os.path.isabs(spec["out"])

    # a claimed (running) job is not queued again
    (q / "P25815_2mjw_B.job").rename(q / "P25815_2mjw_B.running.somehost.4242")
    assert sorted(pockets.enqueue(layout)) == sorted(first)
    assert not (q / "P25815_2mjw_B.job").exists()


def test_enqueue_requeues_failed_and_stale_done(tmp_path):
    layout = _model_run(tmp_path)
    q = layout.pockets / "queue"
    q.mkdir(parents=True)
    (q / "P25815_2mjw_B.failed").write_text('{"exit": 3, "stderr_tail": "x"}')
    (q / "AF-P25815-F1.done").write_text('{"exit": 0}')   # but no _desc.txt on disk
    assert sorted(pockets.enqueue(layout)) == ["AF-P25815-F1", "P25815_2mjw_B"]
    assert sorted(p.name for p in q.iterdir()) == ["AF-P25815-F1.job", "P25815_2mjw_B.job"]


def _finish(layout, job_id, out, sha="abc"):
    """What the host worker leaves behind for one successful job."""
    _put(layout.run_dir / f"{out}_desc.txt", (FAKE_DATA / "desc.txt").read_text())
    shutil.copy(FAKE_DATA / "res_P_1.pdb", layout.run_dir / f"{out}_res_P_1.pdb")
    _put(layout.pockets / "queue" / f"{job_id}.done",
         json.dumps({"exit": 0, "seconds": 1.0, "dogsite_sha256": sha}))


def test_resume_counts_only_missing_models(tmp_path):
    layout = _model_run(tmp_path, missing_fixed=True)
    wu = pockets.units(layout, SimpleNamespace(nprocs=2))
    fixed = (layout.pdbfixer / "P25815" / "2mjw_B.pdb").stat().st_size
    af = (layout.structures / "P25815" / "AF-P25815-F1.pdb").stat().st_size
    assert wu.count == 2 and wu.cost_basis == "bytes" and sorted(wu.sizes) == sorted([fixed, af])

    _finish(layout, "P25815_2mjw_B", "04_pockets/P25815/P25815_2mjw_B.pdb")
    wu = pockets.units(layout, SimpleNamespace(nprocs=2))
    assert (wu.count, wu.sizes) == (1, [af])
    assert pockets.enqueue(layout) == ["AF-P25815-F1"]


def test_wait_warns_once_without_worker(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    jobs = pockets.enqueue(layout)
    models = {m.job_id: m for m in pockets.select_models(layout)[0]}

    def late_worker():
        time.sleep(1.8)
        for j in jobs:
            _finish(layout, j, models[j].out)
    t = threading.Thread(target=late_worker)
    t.start()
    prog = Prog()
    try:
        res = pockets.wait_for_queue(layout, prog, worker_timeout=1, jobs=jobs)
    finally:
        t.join()
    assert res == {"done": 2, "failed": 0}
    assert len(prog.notes) == 1
    _title, message, priority = prog.notes[0]
    assert (message, priority) == (WARNING, 1)
    assert sorted(s for _, s in prog.calls) == sorted(m.size for m in models.values())
    assert (layout.pockets / "queue" / "CLOSED").exists()


def test_wait_no_warning_with_fresh_heartbeat(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    jobs = pockets.enqueue(layout)
    models = {m.job_id: m for m in pockets.select_models(layout)[0]}
    hb = layout.pockets / "queue" / "WORKER.host.1"

    def worker():
        for _ in range(30):
            hb.touch()
            time.sleep(0.05)
        for j in jobs:
            _finish(layout, j, models[j].out)
    t = threading.Thread(target=worker)
    t.start()
    prog = Prog()
    try:
        pockets.wait_for_queue(layout, prog, worker_timeout=1, jobs=jobs)
    finally:
        t.join()
    assert prog.notes == []


def test_run_records_failures_sha_and_drugscore(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path, missing_fixed=True)
    models = {m.job_id: m for m in pockets.select_models(layout)[0]}

    def worker():
        q = layout.pockets / "queue"
        while not (q / "AF-P25815-F1.job").exists():
            time.sleep(0.02)
        _finish(layout, "P25815_2mjw_B", models["P25815_2mjw_B"].out, sha="s1")
        (q / "AF-P25815-F1.failed").write_text(
            json.dumps({"exit": 3, "stderr_tail": "boom\n"}))
    t = threading.Thread(target=worker)
    t.start()
    try:
        counts = pockets.run(layout, SimpleNamespace(nprocs=2, worker_timeout=600), Prog())
    finally:
        t.join()
    assert counts["failed"] == 2 and counts["done"] == 1 and counts["pockets"] == 1
    failed = (layout.pockets / "failed.tsv").read_text().splitlines()
    assert failed[0] == "item\tstage\treason"
    assert sorted(failed[1:]) == [
        "AF-P25815-F1\tpockets\tdogsite exit 3: boom",
        "P53350_2ogq_A\tpockets\tno pdbfixer output"]
    man = json.loads((layout.pockets / "manifest.json").read_text())
    assert man["versions"]["dogsite_sha256"] == "s1"
    df = pd.read_csv(layout.drugscore_csv, index_col=0)
    assert list(df["Pocket Filename"]) == ["P25815_2mjw_B.pdb_res_P_1.pdb"]
    assert list(df["Protein"]) == ["2MJW_B"]


def test_run_fails_when_dogsite_binaries_differ(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    models = {m.job_id: m for m in pockets.select_models(layout)[0]}
    pockets.enqueue(layout)
    _finish(layout, "P25815_2mjw_B", models["P25815_2mjw_B"].out, sha="s1")
    _finish(layout, "AF-P25815-F1", models["AF-P25815-F1"].out, sha="s2")
    with pytest.raises(RuntimeError, match="DoGSiteScorer"):
        pockets.run(layout, SimpleNamespace(nprocs=2, worker_timeout=600), Prog())
    assert not (layout.pockets / "manifest.json").exists()


def test_cli_exposes_worker_timeout():
    from polyphar import cli
    args = cli.build_parser().parse_args(["pockets", "--worker-timeout", "30"])
    assert args.worker_timeout == 30
    assert cli.build_parser().parse_args(["pockets"]).worker_timeout == 600


# --- real DoGSiteScorer (host only) -------------------------------------------------

REPO = Path(__file__).resolve().parent.parent
WORKER = REPO / "bin" / "polyphar-dogsite-worker"


def _run_worker(run_dir, dogsite, *extra):
    import subprocess
    import sys
    return subprocess.run([sys.executable, str(WORKER), "--run-dir", str(run_dir),
                           "--dogsite", dogsite, *extra], capture_output=True, text=True)


@pytest.mark.dogsite
def test_real_dogsite_on_fixture(tmp_path):
    layout = RunLayout(tmp_path)
    _put(layout.pdbfixer / "P25815" / "2mjw_B.pdb",
         (REPO / "polyphar" / "data" / "bench" / "2mjw_B.pdb").read_text())
    _put(layout.structures / "P25815" / "2mjw_B.pdb")
    _write_info(layout, [_row("P25815", "S100P_HUMAN", "2mjw", "B", False)])
    assert pockets.enqueue(layout) == ["P25815_2mjw_B"]
    r = _run_worker(tmp_path, os.environ["DOGSITE_BIN"])
    assert r.returncode == 0, r.stderr
    done = json.loads((layout.pockets / "queue" / "P25815_2mjw_B.done").read_text())
    assert done["dogsite_sha256"] == sha256_file(shutil.which(os.environ["DOGSITE_BIN"])
                                                 or os.environ["DOGSITE_BIN"])
    df = pockets.parse_pockets(layout.pockets, layout.protein_info)
    assert len(df) > 0
    assert all((layout.pockets / "P25815" / f).is_file() for f in df["Pocket Filename"])


@pytest.mark.dogsite
@pytest.mark.regression
def test_matches_2023_pockets(tmp_path):
    net2 = Path(os.environ["ADEP_DATA"]) / "disease_subtypes" / "Liver_cancer" / "net2"
    legacy = pd.read_csv(net2 / "pockets" / "DrugScore.csv", index_col=0)
    picks = [("P53350", "2ogq_A"), ("P53350", "4o56_A"), ("P53350", "4h5x_B")]
    layout = RunLayout(tmp_path)
    rows = []
    for acc, stem in picks:
        src = net2 / "pdbfixer" / acc / f"{stem}.pdb"
        _put(layout.pdbfixer / acc / f"{stem}.pdb", src.read_text())
        _put(layout.structures / acc / f"{stem}.pdb")
        pdb, chain = stem.split("_")
        rows.append(_row(acc, "PLK1_HUMAN", pdb, chain, False))
    _write_info(layout, rows)
    pockets.enqueue(layout)
    r = _run_worker(tmp_path, os.environ["DOGSITE_BIN"], "--nprocs", "3")
    assert r.returncode == 0, r.stderr
    df = pockets.parse_pockets(layout.pockets, layout.protein_info)
    for acc, stem in picks:
        protein = stem.upper()
        new = df[(df["UniProt ID"] == acc) & (df["Protein"] == protein)].set_index("Pocket")
        old = legacy[(legacy["UniProt ID"] == acc) & (legacy["Protein"] == protein)] \
            .set_index("Pocket")
        assert sorted(new.index) == sorted(old.index), stem
        diff = (new["drugScore"] - old.loc[new.index, "drugScore"]).abs().max()
        assert diff <= 1e-4, (stem, diff)
