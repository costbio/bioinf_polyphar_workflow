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
                "Structure method", "SIFTS coverage", "File", "SIFTS segments",
                "SIFTS other accessions"]
LEGACY_HEADER = (
    ",UniProt ID,Gene Name,Protein,Pocket,lig_cov,poc_cov,lig_name,4A_crit,volume,hull,"
    "surface,lid,depth,surf/vol,lid/hull,ellVol,ell c/a,ell b/a,siteAtms,accept,donor,aromat,"
    "hydrophobicity,metal,Cs,Ns,Os,Ss,Xs,negAA,posAA,polarAA,apolarAA,sumAA,ALA,ARG,ASN,ASP,"
    "CSO,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,A,C,G,U,I,N,DA,DC,"
    "DG,DT,DN,UNK,simpleScore,drugScore,Pocket Filename")
WARNING = ("No DoGSiteScorer worker is consuming the queue — "
           "start bin/polyphar-dogsite-worker on the host")


def _row(acc, gene, pdb, chain, af, segments="1:999", others=""):
    if af:
        return [acc, gene, "-", "A", True, "-", "-", "-", f"02_structures/{acc}/AF-{acc}-F1.pdb",
                "", ""]
    return [acc, gene, pdb.upper(), chain, False, "2.0", "X-RAY", "0.9",
            f"02_structures/{acc}/{pdb.lower()}_{chain}.pdb", segments, others]


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
    df = pockets.parse_pockets(layout.pockets, layout.protein_info, "dogsite2")
    assert len(df) == 221
    on_disk = {p.name for p in layout.pockets.rglob("*_res_*.pdb")}
    names = set(df["Pocket Filename"])
    assert "2ogq_A_res_P_1.pdb" in names
    assert "4h5x_B_res_P_2_1.pdb" in names
    assert names == on_disk      # exact case, every file referenced and present
    assert "2ogq_a_res_P_1.pdb" not in names              # the legacy lowercasing

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
    # the base name changed on 2026-11-02 (no accession prefix, no `.pdb` infix, see
    # out_names); every legacy name still maps onto exactly one new name
    def modernise(acc, name):
        base = name[len(acc) + 1:] if name.lower().startswith(f"{acc.lower()}_") else name
        return base.replace(".pdb_res_", "_res_").lower()

    assert [n.lower() for n in merged["Pocket Filename"]] == [
        modernise(a, n) for a, n in zip(merged["UniProt ID"], merged["Pocket Filename_legacy"])]


def test_parse_pockets_missing_file_raises(tmp_path):
    layout = _net2_run(tmp_path)
    (layout.pockets / "P53350" / "4o56_A_res_P_3.pdb").unlink()
    with pytest.raises(FileNotFoundError, match="4o56_A_res_P_3.pdb"):
        pockets.parse_pockets(layout.pockets, layout.protein_info, "dogsite2")


def test_parse_pockets_wrong_case_on_disk_raises(tmp_path):
    # Review Focus 2: a file that exists only under another case is not found,
    # even on a case-insensitive filesystem.
    layout = _net2_run(tmp_path)
    p = layout.pockets / "P53350" / "2ogq_A_res_P_1.pdb"
    p.rename(p.with_name("2ogq_a_res_P_1.pdb"))
    with pytest.raises(FileNotFoundError, match="2ogq_A_res_P_1.pdb"):
        pockets.parse_pockets(layout.pockets, layout.protein_info, "dogsite2")


def test_af_rows_named_without_prefix(tmp_path):
    layout = _model_run(tmp_path)
    models, missing = pockets.select_models(layout, "dogsite2")
    assert missing == []
    by_id = {m.job_id: m for m in models}
    assert set(by_id) == {"P25815_2mjw_B", "AF-P25815-F1"}
    af = by_id["AF-P25815-F1"]
    assert af.model == "02_structures/P25815/AF-P25815-F1.pdb"
    assert af.out == "04_pockets/P25815/AF-P25815-F1"         # no "P25815_" prefix
    pdb = by_id["P25815_2mjw_B"]
    assert pdb.model == "03_pdbfixer/P25815/2mjw_B.pdb"
    assert pdb.out == "04_pockets/P25815/2mjw_B"   # dogsite2 shares the dogsite3 base name

    # parsed: Protein is the AF model id, the filename the actual stem
    layout2 = _net2_run(tmp_path / "net2")
    df = pockets.parse_pockets(layout2.pockets, layout2.protein_info, "dogsite2")
    af_rows = df[df["UniProt ID"] == "O14746"]
    assert set(af_rows["Protein"]) == {"AF-O14746-F1"}
    assert set(af_rows["Gene Name"]) == {"TERT_HUMAN"}
    assert "AF-O14746-F1_res_P_12_1.pdb" in set(af_rows["Pocket Filename"])
    assert set(df.loc[df["UniProt ID"] == "P53350", "Protein"]) == {"2OGQ_A", "4O56_A", "4H5X_B"}


# --- queue (container side) ---------------------------------------------------------

def test_enqueue_idempotent_and_relative_paths(tmp_path):
    layout = _model_run(tmp_path, missing_fixed=True)
    q = layout.pockets / "queue"
    q.mkdir(parents=True)
    (q / "CLOSED").write_text("")                    # left over from an earlier run
    first = pockets.enqueue(layout, "dogsite2")
    assert sorted(first) == ["AF-P25815-F1", "P25815_2mjw_B"]
    assert not (q / "CLOSED").exists()
    snapshot = {p.name: p.read_text() for p in q.iterdir()}
    second = pockets.enqueue(layout, "dogsite2")
    assert sorted(second) == sorted(first)
    assert {p.name: p.read_text() for p in q.iterdir()} == snapshot
    assert sorted(snapshot) == ["AF-P25815-F1.job", "P25815_2mjw_B.job"]   # no temp files
    job = json.loads((q / "P25815_2mjw_B.job").read_text())
    assert job == {"model": "03_pdbfixer/P25815/2mjw_B.pdb",
                   "out": "04_pockets/P25815/2mjw_B",
                   "size": (layout.pdbfixer / "P25815" / "2mjw_B.pdb").stat().st_size,
                   "tool": "dogsite2"}
    af = json.loads((q / "AF-P25815-F1.job").read_text())
    assert af["model"] == "02_structures/P25815/AF-P25815-F1.pdb"
    assert af["out"] == "04_pockets/P25815/AF-P25815-F1"
    for spec in (job, af):
        assert not os.path.isabs(spec["model"]) and not os.path.isabs(spec["out"])

    # a claimed (running) job is not queued again
    (q / "P25815_2mjw_B.job").rename(q / "P25815_2mjw_B.running.somehost.4242")
    assert sorted(pockets.enqueue(layout, "dogsite2")) == sorted(first)
    assert not (q / "P25815_2mjw_B.job").exists()


def test_enqueue_requeues_failed_and_stale_done(tmp_path):
    layout = _model_run(tmp_path)
    q = layout.pockets / "queue"
    q.mkdir(parents=True)
    (q / "P25815_2mjw_B.failed").write_text('{"exit": 3, "stderr_tail": "x"}')
    (q / "AF-P25815-F1.done").write_text('{"exit": 0}')   # but no _desc.txt on disk
    assert sorted(pockets.enqueue(layout, "dogsite2")) == ["AF-P25815-F1", "P25815_2mjw_B"]
    assert sorted(p.name for p in q.iterdir()) == ["AF-P25815-F1.job", "P25815_2mjw_B.job"]


def _finish(layout, job_id, out, sha="abc"):
    """What the host worker leaves behind for one successful job."""
    pocket = FAKE_DATA / "P53350_4rcp_A"
    _put(layout.run_dir / f"{out}_desc.txt", Path(f"{pocket}.desc.txt").read_text())
    shutil.copy(f"{pocket}.res_P_1.pdb", layout.run_dir / f"{out}_res_P_1.pdb")
    _put(layout.pockets / "queue" / f"{job_id}.done",
         json.dumps({"exit": 0, "seconds": 1.0, "dogsite_sha256": sha}))


def test_resume_counts_only_missing_models(tmp_path):
    layout = _model_run(tmp_path, missing_fixed=True)
    wu = pockets.units(layout, SimpleNamespace(nprocs=2, pocket_tool="dogsite2"))
    fixed = (layout.pdbfixer / "P25815" / "2mjw_B.pdb").stat().st_size
    af = (layout.structures / "P25815" / "AF-P25815-F1.pdb").stat().st_size
    assert wu.count == 2 and wu.cost_basis == "bytes" and sorted(wu.sizes) == sorted([fixed, af])

    _finish(layout, "P25815_2mjw_B", "04_pockets/P25815/2mjw_B")
    wu = pockets.units(layout, SimpleNamespace(nprocs=2, pocket_tool="dogsite2"))
    assert (wu.count, wu.sizes) == (1, [af])
    assert pockets.enqueue(layout, "dogsite2") == ["AF-P25815-F1"]


def test_wait_warns_once_without_worker(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    jobs = pockets.enqueue(layout, "dogsite2")
    models = {m.job_id: m for m in pockets.select_models(layout, "dogsite2")[0]}

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


def test_wait_repeats_the_warning_and_fails_once_the_worker_is_dead(tmp_path, monkeypatch):
    # The warning is emitted on a wall-clock interval from inside the poll loop, so a
    # starved loop emits fewer of them. With a 1.5 s deadline and a 0.4 s interval this
    # test failed under a loaded full-suite run and passed alone; the headroom below
    # keeps what it proves (the warning REPEATS) without depending on scheduling.
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    monkeypatch.setattr(pockets, "NO_WORKER_REPEAT_S", 0.2)
    layout = _model_run(tmp_path)
    jobs = pockets.enqueue(layout, "dogsite2")
    prog = Prog()
    t0 = time.monotonic()
    with pytest.raises(RuntimeError, match="no DoGSiteScorer worker heartbeat for"):
        pockets.wait_for_queue(layout, prog, worker_timeout=0.2, jobs=jobs,
                               worker_dead_after=3.0)
    assert 3.0 <= time.monotonic() - t0 < 12
    assert len(prog.notes) >= 3                      # repeated, not once
    assert all(n[1] == WARNING and n[2] == 1 for n in prog.notes)
    assert not (layout.pockets / "queue" / "CLOSED").exists()
    assert sorted(n for n in os.listdir(layout.pockets / "queue")) == sorted(
        f"{j}.job" for j in jobs)                    # jobs stay queued for a later run


def test_cli_exposes_worker_dead_after():
    from polyphar import cli
    assert cli.build_parser().parse_args(["pockets"]).worker_dead_after == 6 * 3600
    args = cli.build_parser().parse_args(["pockets", "--worker-dead-after", "60"])
    assert args.worker_dead_after == 60


def test_wait_no_warning_with_fresh_heartbeat(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    jobs = pockets.enqueue(layout, "dogsite2")
    models = {m.job_id: m for m in pockets.select_models(layout, "dogsite2")[0]}
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
    models = {m.job_id: m for m in pockets.select_models(layout, "dogsite2")[0]}

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
        counts = pockets.run(layout, SimpleNamespace(nprocs=2, worker_timeout=600, pocket_tool="dogsite2"), Prog())
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
    assert list(df["Pocket Filename"]) == ["2mjw_B_res_P_1.pdb"]
    assert list(df["Protein"]) == ["2MJW_B"]


def test_run_fails_when_dogsite_binaries_differ(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    models = {m.job_id: m for m in pockets.select_models(layout, "dogsite2")[0]}
    pockets.enqueue(layout, "dogsite2")
    _finish(layout, "P25815_2mjw_B", models["P25815_2mjw_B"].out, sha="s1")
    _finish(layout, "AF-P25815-F1", models["AF-P25815-F1"].out, sha="s2")
    with pytest.raises(RuntimeError, match="DoGSiteScorer"):
        pockets.run(layout, SimpleNamespace(nprocs=2, worker_timeout=600, pocket_tool="dogsite2"), Prog())
    assert not (layout.pockets / "manifest.json").exists()


def test_cli_exposes_worker_timeout():
    from polyphar import cli
    args = cli.build_parser().parse_args(["pockets", "--worker-timeout", "30"])
    assert args.worker_timeout == 30
    assert cli.build_parser().parse_args(["pockets"]).worker_timeout == 600


# --- fusion / chimeric chains: pockets outside this accession's SIFTS segment(s) ---------

def _atom(serial, resseq, x, chain="A", icode=" ", name="CA", elem="C", resn="ALA"):
    return ("ATOM  %5d %-4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, name, resn, chain, resseq, icode, x, 0.0, 0.0, elem))


def _chimera_run(tmp_path, segments="1:5", others="P22222", icodes=None):
    """One PDB chain of 10 residues (author numbers 1..10, CA at x = 3*i, plus an H atom);
    SIFTS says P11111 is `segments`, the rest belongs to `others`. Pockets (residues
    renumbered 101.. as pdbfixer + DoGSiteScorer would, matched by coordinates):
    P_1 = 1-4 (inside), P_2 = 6-9 (outside), P_3 = 4-7 (half: kept), P_4 = 5-10 (5/6 out)."""
    layout = RunLayout(tmp_path)
    icodes = icodes or {}
    raw = "MODEL        1\n" + "".join(
        _atom(i, i, 3.0 * i, icode=icodes.get(i, " ")) + _atom(100 + i, i, 3.0 * i + 1, name="H",
                                                               elem="H")
        for i in range(1, 11)) + "ENDMDL\nMODEL        2\n" + "".join(
        _atom(i, 50 + i, 3.0 * i) for i in range(1, 11)) + "ENDMDL\nEND\n"
    _put(layout.structures / "P11111" / "9xyz_A.pdb", raw)
    _put(layout.pdbfixer / "P11111" / "9xyz_A.pdb", "ATOM fixed\n")
    _write_info(layout, [_row("P11111", "FUSN_HUMAN", "9xyz", "A", False, segments, others)])
    d = layout.pockets / "P11111"
    pockets_res = {"P_1": range(1, 5), "P_2": range(6, 10), "P_3": range(4, 8),
                   "P_4": range(5, 11)}
    desc = "name\tvolume\tdrugScore\n" + "".join(f"{p}\t100\t0.8\n" for p in pockets_res)
    _put(d / "9xyz_A_desc.txt", desc)
    for p, res in pockets_res.items():
        text = "HEADER File contains complete residues of binding site atoms\n" + "".join(
            _atom(k, 100 + i, 3.0 * i + 0.0004) + _atom(500 + k, 100 + i, 50.0 + i, name="H",
                                                       elem="H")
            for k, i in enumerate(res, 1))
        _put(d / f"9xyz_A_res_{p}.pdb", text)
    return layout


def test_pockets_outside_the_target_segment_are_dropped(tmp_path):
    layout = _chimera_run(tmp_path)
    counts = pockets.run(layout, SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite2"), Prog())
    df = pd.read_csv(layout.drugscore_csv, index_col=0)
    assert sorted(df["Pocket"]) == ["P_1", "P_3"]          # P_3: 50 % outside is kept
    assert counts["outside_segment"] == 2 and counts["failed"] == 0
    assert counts["pockets"] == 2 and counts["segment_unverified"] == 0
    failed = sorted((layout.pockets / "failed.tsv").read_text().splitlines()[1:])
    assert failed == [
        "P11111/9xyz_A_res_P_2.pdb\tpockets\tpocket outside target segment: "
        "100% outside; partner(s) P22222",
        "P11111/9xyz_A_res_P_4.pdb\tpockets\tpocket outside target segment: "
        "83% outside; partner(s) P22222"]


def test_single_protein_chains_are_not_filtered(tmp_path):
    layout = _chimera_run(tmp_path, others="")
    counts = pockets.run(layout, SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite2"), Prog())
    assert counts["outside_segment"] == 0 and counts["pockets"] == 4


@pytest.mark.parametrize("segments", ["1:None", ":5", "1:5;x:9", ""])
def test_unparseable_segment_bounds_keep_pockets_as_unverified(tmp_path, segments):
    layout = _chimera_run(tmp_path, segments=segments)
    counts = pockets.run(layout, SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite2"), Prog())
    assert counts["pockets"] == 4 and counts["outside_segment"] == 0
    assert counts["segment_unverified"] == 4


def test_insertion_code_bounds_are_parsed_and_lenient(tmp_path):
    # bound 5A: residue 5 (no insertion code) at the same number counts as inside
    layout = _chimera_run(tmp_path, segments="1:5A")
    counts = pockets.run(layout, SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite2"), Prog())
    assert counts["outside_segment"] == 2 and counts["segment_unverified"] == 0


def test_pocket_without_any_matching_atom_is_unverified_and_kept(tmp_path):
    layout = _chimera_run(tmp_path)
    p = layout.pockets / "P11111" / "9xyz_A_res_P_2.pdb"
    p.write_text("".join(_atom(k, 7, 900.0 + k) for k in range(3)))
    counts = pockets.run(layout, SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite2"), Prog())
    assert counts["outside_segment"] == 1 and counts["segment_unverified"] == 1
    assert "P_2" in set(pd.read_csv(layout.drugscore_csv, index_col=0)["Pocket"])


def test_protein_info_without_segment_columns_is_an_error(tmp_path):
    layout = _chimera_run(tmp_path)
    info = pd.read_csv(layout.protein_info, dtype=str, keep_default_na=False)
    info.drop(columns=["SIFTS segments", "SIFTS other accessions"]).to_csv(
        layout.protein_info, index=False)
    with pytest.raises(ValueError, match="rerun `polyphar structures`"):
        pockets.run(layout, SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite2"), Prog())


# --- real DoGSiteScorer (host only) -------------------------------------------------

REPO = Path(__file__).resolve().parent.parent
WORKER = REPO / "bin" / "polyphar-dogsite-worker"


def _run_worker(run_dir, dogsite, *extra):
    import subprocess
    import sys
    return subprocess.run([sys.executable, str(WORKER), "--run-dir", str(run_dir),
                           "--dogsite", dogsite, *extra], capture_output=True, text=True)


def _real_tool():
    """The tool of the real binary in $DOGSITE_BIN: $POLYPHAR_POCKET_TOOL, default dogsite2."""
    return os.environ.get("POLYPHAR_POCKET_TOOL") or "dogsite2"


@pytest.mark.dogsite
def test_real_dogsite_on_fixture(tmp_path):
    tool = _real_tool()
    layout = RunLayout(tmp_path)
    _put(layout.pdbfixer / "P25815" / "2mjw_B.pdb",
         (REPO / "polyphar" / "data" / "bench" / "2mjw_B.pdb").read_text())
    _put(layout.structures / "P25815" / "2mjw_B.pdb")
    _write_info(layout, [_row("P25815", "S100P_HUMAN", "2mjw", "B", False)])
    assert pockets.enqueue(layout, tool) == ["P25815_2mjw_B"]
    r = _run_worker(tmp_path, os.environ["DOGSITE_BIN"], "--tool", tool)
    assert r.returncode == 0, r.stderr
    done = json.loads((layout.pockets / "queue" / "P25815_2mjw_B.done").read_text())
    assert done["dogsite_sha256"] == sha256_file(shutil.which(os.environ["DOGSITE_BIN"])
                                                 or os.environ["DOGSITE_BIN"])
    assert done["tool"] == tool and done["tool_version"]
    df = pockets.parse_pockets(layout.pockets, layout.protein_info, tool)
    assert len(df) > 0
    assert all((layout.pockets / "P25815" / f).is_file() for f in df["Pocket Filename"])
    if tool == "dogsite3":
        assert df["Pocket Filename"].iloc[0] == "2mjw_B_P_1_res.pdb"
        assert "volume" in df.columns and "drugScore" not in df.columns
        assert (layout.pockets / "P25815" / "2mjw_B_P_1_res.pdb").read_text().startswith(
            "HEADER File contains complete residues of binding site atoms")


@pytest.mark.dogsite
def test_real_dogsite3_segment_filter(tmp_path):
    """A real DoGSite3 run through the worker and the stage: the fusion segment filter maps
    every real residue file back to the extracted chain by coordinates."""
    if _real_tool() != "dogsite3":
        pytest.skip("needs the DoGSite3 binary")
    chain = REPO / "polyphar" / "data" / "bench" / "5nmm_A.pdb"
    layout = RunLayout(tmp_path)
    _put(layout.structures / "P53350" / "5nmm_A.pdb", chain.read_text())
    _put(layout.pdbfixer / "P53350" / "5nmm_A.pdb", chain.read_text())
    _write_info(layout, [_row("P53350", "PLK1_HUMAN", "5nmm", "A", False, "-999:9999",
                              "P99999")])
    pockets.enqueue(layout, "dogsite3")
    r = _run_worker(tmp_path, os.environ["DOGSITE_BIN"], "--tool", "dogsite3")
    assert r.returncode == 0, r.stderr
    counts = pockets.run(layout, SimpleNamespace(nprocs=1, worker_timeout=600,
                                                 pocket_tool="dogsite3", min_volume=300.0),
                         Prog())
    assert counts["pockets"] > 0 and counts["outside_segment"] == 0
    assert counts["segment_unverified"] == 0      # every pocket mapped by coordinates
    df = pd.read_csv(layout.drugscore_csv, index_col=0)
    assert (df["Selected"] == (df["volume"] >= 300)).all() and set(df["Tool"]) == {"dogsite3"}


@pytest.mark.dogsite
@pytest.mark.regression
def test_matches_2023_pockets(tmp_path):
    if _real_tool() != "dogsite2":
        pytest.skip("the 2023 pockets are DoGSiteScorer 2 output (POLYPHAR_POCKET_TOOL=dogsite2)")
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
    pockets.enqueue(layout, "dogsite2")
    r = _run_worker(tmp_path, os.environ["DOGSITE_BIN"], "--nprocs", "3")
    assert r.returncode == 0, r.stderr
    df = pockets.parse_pockets(layout.pockets, layout.protein_info, "dogsite2")
    for acc, stem in picks:
        protein = stem.upper()
        new = df[(df["UniProt ID"] == acc) & (df["Protein"] == protein)].set_index("Pocket")
        old = legacy[(legacy["UniProt ID"] == acc) & (legacy["Protein"] == protein)] \
            .set_index("Pocket")
        assert sorted(new.index) == sorted(old.index), stem
        diff = (new["drugScore"] - old.loc[new.index, "drugScore"]).abs().max()
        assert diff <= 1e-4, (stem, diff)


def test_fusion_bounds_recovered_from_seqres_drop_the_partner_pocket(tmp_path):
    """SIFTS gives P11111's segment on the fusion chain with empty PDB_BEG/PDB_END (it used
    to stay segment_unverified, every pocket kept). Its SEQRES range 3..7 maps through the
    cif's _pdbx_poly_seq_scheme to author 1..5 (seq 1-2 are an unobserved tag), so the
    pockets in the partner's part (P_2, P_4) are now dropped."""
    from polyphar.stages import structures
    from test_stage_structures import scheme_cif

    residues = [(1, None, ""), (2, None, "")] + [(s, s - 2, "") for s in range(3, 13)]
    cd = structures._cif_dict(scheme_cif(tmp_path / "9xyz.cif", {"A": residues}))
    segments, source = structures.author_segments("None:None", "3:7", cd, "A")
    assert (segments, source) == ("1:5", "seqres")
    layout = _chimera_run(tmp_path / "run", segments=segments)
    counts = pockets.run(layout, SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite2"), Prog())
    assert sorted(pd.read_csv(layout.drugscore_csv, index_col=0)["Pocket"]) == ["P_1", "P_3"]
    assert counts["outside_segment"] == 2 and counts["segment_unverified"] == 0


# --- DoGSite3 (the sensitivity analysis) -----------------------------------------------

FAKE3 = FIX / "fake_dogsite3_data"
DS3 = SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite3", min_volume=300.0,
                      drugscore_cutoff=0.5)
DS2 = SimpleNamespace(nprocs=1, worker_timeout=600, pocket_tool="dogsite2", min_volume=300.0,
                      drugscore_cutoff=0.5)


def _finish3(layout, job_id, out, label="O60218_4gqg_A", sha="abc", version="1.3.0"):
    """What the host worker leaves behind for one successful DoGSite3 job: a real DoGSite3
    output set renamed to `out`."""
    src = FAKE3 / label
    for f in sorted(os.listdir(src)):
        if f.startswith(label + "_"):
            _put(layout.run_dir / f"{out}{f[len(label):]}", (src / f).read_text())
    _put(layout.pockets / "queue" / f"{job_id}.done",
         json.dumps({"exit": 0, "seconds": 1.0, "dogsite_sha256": sha, "tool": "dogsite3",
                     "tool_version": version}))


def test_default_tool_is_dogsite2_and_dogsite3_uses_the_calibrated_volume(monkeypatch):
    """DoGSiteScorer 2 is the primary tool; DoGSite3 (sensitivity analysis) keeps pockets
    with volume >= 120 A^3 (scripts/calibrate_dogsite3_volume.py)."""
    from polyphar import cli
    monkeypatch.delenv("POLYPHAR_POCKET_TOOL", raising=False)
    args = cli.build_parser().parse_args(["pockets"])
    assert (args.pocket_tool, args.min_volume, args.drugscore_cutoff) == ("dogsite2", 120.0, 0.5)
    args = cli.build_parser().parse_args(["run", "--pocket-tool", "dogsite3",
                                          "--min-volume", "250"])
    assert (args.pocket_tool, args.min_volume) == ("dogsite3", 250.0)
    monkeypatch.setenv("POLYPHAR_POCKET_TOOL", "dogsite3")
    assert cli.build_parser().parse_args(["pockets"]).pocket_tool == "dogsite3"
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["pockets", "--pocket-tool", "fpocket"])


def test_dogsite3_names_and_job_files(tmp_path):
    layout = _model_run(tmp_path)
    models = {m.job_id: m for m in pockets.select_models(layout, "dogsite3")[0]}
    assert set(models) == {"P25815_2mjw_B", "AF-P25815-F1"}      # job ids as before
    assert models["P25815_2mjw_B"].out == "04_pockets/P25815/2mjw_B"   # no UniProt prefix
    assert models["AF-P25815-F1"].out == "04_pockets/P25815/AF-P25815-F1"
    pockets.enqueue(layout, "dogsite3")
    job = json.loads((layout.pockets / "queue" / "P25815_2mjw_B.job").read_text())
    assert job == {"model": "03_pdbfixer/P25815/2mjw_B.pdb", "out": "04_pockets/P25815/2mjw_B",
                   "size": (layout.pdbfixer / "P25815" / "2mjw_B.pdb").stat().st_size,
                   "tool": "dogsite3"}
    assert pockets.pocket_filename("dogsite3", "2mjw_B", "P_1_2") == "2mjw_B_P_1_2_res.pdb"
    assert pockets.pocket_filename("dogsite2", "2mjw_B", "P_1_2") == "2mjw_B_res_P_1_2.pdb"


def test_dogsite3_chain_shared_by_accessions_keeps_the_accession(tmp_path):
    """A chimeric chain extracted for two accessions would give the same short name in both
    folders (PocketMatch keeps only base names): those keep the <UniProt>_ prefix."""
    layout = RunLayout(tmp_path)
    rows = [_row("P11111", "A_HUMAN", "9xyz", "A", False, "1:5", "P22222"),
            _row("P22222", "B_HUMAN", "9xyz", "A", False, "6:10", "P11111"),
            _row("P22222", "B_HUMAN", "1abc", "B", False)]
    for acc, stem in (("P11111", "9xyz_A"), ("P22222", "9xyz_A"), ("P22222", "1abc_B")):
        _put(layout.pdbfixer / acc / f"{stem}.pdb")
    _write_info(layout, rows)
    outs = sorted(m.out for m in pockets.select_models(layout, "dogsite3")[0])
    assert outs == ["04_pockets/P11111/P11111_9xyz_A", "04_pockets/P22222/1abc_B",
                    "04_pockets/P22222/P22222_9xyz_A"]
    assert len({m.job_id for m in pockets.select_models(layout, "dogsite3")[0]}) == 3


def test_dogsite3_longest_realistic_names_fit_pocketmatch():
    from polyphar.stages.similarity import PM_NAME_LIMIT
    # 10-character accession, 2-letter chain, pocket 12 subpocket 1
    assert len(pockets.pocket_filename("dogsite3", "7a8b_AB", "P_12_1")) <= PM_NAME_LIMIT
    assert len(pockets.pocket_filename("dogsite3", "AF-A0A0B4J2F0-F1", "P_12_1")) <= PM_NAME_LIMIT
    # the long form is only used for chains shared by accessions (rare); it is 33 here and
    # the similarity name guard would report it
    assert len(pockets.pocket_filename("dogsite3", "A0A0B4J2F0_7a8b_AB", "P_12_1")) == 33


def test_both_tools_share_the_base_name_and_fit_a_long_chain_id(tmp_path):
    """Decided 2026-11-02: one base name for both tools, `<pdb>_<chain>`. 2023's dogsite2
    name `<UniProt>_<pdb>_<chain>.pdb_res_<pocket>.pdb` is 33 characters for PLK1's 9r1w
    chain AAA (a current PDB entry's multi-character author chain ID), which PocketMatch
    would truncate and merge with another pocket."""
    from polyphar.stages.similarity import PM_NAME_LIMIT
    layout = RunLayout(tmp_path)
    _write_info(layout, [_row("P53350", "PLK1_HUMAN", "9r1w", "AAA", False)])
    for tool in ("dogsite2", "dogsite3"):
        assert pockets.out_names(pd.read_csv(layout.protein_info), tool) == ["9r1w_AAA"]
        name = pockets.pocket_filename(tool, "9r1w_AAA", "P_1_1")
        assert len(name) <= PM_NAME_LIMIT, (tool, name)
        assert pockets._model_of(tool, name) == "9r1w_AAA"     # round-trips
    assert len("P53350_9r1w_AAA.pdb_res_P_1_1.pdb") == 33      # what 2023 would have written


def test_dogsite3_run_writes_tool_selected_and_rule(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    models = {m.job_id: m for m in pockets.select_models(layout, "dogsite3")[0]}

    def worker():
        q = layout.pockets / "queue"
        while not (q / "AF-P25815-F1.job").exists():
            time.sleep(0.02)
        _finish3(layout, "P25815_2mjw_B", models["P25815_2mjw_B"].out, "O60218_4gqg_A")
        _finish3(layout, "AF-P25815-F1", models["AF-P25815-F1"].out, "O60218_4ga8_A")
    t = threading.Thread(target=worker)
    t.start()
    try:
        counts = pockets.run(layout, DS3, Prog())
    finally:
        t.join()
    df = pd.read_csv(layout.drugscore_csv, index_col=0)
    assert list(df.columns[:5]) == ["UniProt ID", "Gene Name", "Protein", "Pocket", "lig_cov"]
    assert list(df.columns[-3:]) == ["Pocket Filename", "Tool", "Selected"]
    assert "drugScore" not in df.columns and "volume" in df.columns
    assert set(df["Tool"]) == {"dogsite3"}
    assert (df["Selected"] == (df["volume"] >= 300)).all()
    sel = df[df["Selected"]]
    assert sorted(sel["Pocket Filename"]) == ["2mjw_B_P_1_1_res.pdb", "2mjw_B_P_1_res.pdb",
                                              "AF-P25815-F1_P_1_res.pdb"]
    on_disk = {p.name for p in (layout.pockets / "P25815").glob("*_res.pdb")}
    assert set(df["Pocket Filename"]) == on_disk
    assert counts["pockets"] == 13 and counts["selected"] == 3
    man = json.loads((layout.pockets / "manifest.json").read_text())
    assert man["versions"]["pocket_tool"] == "dogsite3"
    assert man["versions"]["dogsite_version"] == "1.3.0"
    assert man["selection_rule"] == {"tool": "dogsite3", "column": "volume", "min": 300.0,
                                     "text": "volume >= 300 A^3"}
    assert man["args"]["pocket_tool"] == "dogsite3" and man["args"]["min_volume"] == 300.0


def test_dogsite2_run_marks_selected_by_drugscore(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _net2_run(tmp_path)
    for m in pockets.select_models(layout, "dogsite2")[0]:
        _put(layout.pockets / "queue" / f"{m.job_id}.done",
             json.dumps({"exit": 0, "dogsite_sha256": "x"}))      # legacy .done: no tool
    pockets.run(layout, DS2, Prog())
    df = pd.read_csv(layout.drugscore_csv, index_col=0)
    assert len(df) == 221 and set(df["Tool"]) == {"dogsite2"}
    assert (df["Selected"] == (df["drugScore"] >= 0.5)).all() and df["Selected"].any()
    man = json.loads((layout.pockets / "manifest.json").read_text())
    assert man["selection_rule"]["text"] == "drugScore >= 0.5"


def test_refuses_to_mix_tools_in_one_run_folder(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    # legacy outputs: a DoGSiteScorer 2 .done without a tool field
    m = {x.job_id: x for x in pockets.select_models(layout, "dogsite2")[0]}
    _finish(layout, "P25815_2mjw_B", m["P25815_2mjw_B"].out)
    with pytest.raises(RuntimeError, match=r"dogsite2.*--pocket-tool dogsite2"):
        pockets.units(layout, DS3)
    with pytest.raises(RuntimeError, match="dogsite2"):
        pockets.run(layout, DS3, Prog())
    assert not list((layout.pockets / "queue").glob("*.job"))   # nothing queued
    # the recorded tool also comes from the TOOL marker once a run started
    layout2 = _model_run(tmp_path / "b")
    layout2.pockets.mkdir(parents=True)
    (layout2.pockets / "TOOL").write_text("dogsite3\n")
    assert pockets.recorded_tool(layout2) == "dogsite3"
    with pytest.raises(RuntimeError, match="dogsite3"):
        pockets.units(layout2, DS2)
    assert pockets.units(layout2, DS3).count == 2


def test_selection_change_reruns_the_stage(tmp_path, monkeypatch):
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = _model_run(tmp_path)
    for m in pockets.select_models(layout, "dogsite3")[0]:
        _finish3(layout, m.job_id, m.out)
    pockets.run(layout, DS3, Prog())
    assert pockets.needs_rerun(layout, DS3) is None
    stricter = SimpleNamespace(**{**vars(DS3), "min_volume": 400.0})
    assert "volume >= 300" in pockets.needs_rerun(layout, stricter)
    pockets.run(layout, stricter, Prog())
    df = pd.read_csv(layout.drugscore_csv, index_col=0)
    assert (df["Selected"] == (df["volume"] >= 400)).all() and df["Selected"].sum() == 2


def test_fusion_filter_works_on_real_dogsite3_residue_files(tmp_path):
    """The segment filter maps DoGSite3 residue files back to author numbers by coordinates:
    a real DoGSite3 pocket of 4gqg_A is inside a segment covering the whole chain, and
    outside a segment that excludes every residue."""
    src = FAKE3 / "O60218_4gqg_A"
    for segments, n_out in (("-999:9999", 0), ("9000:9999", 8)):
        layout = RunLayout(tmp_path / segments.replace(":", "_"))
        # the "extracted chain" is made of the residue files' atoms (same coordinates)
        atoms = "".join(l for f in sorted(src.glob("*_res.pdb"))
                        for l in f.read_text().splitlines(True) if l.startswith("ATOM"))
        _put(layout.structures / "O60218" / "4gqg_A.pdb", atoms)
        _put(layout.pdbfixer / "O60218" / "4gqg_A.pdb")
        _write_info(layout, [_row("O60218", "AK1BA_HUMAN", "4gqg", "A", False, segments,
                                  "P99999")])
        m = pockets.select_models(layout, "dogsite3")[0][0]
        assert m.out == "04_pockets/O60218/4gqg_A"
        _finish3(layout, m.job_id, m.out, "O60218_4gqg_A")
        counts = pockets.run(layout, DS3, Prog())
        assert counts["outside_segment"] == n_out and counts["segment_unverified"] == 0


def test_stage_fails_on_selected_names_pocketmatch_cannot_keep(tmp_path, monkeypatch):
    """The similarity name guard also runs at the end of pockets, on the selected pockets:
    a chimeric chain of a 10-character accession can give a 33-character name."""
    monkeypatch.setattr(pockets, "POLL_S", 0.05)
    layout = RunLayout(tmp_path)
    rows = [_row("A0A0B4J2F0", "X_HUMAN", "7a8b", "AB", False, "1:999", "P22222"),
            _row("P22222", "Y_HUMAN", "7a8b", "AB", False, "1:999", "A0A0B4J2F0")]
    for acc in ("A0A0B4J2F0", "P22222"):
        _put(layout.pdbfixer / acc / "7a8b_AB.pdb")
        _put(layout.structures / acc / "7a8b_AB.pdb")
    _write_info(layout, rows)
    for m in pockets.select_models(layout, "dogsite3")[0]:
        _put(layout.run_dir / f"{m.out}_desc.txt", "name\tvolume\nP_12_1\t400\nP_13_1\t100\n")
        for p in ("P_12_1", "P_13_1"):
            _put(layout.run_dir / f"{m.out}_{p}_res.pdb", "HEADER pocket\n")
        _put(layout.pockets / "queue" / f"{m.job_id}.done",
             json.dumps({"exit": 0, "tool": "dogsite3", "dogsite_sha256": "x"}))
    with pytest.raises(ValueError, match=r"longer than 32.*A0A0B4J2F0_7a8b_AB_P_12_1_res\.pdb"):
        pockets.run(layout, DS3, Prog())
    failed = (layout.pockets / "failed.tsv").read_text()
    assert "A0A0B4J2F0/A0A0B4J2F0_7a8b_AB_P_12_1_res.pdb\tpockets\tpocket name longer" in failed
    assert "P_13_1" not in failed                    # not selected: never reaches PocketMatch
    assert not (layout.pockets / "manifest.json").exists()


# --- pockets whose lining was never modelled -------------------------------------------

def _pocket_pdb(path, residues):
    """residues = [(resname, [atom names])] -> a pocket PDB."""
    lines, i = [], 0
    for r, (name, atoms) in enumerate(residues, start=1):
        for atom in atoms:
            i += 1
            lines.append(f"ATOM  {i:5d}  {atom:<3} {name} A{r:4d}    "
                         f"{r:8.3f}{i:8.3f}{r:8.3f}  1.00  0.00           {atom[0]}\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "END\n")
    return path


_REAL = ("LEU", ["N", "CA", "C", "O", "CB", "CG", "CD1"])
_UNK = ("UNK", ["N", "CA", "C", "O", "CB"])


def test_a_pocket_lined_by_unmodelled_backbone_is_dropped_with_a_reason(tmp_path):
    """The 6BQ1 shape: every lining residue is traced main chain, so there is no
    chemistry to dock against and PocketMatch hangs on it."""
    layout = RunLayout(tmp_path)
    _pocket_pdb(layout.pockets / "P1" / "good_res_P_1.pdb", [_REAL] * 8)
    _pocket_pdb(layout.pockets / "P2" / "unk_res_P_1.pdb", [_UNK] * 8)
    df = pd.DataFrame([{"UniProt ID": "P1", "Pocket Filename": "good_res_P_1.pdb"},
                       {"UniProt ID": "P2", "Pocket Filename": "unk_res_P_1.pdb"}])

    kept, dropped = pockets.filter_unmodelled(layout, df)
    assert kept["Pocket Filename"].tolist() == ["good_res_P_1.pdb"]
    assert len(dropped) == 1
    item, reason = dropped[0]
    assert item == "P2/unk_res_P_1.pdb"
    assert pockets.REASON_UNMODELLED in reason and "8 of 8 residues" in reason


def test_a_pocket_with_a_minority_of_unmodelled_residues_is_kept(tmp_path):
    """Partial tracing is ordinary; only a lining that is mostly unmodelled is not a site."""
    layout = RunLayout(tmp_path)
    _pocket_pdb(layout.pockets / "P1" / "mixed_res_P_1.pdb", [_REAL] * 6 + [_UNK] * 2)
    df = pd.DataFrame([{"UniProt ID": "P1", "Pocket Filename": "mixed_res_P_1.pdb"}])
    kept, dropped = pockets.filter_unmodelled(layout, df)
    assert len(kept) == 1 and dropped == []


def test_filter_unmodelled_on_an_empty_frame_does_nothing(tmp_path):
    kept, dropped = pockets.filter_unmodelled(RunLayout(tmp_path), pd.DataFrame())
    assert kept.empty and dropped == []
