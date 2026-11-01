import json
import multiprocessing
import os
import shutil
import time
import unittest.mock
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from polyphar import cli, pool
from polyphar.common import RunLayout
from polyphar.stages import similarity as sim

FIX = Path(__file__).parent / "fixtures"
PM_BIN = os.environ.get("PM_BIN_DIR", "")


class Prog:
    def __init__(self):
        self.calls = []

    def advance(self, n=1, size=None):
        self.calls.append(n)


def _args(nprocs=1, pm_bin_dir=None):
    return SimpleNamespace(nprocs=nprocs, pm_bin_dir=pm_bin_dir or PM_BIN or None)


def _run_folder(tmp_path, files):
    """files: {(uniprot, file name): source path}. Writes 04_pockets/<acc>/<name> and a
    DrugScore.csv (every pocket drugScore 0.9) plus a Protein_info.csv. Returns layout."""
    layout = RunLayout(tmp_path)
    rows, accs = [], []
    for (acc, name), src in files.items():
        dst = layout.pockets / acc / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        rows.append([acc, f"G{acc}_HUMAN", Path(name).stem.split(".")[0].upper(), "P_1", 0.9,
                     name])
        accs.append(acc)
    pd.DataFrame(rows, columns=["UniProt ID", "Gene Name", "Protein", "Pocket", "drugScore",
                                "Pocket Filename"]).to_csv(layout.drugscore_csv)
    layout.structures.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"UniProt ID": sorted(set(accs)),
                  "Gene Name": [f"G{a}_HUMAN" for a in sorted(set(accs))]}
                 ).to_csv(layout.protein_info, index=False)
    return layout


# --- name guard (Review Focus 4) -----------------------------------------------------

def test_name_guard_rejects_long_and_colliding(tmp_path):
    ok = [tmp_path / "P1" / ("x" * 28 + ".pdb"), tmp_path / "P2" / "AF-P06493-F1.pdb_res_P_1.pdb"]
    sim.check_pocket_names(ok)   # exactly 32 characters is allowed

    long_name = "P53350_4o56_A.pdb_res_P_12_10.pdb"   # 33 characters
    with pytest.raises(ValueError, match=r"longer than 32.*P53350_4o56_A\.pdb_res_P_12_10\.pdb"):
        sim.check_pocket_names(ok + [tmp_path / "P53350" / long_name])

    dup = "P53350_4o56_A.pdb_res_P_1.pdb"
    with pytest.raises(ValueError, match=r"collide.*P53350_4o56_A\.pdb_res_P_1\.pdb"):
        sim.check_pocket_names([tmp_path / "A" / dup, tmp_path / "B" / dup])

    # names PocketMatch would pad identically (trailing '_' is padding) also collide
    with pytest.raises(ValueError, match="collide"):
        sim.check_pocket_names([tmp_path / "A" / "p.pdb", tmp_path / "B" / "p.pdb_"])


def test_name_guard_fails_before_pocketmatch_runs(tmp_path):
    src = FIX / "pockets_net2_P53350" / "P53350" / "P53350_4o56_A.pdb_res_P_1.pdb"
    layout = _run_folder(tmp_path, {("P53350", "P53350_4o56_A.pdb_res_P_1.pdb"): src,
                                    ("Q99999", "Q99999_1abc_A.pdb_res_P_12_10_1.pdb"): src})
    (layout.pockets / "P53350" / "P53350_4o56_A.pdb_res_P_1.pdb").unlink()   # also missing
    # a bogus binary dir: the guard must fire before anything looks for PocketMatch
    with pytest.raises(ValueError, match="longer than 32"):
        sim.run(layout, _args(2, pm_bin_dir=tmp_path / "no_pocketmatch"), Prog())
    assert not (layout.similarity / "PocketMatch_score.txt").exists()
    assert not (layout.similarity / "manifest.json").exists()
    assert (layout.similarity / "failed.tsv").read_text().splitlines() == [
        "item\tstage\treason",
        "P53350/P53350_4o56_A.pdb_res_P_1.pdb\tsimilarity\tpocket file missing",
        "Q99999/Q99999_1abc_A.pdb_res_P_12_10_1.pdb\tsimilarity\t"
        "pocket name longer than 32 characters"]


def test_name_guard_in_units_and_dry_run(tmp_path):
    src = FIX / "pockets_net2_P53350" / "P53350" / "P53350_4o56_A.pdb_res_P_1.pdb"
    dup = "P53350_4o56_A.pdb_res_P_1.pdb"
    layout = _run_folder(tmp_path, {("P53350", dup): src, ("Q99999", dup): src})
    with pytest.raises(ValueError, match="collide in PocketMatch") as dry:
        cli.main(["similarity", "--run-dir", str(tmp_path), "--dry-run", "--no-notify"])
    assert not layout.similarity.exists()   # --dry-run writes nothing
    with pytest.raises(ValueError) as real:   # a real run fails in units(), before run()
        cli.main(["similarity", "--run-dir", str(tmp_path), "--no-notify"])
    assert str(real.value) == str(dry.value)
    assert (layout.similarity / "failed.tsv").read_text().splitlines() == [
        "item\tstage\treason",
        f"P53350/{dup}\tsimilarity\tpocket name collides in PocketMatch with Q99999/{dup}",
        f"Q99999/{dup}\tsimilarity\tpocket name collides in PocketMatch with P53350/{dup}"]


# --- inputs --------------------------------------------------------------------------

def test_select_pockets_cutoff_order_and_exact_names(tmp_path):
    src = FIX / "pockets_net2_P53350" / "P53350" / "P53350_4o56_A.pdb_res_P_1.pdb"
    layout = RunLayout(tmp_path)
    names = [("P53350", "P53350_4o56_A.pdb_res_P_2.pdb", 0.7),
             ("O14746", "AF-O14746-F1.pdb_res_P_1.pdb", 0.5),
             ("P53350", "P53350_2ogq_A.pdb_res_P_1.pdb", 0.49),   # below the cutoff
             ("P53350", "P53350_2ogq_A.pdb_res_P_3.pdb", 0.9)]
    for acc, name, _ in names:
        (layout.pockets / acc).mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, layout.pockets / acc / name)
    pd.DataFrame([[a, "G", "X", "P_1", s, n] for a, n, s in names],
                 columns=["UniProt ID", "Gene Name", "Protein", "Pocket", "drugScore",
                          "Pocket Filename"]).to_csv(layout.drugscore_csv)
    got = sim.select_pockets(layout.pockets)
    assert [p.relative_to(layout.pockets).as_posix() for p in got] == [
        "O14746/AF-O14746-F1.pdb_res_P_1.pdb",
        "P53350/P53350_2ogq_A.pdb_res_P_3.pdb",
        "P53350/P53350_4o56_A.pdb_res_P_2.pdb"]

    # exact on-disk name (Review Focus 2): a case mismatch is missing on any filesystem
    (layout.pockets / "P53350" / "P53350_4o56_A.pdb_res_P_2.pdb").rename(
        layout.pockets / "P53350" / "P53350_4o56_a.pdb_res_P_2.pdb")
    table, missing = sim.pocket_table(layout.pockets)
    assert missing == ["P53350/P53350_4o56_A.pdb_res_P_2.pdb"]
    assert list(table["Pocket Filename"]) == ["AF-O14746-F1.pdb_res_P_1.pdb",
                                              "P53350_2ogq_A.pdb_res_P_3.pdb"]
    assert len(sim.select_pockets(layout.pockets)) == 2


def test_missing_drugscore_is_clear(tmp_path):
    with pytest.raises(FileNotFoundError, match="polyphar pockets"):
        sim.units(RunLayout(tmp_path), _args())


# --- scores --------------------------------------------------------------------------

EXPECTED_POCKETSIM = """\
,UniProt ID,Gene Name,Protein,Pocket,Pocket_Path,UniProt ID_2,Gene Name_2,Protein_2,Pocket_2,Pocket_Path_2,Pmin,Pmax
0,P06493,CDK1_HUMAN,AF-P06493-F1,P_1,04_pockets/P06493/AF-P06493-F1.pdb_res_P_1.pdb,P53350,PLK1_HUMAN,4O56_A,P_1,04_pockets/P53350/P53350_4o56_A.pdb_res_P_1.pdb,0.853313,0.805152
3,P48431,SOX2_HUMAN,6T90_A,P_3,04_pockets/P48431/P48431_6t90_A.pdb_res_P_3.pdb,P53350,PLK1_HUMAN,3P36_A,P_1_1,04_pockets/P53350/P53350_3p36_A.pdb_res_P_1_1.pdb,0.643077,0.135714
2,O60218,AK1BA_HUMAN,4GAB_A,P_1,04_pockets/O60218/O60218_4gab_A.pdb_res_P_1.pdb,P49450,CENPA_HUMAN,6MUP_A,P_5,04_pockets/P49450/P49450_6mup_A.pdb_res_P_5.pdb,0.943333,0.028673
1,P06493,CDK1_HUMAN,AF-P06493-F1,P_1,04_pockets/P06493/AF-P06493-F1.pdb_res_P_1.pdb,P53350,PLK1_HUMAN,4O56_A,P_5,04_pockets/P53350/P53350_4o56_A.pdb_res_P_5.pdb,0.714286,0.005018
"""


def test_parse_scores_schema_and_filters(tmp_path):
    """8 net2 lines (names in on-disk chain case): kept AF-PDB and PDB-PDB pairs, a
    same-UniProt pair, a Pmin == 1 pair, a Pmin == Pmax == 1 pair and a NULL____ row
    (synthetic: no 2023 output has one; the format is what Step3-PM_serial prints)."""
    pockets = pd.read_csv(FIX / "pm_scores_sample_DrugScore.csv", index_col=0)
    df = sim.parse_scores(FIX / "pm_scores_sample.txt", pockets,
                          FIX / "pm_scores_sample_Protein_info.csv")
    assert list(df.columns) == sim.POCKETSIM_COLUMNS
    out = tmp_path / "Pocketsim.csv"
    df.to_csv(out)
    assert out.read_text() == EXPECTED_POCKETSIM


def test_parse_scores_rejects_unknown_names(tmp_path):
    pockets = pd.read_csv(FIX / "pm_scores_sample_DrugScore.csv", index_col=0)
    with pytest.raises(ValueError, match="P53350_4o56_A.pdb_res_P_5.pdb"):
        sim.parse_scores(FIX / "pm_scores_sample.txt",
                         pockets[pockets["Pocket"] != "P_5"],
                         FIX / "pm_scores_sample_Protein_info.csv")


# --- block plan ----------------------------------------------------------------------

@pytest.mark.parametrize("n,nprocs", [(0, 4), (1, 4), (2, 16), (5, 1), (16, 4), (16, 16),
                                      (17, 3), (100, 16), (1182, 16), (7, 21)])
def test_plan_blocks_cover_all(n, nprocs):
    blocks = sim.plan_blocks(n, nprocs)
    # contiguous, non-empty (when n > 0), covering range(n) in order
    assert blocks[0][0] == 0 and blocks[-1][1] == n
    assert all(a[1] == b[0] for a, b in zip(blocks, blocks[1:]))
    if n:
        assert all(e > s for s, e in blocks)
    k = len(blocks)
    if nprocs == 1 or n < 2:
        assert k == 1
    else:
        want = next(k for k in range(1, 1000) if k * (k - 1) // 2 >= nprocs)
        assert k == min(want, n)
    # every pair a<b comes from exactly one job
    block_of = [i for i, (s, e) in enumerate(blocks) for _ in range(s, e)]
    jobs = sim.plan_jobs(blocks)
    seen = {}
    for job in jobs:
        members = [a for b in job.blocks for a in range(*blocks[b])]
        for x in range(len(members)):
            for y in range(x + 1, len(members)):
                a, b = members[x], members[y]
                if job.keeps(block_of[a], block_of[b]):
                    seen[(a, b)] = seen.get((a, b), 0) + 1
        assert job.pairs(blocks) == sum(
            1 for x in range(len(members)) for y in range(x + 1, len(members))
            if job.keeps(block_of[members[x]], block_of[members[y]]))
    assert seen == {(a, b): 1 for a in range(n) for b in range(a + 1, n)}
    assert sum(j.pairs(blocks) for j in jobs) == n * (n - 1) // 2


def test_job_output_is_validated(tmp_path):
    """Every pair once, in serial order, under the expected names (a bare pocket name
    without '/' makes cabbage-core record garbage, which must not pass silently)."""
    names = ["a.pdb", "b.pdb", "c.pdb"]
    job = sim.Job((0,), (0,))
    pad = lambda n: n.ljust(32, "_")
    good = [f"{pad(x)} {pad(y)} 0.5 0.4 1 2 3 4\n"
            for x, y in [("a.pdb", "b.pdb"), ("a.pdb", "c.pdb"), ("b.pdb", "c.pdb")]]
    f = tmp_path / "s.txt"
    f.write_text("".join(good))
    assert list(sim._job_lines(f, names, [0, 0, 0], job)) == good
    for bad in (good[:2], good + good[:1], [good[1], good[0], good[2]],
                ["\x14\x9d\xf9 " + good[0][33:]] + good[1:]):
        f.write_bytes("".join(bad).encode("latin-1"))
        with pytest.raises(sim.PocketMatchError, match="unexpected|ends early"):
            list(sim._job_lines(f, names, [0, 0, 0], job))


def _fake_pm_bin(tmp_path):
    d = tmp_path / "fake_pm"
    d.mkdir()
    for b in sim.PM_BINARIES:
        (d / b).write_text("#!/bin/sh\nexit 99\n")
        (d / b).chmod(0o755)
    return d


def test_dead_workers_are_not_retried_forever(tmp_path, monkeypatch):
    """Every worker dies: after the first solo probes all die (nothing ever finished), the
    remaining pockets are failed without running instead of one pool each."""
    src = FIX / "pockets_net2_P53350" / "P53350" / "P53350_4o56_A.pdb_res_P_1.pdb"
    layout = _run_folder(tmp_path, {("P53350", f"P53350_4o56_A.pdb_res_P_{i}.pdb"): src
                                    for i in range(1, 9)})
    monkeypatch.setenv(sim._CRASH_ENV, "*")
    t0 = time.monotonic()
    counts = sim.run(layout, _args(2, pm_bin_dir=_fake_pm_bin(tmp_path)), Prog())
    assert time.monotonic() - t0 < 60
    assert counts["cabbage_failed"] == counts["failed"] == 8 and counts["pairs"] == 0
    rows = (layout.similarity / "failed.tsv").read_text().splitlines()[1:]
    assert len(rows) == 8
    assert all(r.split("\t")[2].startswith(sim.REASON_DIED) for r in rows), rows
    probed = [r for r in rows if not r.endswith("not retried, every pool died")]
    assert 1 <= len(probed) <= 5, rows   # solo runs: two pools' worth + one fresh probe
    assert not multiprocessing.active_children()


def _echo_pm_bin(tmp_path):
    d = _fake_pm_bin(tmp_path)
    (d / sim.CABBAGE).write_text('#!/bin/sh\necho "record $1"\n')
    return d


def _script_pm_bin(tmp_path, cabbage_sh):
    d = _fake_pm_bin(tmp_path)
    (d / sim.CABBAGE).write_text("#!/bin/sh\n" + cabbage_sh)
    return d


@pytest.mark.parametrize("script,cached,text", [
    ("exit 1\n", True, "failed (exit 1)"),
    ("kill -SEGV $$\n", True, "crashed with signal 11 (SIGSEGV)"),   # the tool's own crash
    ("kill -KILL $$\n", False, "killed by signal 9 (SIGKILL, exit -9)"),   # OOM killer
    ("kill -TERM $$\n", False, "killed by signal 15 (SIGTERM, exit -15)"),
    ("exit 137\n", False, "killed by signal 9 (SIGKILL, exit 137)"),
    ("exit 130\n", False, "killed by signal 2 (SIGINT, exit 130)"),
])
def test_only_the_tools_own_verdict_is_cached(tmp_path, script, cached, text):
    pm = _script_pm_bin(tmp_path, script)
    pocket = tmp_path / "P1" / "a.pdb"
    pocket.parent.mkdir()
    pocket.write_text("x\n")
    rec, fl = tmp_path / "a.cabbage", tmp_path / "a.failed"
    reason = sim._cabbage(pocket, pm, rec, fl)
    assert text in reason
    assert fl.exists() == cached and not rec.exists()
    assert ("not cached, retried next run" in reason) == (not cached)
    assert not list(tmp_path.glob("*.tmp"))


def test_failure_while_stopping_is_not_cached(tmp_path, monkeypatch):
    pm = _script_pm_bin(tmp_path, "exit 1\n")
    pocket = tmp_path / "P1" / "a.pdb"
    pocket.parent.mkdir()
    pocket.write_text("x\n")
    monkeypatch.setattr(pool, "_STOP", True)   # this worker is being stopped
    assert "failed (exit 1)" in sim._cabbage(pocket, pm, tmp_path / "r", tmp_path / "f")
    assert not (tmp_path / "f").exists()


def test_killed_pocket_is_retried_next_run(tmp_path, monkeypatch):
    """A pocket whose cabbage-core is killed from outside is dropped for this run only:
    the stage does not count as finished and the pocket is counted again."""
    src = FIX / "pockets_net2_P53350" / "P53350" / "P53350_4o56_A.pdb_res_P_1.pdb"
    layout = _run_folder(tmp_path, {("P1", "a.pdb"): src, ("P2", "b.pdb"): src})
    flag = tmp_path / "oom"
    flag.touch()
    monkeypatch.setenv("FAKE_OOM_FLAG", str(flag))
    pm = _script_pm_bin(tmp_path, '[ "$1" = ./b.pdb ] && [ -e "$FAKE_OOM_FLAG" ] && kill -KILL $$\n'
                                  'echo "record $1"\n')
    counts = sim.run(layout, _args(1, pm_bin_dir=pm), Prog())
    assert counts["cabbage_failed"] == counts["failed"] == 1 and counts["pockets"] == 1
    rows = (layout.similarity / "failed.tsv").read_text().splitlines()
    assert len(rows) == 2 and rows[1].startswith("P2/b.pdb\tsimilarity\t")
    assert "killed by signal 9 (SIGKILL" in rows[1]
    assert not list((layout.similarity / "cabbage").glob("*.failed"))
    flag.unlink()
    assert sim.units(layout, _args(1, pm_bin_dir=pm)).count == 1   # a-b counted again


def test_broken_pool_isolates_only_repeat_suspects(tmp_path, monkeypatch):
    """A crashing task sends the tasks in flight back to a fresh full-size pool; only
    tasks in flight during a second break run alone."""
    pm = _echo_pm_bin(tmp_path)
    pockets, tasks = [], []
    (tmp_path / "cache").mkdir()
    for i in range(24):
        p = tmp_path / "P1" / f"p{i:02d}.pdb"
        p.parent.mkdir(exist_ok=True)
        p.write_text("x\n")
        rec, fl = tmp_path / "cache" / f"{i}.cabbage", tmp_path / "cache" / f"{i}.failed"
        tasks.append((p.name, sim._cabbage_worker, (str(p), str(pm), str(rec), str(fl))))
        pockets.append(p)
    monkeypatch.setenv(sim._CRASH_ENV, "p05.pdb")
    got = []
    out = sim._run_pool(tasks, 4, lambda item, r: got.append(item))
    assert out.failures == [("p05.pdb", sim.REASON_DIED)]
    assert "p05.pdb" in out.isolated and len(out.isolated) <= 4   # at most one pool's worth
    assert sorted(got) == sorted(p.name for p in pockets if p.name != "p05.pdb")
    assert (tmp_path / "cache" / "0.cabbage").read_text() == "record ./p00.pdb\n"
    assert not multiprocessing.active_children()


# --- PocketMatch (inside the image) --------------------------------------------------

def _sample_layout(tmp_path):
    pockets = sorted((Path(PM_BIN) / "Sample_pockets").glob("*.pdb"))
    files = {}
    for i, p in enumerate(pockets):   # spread over 3 accessions, names kept
        files[(f"P{i % 3:05d}", p.name)] = p
    return _run_folder(tmp_path, files)


@pytest.mark.pocketmatch
def test_run_pocketmatch_reproduces_upstream_sample(tmp_path):
    pockets = sorted((Path(PM_BIN) / "Sample_pockets").glob("*.pdb"))
    out = sim.run_pocketmatch(pockets, PM_BIN, tmp_path)
    ref = Path(PM_BIN) / "sample_outputs" / "PocketMatch_score.txt"
    assert out.read_bytes() == ref.read_bytes()


@pytest.mark.pocketmatch
def test_parallel_equals_serial_on_sample(tmp_path):
    outs = {}
    for nprocs in (1, 4):
        layout = _sample_layout(tmp_path / f"n{nprocs}")
        counts = sim.run(layout, _args(nprocs), Prog())
        assert counts["failed"] == 0 and counts["pairs"] == 120
        outs[nprocs] = layout.similarity
    for name in ("PocketMatch_score.txt", "Pocketsim.csv"):
        assert (outs[4] / name).read_bytes() == (outs[1] / name).read_bytes(), name
    assert len((outs[1] / "PocketMatch_score.txt").read_text().splitlines()) == 120
    m = json.loads((outs[4] / "manifest.json").read_text())
    assert m["stage"] == "similarity" and m["args"]["nprocs"] == 4
    assert m["versions"]["pocketmatch_commit"]
    assert set(m["versions"]["pocketmatch_sha256"]) == set(sim.PM_BINARIES)
    assert not list(outs[4].rglob("*.tmp"))
    assert not (outs[4] / "jobs").exists()
    assert len(list((outs[4] / "cabbage").glob("*.cabbage"))) == 16


@pytest.mark.pocketmatch
def test_finished_stage_is_not_recomputed(tmp_path, capsys):
    """A second run on the same inputs reports 0 units, runs no PocketMatch and leaves
    every file as it was; any nprocs. Changed inputs invalidate it."""
    layout = _sample_layout(tmp_path)
    first = sim.run(layout, _args(4), Prog())
    files = {p: p.stat().st_mtime_ns for p in layout.similarity.rglob("*")}
    assert sim.units(layout, _args(4)).count == 0
    assert sim.units(layout, _args(1)).count == 0   # independent of nprocs

    def boom(*a, **kw):
        raise AssertionError("PocketMatch must not run")
    with unittest.mock.patch.object(sim, "_similarity", boom), \
            unittest.mock.patch.object(sim, "run_pocketmatch", boom), \
            unittest.mock.patch.object(sim, "cabbage_record", boom):
        assert sim.run(layout, _args(1), Prog()) == first
        capsys.readouterr()
        assert cli.main(["similarity", "--run-dir", str(tmp_path), "--no-notify",
                         "--nprocs", "4"]) == 0
    assert "similarity: 0 units (pair)" in capsys.readouterr().out
    assert {p: p.stat().st_mtime_ns for p in layout.similarity.rglob("*")} == files

    # a truncated output, changed Protein_info, or a changed pocket invalidates it
    score = layout.similarity / "PocketMatch_score.txt"
    data = score.read_bytes()
    score.write_bytes(data[:-1])
    assert sim.units(layout, _args(4)).count == 120
    score.write_bytes(data)
    assert sim.units(layout, _args(4)).count == 0
    info = layout.protein_info.read_text()
    layout.protein_info.write_text(info + "\n")
    assert sim.units(layout, _args(4)).count == 120
    layout.protein_info.write_text(info)
    assert sim.units(layout, _args(4)).count == 0
    victim = sim.select_pockets(layout.pockets)[3]
    victim.write_text(victim.read_text() + "REMARK changed\n")
    assert sim.units(layout, _args(4)).count == 120
    again = sim.run(layout, _args(4), Prog())
    assert again["skipped_existing"] == 0 and again["pairs"] == 120
    assert sim.units(layout, _args(4)).count == 0


@pytest.mark.pocketmatch
def test_resume_skips_done_jobs_and_counts_remaining(tmp_path):
    layout = _sample_layout(tmp_path)
    args = _args(4)

    class Boom(Prog):
        def advance(self, n=1, size=None):
            super().advance(n)
            raise KeyboardInterrupt

    assert sim.units(layout, args).count == 120
    with pytest.raises(KeyboardInterrupt):
        sim.run(layout, args, Boom())
    deadline = time.monotonic() + 60
    while multiprocessing.active_children() and time.monotonic() < deadline:
        time.sleep(0.1)
    assert not multiprocessing.active_children(), "worker processes still alive after 60 s"
    assert not (layout.similarity / "manifest.json").exists()
    assert not (layout.similarity / "PocketMatch_score.txt").exists()
    jobs = layout.similarity / "jobs"
    done = sorted(jobs.glob("*.done"))
    assert done   # at least the job that reported back
    remaining = sim.units(layout, args).count
    assert remaining < 120
    assert remaining == 120 - sum(json.loads(p.read_text())["pairs"] for p in done)

    before = {p.with_suffix(".txt"): p.with_suffix(".txt").stat().st_mtime_ns for p in done}
    rewritten = []
    real_job = sim._run_jobs

    def spy(plan, *a, **kw):   # finished jobs are not rerun
        rewritten.extend(j.tag for j in plan.jobs if j.tag not in plan.done)
        assert {t.stem for t in before} == set(plan.done)
        return real_job(plan, *a, **kw)
    prog = Prog()
    with unittest.mock.patch.object(sim, "_run_jobs", spy):
        counts = sim.run(layout, args, prog)
    assert not {t.stem for t in before} & set(rewritten)
    assert sum(prog.calls) == remaining
    assert counts["skipped_existing"] == 120 - remaining and counts["failed"] == 0
    assert not jobs.exists()   # deleted once Pocketsim.csv is written

    serial = _sample_layout(tmp_path / "serial")
    sim.run(serial, _args(1), Prog())
    assert ((layout.similarity / "PocketMatch_score.txt").read_bytes()
            == (serial.similarity / "PocketMatch_score.txt").read_bytes())


def _interrupted_run(layout, args):
    class Boom(Prog):
        def advance(self, n=1, size=None):
            if len(self.calls) == 2:
                raise KeyboardInterrupt
            super().advance(n)
    with pytest.raises(KeyboardInterrupt):
        sim.run(layout, args, Boom())
    deadline = time.monotonic() + 60
    while multiprocessing.active_children() and time.monotonic() < deadline:
        time.sleep(0.1)


@pytest.mark.pocketmatch
def test_done_markers_are_checked(tmp_path):
    """A .done marker counts only with its fingerprint (pocket contents, PocketMatch
    binaries) and the size of its .txt as written."""
    layout = _sample_layout(tmp_path)
    args = _args(4)
    _interrupted_run(layout, args)
    jobs = layout.similarity / "jobs"
    done = sorted(jobs.glob("*.done"))
    assert len(done) >= 2
    remaining = sim.units(layout, args).count
    assert remaining < 120

    # a host crash can leave a short .txt behind a marker: that job is redone
    txt = done[0].with_suffix(".txt")
    txt.write_bytes(txt.read_bytes()[:-10])
    pairs0 = json.loads(done[0].read_text())["pairs"]
    assert sim.units(layout, args).count == remaining + pairs0

    # other PocketMatch binaries (an image update) invalidate every job
    paths = sim.select_pockets(layout.pockets)
    digests = [sim.sha256_file(p) for p in paths]
    bins = sim._bin_hashes(Path(PM_BIN))
    assert sim._plan(paths, 4, layout.similarity, digests, bins).done
    other = {**bins, sim.PM_SERIAL: "0" * 64}
    assert sim._plan(paths, 4, layout.similarity, digests, other).done == {}

    # changed pocket contents invalidate every job holding them
    for p in paths:
        p.write_text(p.read_text() + "REMARK changed\n")
    assert sim.units(layout, args).count == 120


@pytest.mark.pocketmatch
def test_rejected_and_missing_pockets_are_dropped(tmp_path):
    """A pocket cabbage-core rejects and a listed pocket missing on disk are recorded in
    failed.tsv and left out; the stage completes. A rerun takes cabbage-core's verdicts
    from the cache instead of running it again."""
    s = Path(PM_BIN) / "Sample_pockets"
    layout = _run_folder(tmp_path, {("P1", "a.pdb"): s / "1g64.pdb", ("P2", "b.pdb"): s / "1g64.pdb",
                                    ("P3", "c.pdb"): s / "1i9c.pdb", ("P4", "d.pdb"): s / "1iwb.pdb"})
    one_atom = next(l for l in (s / "1g64.pdb").read_text().splitlines() if l.startswith("ATOM"))
    (layout.pockets / "P2" / "b.pdb").write_text(one_atom + "\n")   # cabbage-core aborts
    (layout.pockets / "P4" / "d.pdb").unlink()
    assert sim.units(layout, _args(1)).count == 3   # a, b, c: b's fate is not known yet

    counts = sim.run(layout, _args(1), Prog())
    assert counts["pockets"] == 2 and counts["pairs"] == 1 and counts["rows"] == 1
    assert counts["missing"] == 1 and counts["cabbage_failed"] == 1 and counts["failed"] == 2
    rows = (layout.similarity / "failed.tsv").read_text().splitlines()
    assert rows[:2] == ["item\tstage\treason", "P4/d.pdb\tsimilarity\tpocket file missing"]
    assert rows[2].startswith("P2/b.pdb\tsimilarity\tStep0-cabbage-core ./b.pdb crashed with "
                              "signal 6 (SIGABRT)")
    assert "malloc" in rows[2] and len(rows) == 3   # the tool's stderr tail
    assert (layout.similarity / "manifest.json").is_file()
    assert [l.split()[:2] for l in (layout.similarity / "PocketMatch_score.txt").read_text()
            .splitlines()] == [["a.pdb" + "_" * 27, "c.pdb" + "_" * 27]]

    cache = layout.similarity / "cabbage"
    entries = {p.name: p.stat().st_mtime_ns for p in cache.iterdir()}
    assert sorted(n.rsplit(".", 1)[1] for n in entries if n != sim.CABBAGE_STAMP) == [
        "cabbage", "cabbage", "failed"]
    assert sim.units(layout, _args(1)).count == 0   # finished on these inputs
    (layout.similarity / "Pocketsim.csv").unlink()   # force a recompute
    assert sim.units(layout, _args(1)).count == 1   # b is known to fail now
    again = sim.run(layout, _args(1), Prog())
    assert again == counts
    assert {p.name: p.stat().st_mtime_ns for p in cache.iterdir()} == entries
    assert (layout.similarity / "failed.tsv").read_text().splitlines() == rows


@pytest.mark.pocketmatch
def test_worker_crash_is_attributed(tmp_path, monkeypatch):
    layout = _sample_layout(tmp_path)
    monkeypatch.setenv(sim._CRASH_ENV, "1req.pdb")   # inherited by spawned workers
    counts = sim.run(layout, _args(4), Prog())
    assert counts["failed"] == 1 and counts["pockets"] == 15 and counts["pairs"] == 105
    rows = (layout.similarity / "failed.tsv").read_text().splitlines()
    assert rows == ["item\tstage\treason",
                    "P00000/1req.pdb\tsimilarity\t" + sim.REASON_DIED]
    monkeypatch.delenv(sim._CRASH_ENV)   # a death is not cached: the rerun retries it
    assert sim.run(layout, _args(4), Prog())["failed"] == 0


def test_cli_dry_run(tmp_path, capsys):
    src = FIX / "pockets_net2_P53350" / "P53350" / "P53350_4o56_A.pdb_res_P_1.pdb"
    layout = _run_folder(tmp_path, {("P53350", f"P53350_4o56_A.pdb_res_P_{i}.pdb"): src
                                    for i in range(1, 6)})
    assert cli.main(["similarity", "--run-dir", str(tmp_path), "--dry-run", "--no-notify",
                     "--nprocs", "4"]) == 0
    assert "similarity: 10 units (pair)" in capsys.readouterr().out
    assert not layout.similarity.exists()


# --- regression: Liver_cancer/net2 (2023) --------------------------------------------

def _net2():
    return Path(os.environ["ADEP_DATA"]) / "disease_subtypes" / "Liver_cancer" / "net2"


@pytest.fixture(scope="module")
def net2_run(tmp_path_factory):
    """A run folder holding only the net2 pockets with drugScore >= 0.5, copied under
    their ON-DISK names, and a DrugScore.csv whose Pocket Filename matches those names
    (the legacy file lowercases the chain; it only resolved on case-insensitive NTFS).
    Legacy AlphaFold `Protein` AF-<acc> is written in the current form AF-<acc>-F1.
    9 rows (O95235) name files that are not on disk (the legacy script skipped missing
    files silently); they stay in DrugScore.csv for the stage to report and drop."""
    net2 = _net2()
    base = tmp_path_factory.mktemp("net2")
    layout = RunLayout(base)
    ds = pd.read_csv(net2 / "pockets" / "DrugScore.csv", index_col=0, dtype={"Pocket": str})
    ds = ds[ds["drugScore"] >= 0.5].copy()
    listing = {}
    names = []
    ds["legacy"] = ds["Pocket Filename"]
    for acc, fname in zip(ds["UniProt ID"], ds["Pocket Filename"]):
        d = net2 / "pockets" / acc
        if acc not in listing:
            listing[acc] = {n.lower(): n for n in os.listdir(d)}
        real = listing[acc].get(fname.lower())
        if real is not None:
            (layout.pockets / acc).mkdir(parents=True, exist_ok=True)
            shutil.copyfile(d / real, layout.pockets / acc / real)
        names.append(real)
    ds["Pocket Filename"] = names
    absent = ds["Pocket Filename"].isna()
    ds.loc[absent, "Pocket Filename"] = ds.loc[absent, "legacy"]   # left for the stage to report
    ds = ds.drop(columns="legacy")
    af = ds["Protein"].str.match(r"^AF-[^-]+$")
    ds.loc[af, "Protein"] = ds.loc[af, "Protein"] + "-F1"
    ds.to_csv(layout.drugscore_csv)
    layout.structures.mkdir(parents=True)
    shutil.copyfile(net2 / "structures" / "Protein_info.csv", layout.protein_info)
    return layout


def _score_map(path):
    """unordered pair (lowercased names, trailing '_' stripped) -> (Pmin, Pmax)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            a, b = f[0].rstrip("_").lower(), f[1].rstrip("_").lower()
            key = (a, b) if a <= b else (b, a)
            assert key not in out, f"duplicate pair {key}"
            out[key] = (f[2], f[3])
    return out


KEY = ["UniProt ID", "Protein", "Pocket", "UniProt ID_2", "Protein_2", "Pocket_2"]


def _pocketsim_keys(path):
    """Multiset of unordered rows (UniProt ID, Protein, Pocket) x2 + (Pmin, Pmax), Protein
    upper-cased. The pair orientation follows the pocket order, which is now sorted by
    (UniProt ID, Pocket Filename) instead of the legacy directory listing."""
    from collections import Counter

    df = pd.read_csv(path, index_col=0, dtype=str)
    c = Counter()
    for r in df[KEY + ["Pmin", "Pmax"]].itertuples(index=False):
        a = (r[0], r[1].upper(), r[2])
        b = (r[3], r[4].upper(), r[5])
        c[(min(a, b), max(a, b), float(r[6]), float(r[7]))] += 1
    return c


@pytest.mark.pocketmatch
@pytest.mark.regression
@pytest.mark.parametrize("nprocs", [1, 16])
def test_reproduces_net2_2023(net2_run, tmp_path, nprocs):
    run_dir = tmp_path / "run"
    shutil.copytree(net2_run.run_dir, run_dir)
    layout = RunLayout(run_dir)
    t0 = time.monotonic()
    counts = sim.run(layout, _args(nprocs), Prog())
    wall = time.monotonic() - t0
    print(f"\nnet2 similarity nprocs={nprocs}: {counts['pockets']} pockets, "
          f"{counts['pairs']} pairs, wall {wall:.1f} s")

    rows = (layout.similarity / "failed.tsv").read_text().splitlines()[1:]
    assert len(rows) == 9 and counts["missing"] == counts["failed"] == 9
    assert all(r.startswith("O95235/AF-O95235-F1.pdb_res_P_") and
               r.endswith("\tsimilarity\tpocket file missing") for r in rows), rows
    assert counts["pockets"] == 1182 and counts["cabbage_failed"] == 0

    ref_dir = _net2() / "PM_net2"
    got, ref = _score_map(layout.similarity / "PocketMatch_score.txt"), _score_map(
        ref_dir / "PocketMatch_score.txt")
    assert len(got) == len(ref)
    assert got.keys() == ref.keys()
    diff = [k for k in ref if got[k] != ref[k]]
    assert not diff, f"{len(diff)} pairs differ, e.g. {diff[:5]}"

    got_rows = _pocketsim_keys(layout.similarity / "Pocketsim.csv")
    ref_rows = _pocketsim_keys(ref_dir / "Pocketsim_net2.csv")
    assert sum(got_rows.values()) == sum(ref_rows.values())
    assert got_rows == ref_rows, (
        f"only new: {list((got_rows - ref_rows).items())[:5]}; "
        f"only legacy: {list((ref_rows - got_rows).items())[:5]}")

    # a second invocation on the same folder recognises the finished stage
    files = {p: p.stat().st_mtime_ns for p in layout.similarity.rglob("*")}
    t0 = time.monotonic()
    assert sim.units(layout, _args(nprocs)).count == 0
    assert sim.run(layout, _args(nprocs), Prog()) == counts
    print(f"net2 similarity nprocs={nprocs} second invocation (units + run): "
          f"{time.monotonic() - t0:.1f} s")
    assert {p: p.stat().st_mtime_ns for p in layout.similarity.rglob("*")} == files
