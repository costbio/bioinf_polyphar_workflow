"""Stage 05: pocket similarity with upstream PocketMatch v2.0, serial or block-pair parallel.

Inputs: 04_pockets/DrugScore.csv, rows with drugScore >= 0.5; each pocket file is
04_pockets/<UniProt ID>/<Pocket Filename>, found by its exact on-disk name. The pocket
list is ordered by (UniProt ID, Pocket Filename) so every run sees the same order. A
listed file that is missing is recorded in failed.tsv (`pocket file missing`) and dropped.
PocketMatch records a file's base name truncated to 32 characters, so names longer than
32 or colliding fail the stage before PocketMatch runs (also in --dry-run), listed in
failed.tsv.

PocketMatch (the upstream prebuilt 32-bit binaries in $PM_BIN_DIR) is run as the upstream
pipeline does: `Step0-cabbage-core ./<pocket>` (from the pocket's directory) for each
pocket into outfile.cabbage, then `Step0-END-FILE`, then `Step3-PM_serial outfile.cabbage`,
every exit code checked, in private temp dirs. Its PocketMatch_pairs.txt (hundreds of
bytes per pair, not used) goes to /dev/null.

Phase 1, cabbage: cabbage-core's output depends only on the pocket file (content and name),
so it runs once per pocket, in parallel, into the persistent cache 05_similarity/cabbage/
<sha256>_<name>.cabbage (atomic writes; the cache is cleared when the cabbage-core binary
changes). A pocket cabbage-core rejects is recorded in failed.tsv with the tool's stderr
tail (and cached as <key>.failed, so a resume does not rerun it) and dropped.

Phase 2, PocketMatch: the surviving list is cut into k contiguous blocks B_0..B_{k-1}, k
the smallest with k(k-1)/2 >= nprocs. Job (i, j), i < j, concatenates the records of
B_i ++ B_j, runs PM_serial and keeps the cross pairs; the within-block pairs of B_i are
kept by exactly one job (the ring neighbour (i, i+1 mod k)). Each job writes
05_similarity/jobs/<i>-<j>.txt and then a <i>-<j>.done marker holding the .txt size and
a fingerprint (its pockets' names and sha256, the PocketMatch binaries' sha256); a rerun
skips jobs whose marker matches. A failed job is fatal. The job files are merged in the
order a serial run emits pairs ((index(name1), index(name2)) in the ordered list), so the
merged PocketMatch_score.txt is byte-identical for every nprocs. jobs/ is deleted once
PocketMatch_score.txt and Pocketsim.csv are written.

Outputs: 05_similarity/PocketMatch_score.txt, Pocketsim.csv (legacy columns and filters;
Pocket_Path relative to the run folder), failed.tsv, manifest.json (written last).
"""

from __future__ import annotations

import hashlib
import heapq
import json
import os
import signal
import subprocess
import shutil
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

from polyphar.common import PM_NAME_LIMIT, WorkUnits, append_failure, sha256_file, write_manifest

STAGE = "similarity"
DRUGSCORE_CUTOFF = 0.5
POCKETS_REL = "04_pockets"   # Pocket_Path prefix: relative to the run folder (Ruling 3)
CABBAGE, END_FILE, PM_SERIAL = "Step0-cabbage-core", "Step0-END-FILE", "Step3-PM_serial"
PM_BINARIES = (CABBAGE, END_FILE, PM_SERIAL)
SCORE_FILE = "PocketMatch_score.txt"
JOBS = "jobs"
CABBAGE_DIR = "cabbage"
CABBAGE_STAMP = "BINARY_SHA256"
NULL = "NULL____"
POCKETSIM_COLUMNS = ["UniProt ID", "Gene Name", "Protein", "Pocket", "Pocket_Path",
                     "UniProt ID_2", "Gene Name_2", "Protein_2", "Pocket_2", "Pocket_Path_2",
                     "Pmin", "Pmax"]
REASON_DIED = "worker process died (BrokenProcessPool)"
REASON_MISSING = "pocket file missing"
SOLO_RETRY_PROBES = 3
TERMINATE_GRACE_S = 15.0


def add_arguments(parser) -> None:
    parser.add_argument("--pm-bin-dir", default=None,
                        help="directory holding the PocketMatch binaries (default: $PM_BIN_DIR)")


class PocketMatchError(RuntimeError):
    """A PocketMatch step failed or wrote unexpected output (a plain message, so it
    always pickles across the process boundary)."""


class JobsFailed(RuntimeError):
    def __init__(self, failures: list[tuple[str, str]]):
        self.failures = failures
        first = "; ".join(f"{item}: {reason}" for item, reason in failures[:3])
        super().__init__(f"{len(failures)} PocketMatch job(s) failed: {first}")


# --- inputs --------------------------------------------------------------------------

def pocket_table(pockets_dir, cutoff: float = DRUGSCORE_CUTOFF):
    """(table, missing): the DrugScore.csv rows with drugScore >= cutoff whose pocket file
    exists, one per file, sorted by (UniProt ID, Pocket Filename); and the listed files that
    do not exist, as `<UniProt ID>/<Pocket Filename>` (compared by exact name, so a case
    mismatch counts as missing on any filesystem)."""
    import pandas as pd

    pockets_dir = Path(pockets_dir)
    path = pockets_dir / "DrugScore.csv"
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; run `polyphar pockets` first")
    text_cols = ["UniProt ID", "Gene Name", "Protein", "Pocket", "Pocket Filename"]
    df = pd.read_csv(path, index_col=0, dtype={c: str for c in text_cols})
    df[text_cols] = df[text_cols].fillna("")
    df = df[pd.to_numeric(df["drugScore"]) >= cutoff]
    df = df.drop_duplicates(["UniProt ID", "Pocket Filename"])
    df = df.sort_values(["UniProt ID", "Pocket Filename"], kind="mergesort")
    df = df.reset_index(drop=True)
    listing: dict[str, set[str]] = {}
    present, missing = [], []
    for acc, name in zip(df["UniProt ID"], df["Pocket Filename"]):
        if acc not in listing:
            d = pockets_dir / acc
            listing[acc] = set(os.listdir(d)) if d.is_dir() else set()
        present.append(name in listing[acc])
        if not present[-1]:
            missing.append(f"{acc}/{name}")
    return df[present].reset_index(drop=True), missing


def _paths(pockets_dir, table) -> list[Path]:
    return [Path(pockets_dir) / a / f for a, f in zip(table["UniProt ID"], table["Pocket Filename"])]


def select_pockets(pockets_dir, cutoff: float = DRUGSCORE_CUTOFF) -> list[Path]:
    """Existing pocket files with drugScore >= cutoff, in the stage's deterministic order."""
    return _paths(pockets_dir, pocket_table(pockets_dir, cutoff)[0])


def _name_problems(pockets) -> tuple[list[tuple[str, str]], str]:
    """([(item, reason)], summary) for names PocketMatch cannot keep apart: longer than
    32 characters (it truncates them), colliding, or with whitespace / a trailing '_'
    (it pads names with '_')."""
    pockets = [Path(p) for p in pockets]
    rows: list[tuple[str, str]] = []
    long = [p for p in pockets if len(p.name) > PM_NAME_LIMIT]
    rows += [(_item(p), f"pocket name longer than {PM_NAME_LIMIT} characters") for p in long]
    seen: dict[str, list[Path]] = {}
    for p in pockets:
        seen.setdefault(p.name[:PM_NAME_LIMIT].rstrip("_"), []).append(p)
    collide = {k: v for k, v in seen.items() if len(v) > 1}
    for v in collide.values():
        rows += [(_item(p), "pocket name collides in PocketMatch with "
                  + ", ".join(_item(q) for q in v if q is not p)) for p in v]
    bad = [p for p in pockets if p.name != p.name.strip() or any(c.isspace() for c in p.name)
           or p.name.endswith("_")]
    rows += [(_item(p), "pocket name with whitespace or a trailing '_'") for p in bad]
    problems = []
    if long:
        names = sorted({p.name for p in long})
        problems.append(f"{len(names)} pocket name(s) longer than {PM_NAME_LIMIT} characters "
                        f"(PocketMatch truncates them): " + ", ".join(names[:20]))
    if collide:
        problems.append(f"{len(collide)} pocket name(s) collide in PocketMatch: " + "; ".join(
            f"{k}: {', '.join(map(str, v))}" for k, v in list(collide.items())[:20]))
    if bad:
        problems.append("pocket name(s) with whitespace or a trailing '_': "
                        + ", ".join(sorted({p.name for p in bad})[:20]))
    return rows, "; ".join(problems)


def check_pocket_names(pockets) -> None:
    """Raise ValueError, listing the names, when a name is longer than 32 characters,
    cannot be read back (whitespace, trailing '_'), or two names collide in PocketMatch."""
    _, summary = _name_problems(pockets)
    if summary:
        raise ValueError(summary)


def _item(path) -> str:
    path = Path(path)
    return f"{path.parent.name}/{path.name}"


# --- PocketMatch ---------------------------------------------------------------------

def _pm_bin_dir(args) -> Path:
    d = getattr(args, "pm_bin_dir", None) or os.environ.get("PM_BIN_DIR")
    if not d:
        raise FileNotFoundError("PocketMatch binaries not configured: set PM_BIN_DIR or pass "
                                "--pm-bin-dir (the polyphar image has /opt/pocketmatch)")
    d = Path(d)
    missing = [b for b in PM_BINARIES if not os.access(d / b, os.X_OK)]
    if missing:
        raise FileNotFoundError(f"PocketMatch binaries missing or not executable in {d}: "
                                + ", ".join(missing))
    return d


def _check(proc, what: str) -> None:
    if proc.returncode != 0:
        tail = (proc.stderr or b"")[-500:].decode(errors="replace").strip()
        raise PocketMatchError(f"{what} failed (exit {proc.returncode})"
                               + (f": {tail}" if tail else ""))


def cabbage_record(pocket, pm_bin_dir, out) -> None:
    """`Step0-cabbage-core ./<base name>` run from the pocket's directory, into `out`.
    cabbage-core records the text after the last '/' as the name (a bare name without
    '/' records garbage); `./` keeps the argument short however deep the run folder is."""
    p = Path(pocket)
    with open(out, "wb") as fh:
        proc = subprocess.run([str(Path(pm_bin_dir) / CABBAGE), f"./{p.name}"], cwd=p.parent,
                              stdout=fh, stderr=subprocess.PIPE)
    _check(proc, f"{CABBAGE} ./{p.name}")


def run_pocketmatch(pockets, pm_bin_dir, work_dir, records=None) -> Path:
    """Upstream PocketMatch all-vs-all on `pockets` (in that order) inside work_dir:
    Step0-cabbage-core per pocket into outfile.cabbage, Step0-END-FILE, Step3-PM_serial.
    `records`, if given, are the pockets' cabbage-core outputs computed earlier (by
    cabbage_record), concatenated instead of rerunning cabbage-core. Returns
    work_dir/PocketMatch_score.txt."""
    pm, work = Path(pm_bin_dir), Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    cabbage = work / "outfile.cabbage"
    with open(cabbage, "wb") as fh:
        for i, p in enumerate(pockets):
            rec = Path(records[i]) if records is not None else work / "pocket.cabbage"
            if records is None:
                cabbage_record(p, pm, rec)
            with open(rec, "rb") as src:
                shutil.copyfileobj(src, fh)
        fh.flush()
        _check(subprocess.run([str(pm / END_FILE)], cwd=work, stdout=fh,
                              stderr=subprocess.PIPE), END_FILE)
    pairs = work / "PocketMatch_pairs.txt"
    pairs.unlink(missing_ok=True)
    os.symlink(os.devnull, pairs)
    score = work / SCORE_FILE
    score.unlink(missing_ok=True)
    proc = subprocess.run([str(pm / PM_SERIAL), cabbage.name], cwd=work,
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    _check(proc, PM_SERIAL)
    if not score.is_file():
        raise PocketMatchError(f"{PM_SERIAL} wrote no {SCORE_FILE}")
    return score


# --- block-pair plan -----------------------------------------------------------------

def plan_blocks(n: int, nprocs: int) -> list[tuple[int, int]]:
    """Contiguous blocks [(start, end)) over range(n): k = 1 for nprocs = 1 (or n < 2),
    else the smallest k with k(k-1)/2 >= nprocs, at most n (no empty block)."""
    if nprocs <= 1 or n < 2:
        return [(0, n)]
    k = 2
    while k * (k - 1) // 2 < nprocs:
        k += 1
    k = min(k, n)
    return [(i * n // k, (i + 1) * n // k) for i in range(k)]


@dataclass(frozen=True)
class Job:
    blocks: tuple[int, ...]   # (i,) for a single block, else (i, j) with i < j
    within: tuple[int, ...]   # blocks whose within-block pairs this job keeps

    @property
    def tag(self) -> str:
        i, j = self.blocks[0], self.blocks[-1]
        return f"{i:03d}-{j:03d}"

    def keeps(self, block_a: int, block_b: int) -> bool:
        return block_a != block_b or block_a in self.within

    def pairs(self, blocks) -> int:
        sizes = [blocks[b][1] - blocks[b][0] for b in self.blocks]
        cross = sizes[0] * sizes[1] if len(sizes) == 2 else 0
        return cross + sum((blocks[b][1] - blocks[b][0]) * (blocks[b][1] - blocks[b][0] - 1) // 2
                           for b in self.within)


def plan_jobs(blocks) -> list[Job]:
    """One job per block pair i < j; within-block pairs of B_i belong to the job
    (i, i+1 mod k), so every pair is computed by exactly one job."""
    k = len(blocks)
    if k == 1:
        return [Job((0,), (0,))]
    owner = {}
    for i in range(k):
        owner.setdefault(tuple(sorted((i, (i + 1) % k))), []).append(i)
    return [Job((i, j), tuple(owner.get((i, j), ()))) for i in range(k) for j in range(i + 1, k)]


@dataclass(frozen=True)
class _Plan:
    paths: list[Path]
    blocks: list[tuple[int, int]]
    jobs: list[Job]
    fingerprints: dict[str, str]   # job tag -> fingerprint
    done: dict[str, int]           # job tag -> kept pairs, for jobs already finished

    @property
    def total(self) -> int:
        n = len(self.paths)
        return n * (n - 1) // 2

    @property
    def remaining(self) -> int:
        return sum(j.pairs(self.blocks) for j in self.jobs if j.tag not in self.done)


def _members(job: Job, blocks) -> list[int]:
    return [a for b in job.blocks for a in range(*blocks[b])]


def _read_done(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def _plan(paths: list[Path], nprocs: int, work_dir: Path, digests: list[str],
          bins: dict | None) -> _Plan:
    """Blocks, jobs and fingerprints for `paths`; a job is done when its .done marker has
    the same fingerprint and the size of its .txt as written. `bins` are the PocketMatch
    binaries' sha256 (None when unknown: then no job counts as done)."""
    blocks = plan_blocks(len(paths), nprocs)
    jobs = plan_jobs(blocks) if len(paths) >= 2 else []
    jobs_dir = Path(work_dir) / JOBS
    fps, done = {}, {}
    for job in jobs:
        payload = {"blocks": [[[paths[a].name, digests[a]] for a in range(*blocks[b])]
                              for b in job.blocks],
                   "within": list(job.within), "pocketmatch_sha256": bins}
        fp = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        fps[job.tag] = fp
        marker = _read_done(jobs_dir / f"{job.tag}.done")
        try:
            size = (jobs_dir / f"{job.tag}.txt").stat().st_size
        except OSError:
            continue
        if bins is not None and marker.get("fingerprint") == fp and marker.get("bytes") == size:
            done[job.tag] = int(marker["pairs"])
    return _Plan(list(paths), blocks, jobs, fps, done)


# --- one block-pair job (runs in a worker process) -----------------------------------

_BUSY = False   # worker process: inside a task (temp files / a PocketMatch child to clean up)
_STOP = False   # worker process: SIGTERM received while busy


def _worker_init() -> None:
    # SIGTERM (the parent stopping the pool): an idle worker exits at once. A busy one
    # unwinds through subprocess.run, which kills the PocketMatch child, and through the
    # temp dir / temp file cleanup, then exits in _task. Further SIGTERMs (the pool's
    # terminate_broken sends one) do not cut that cleanup short.
    def handler(signum, frame):
        global _STOP
        if not _BUSY:
            os._exit(128 + signum)
        if not _STOP:
            _STOP = True
            raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, handler)


def _task(fn):
    """Run a worker task: mark the process busy, carry failures across the process
    boundary as PocketMatchError (always pickles), exit once stopped by SIGTERM."""
    global _BUSY
    _BUSY = True
    try:
        return fn()
    except PocketMatchError:
        raise
    except Exception as exc:
        raise PocketMatchError(f"{type(exc).__name__}: {exc}") from None
    except BaseException:
        os._exit(143)   # SIGTERM / Ctrl-C: the cleanup inside fn has run
    finally:
        _BUSY = False
        if _STOP:
            os._exit(143)   # the parent is stopping the pool; never go back to idle


def _atomic_write_lines(path: Path, lines, encoding: str = "latin-1") -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline="") as fh:
            for line in lines:
                fh.write(line)
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _job_lines(score: Path, names: list[str], block_of: list[int], job: Job):
    """Kept lines of a job's score file, checking that PocketMatch emitted every pair
    of `names` once, in serial order, under the expected names."""
    m = len(names)
    x, y = 0, 1
    with open(score, encoding="latin-1", newline="") as fh:
        for line in fh:
            f = line.split(None, 2)
            if y >= m or len(f) < 3 or f[0].rstrip("_") != names[x] \
                    or f[1].rstrip("_") != names[y]:
                want = f"{names[x]} {names[y]}" if y < m else "end of file"
                raise PocketMatchError(f"unexpected {SCORE_FILE} line {line[:80]!r} "
                                       f"(expected {want})")
            if not line.endswith("\n"):
                line += "\n"
            if job.keeps(block_of[x], block_of[y]):
                yield line
            y += 1
            if y == m:
                x, y = x + 1, x + 2
    if m >= 2 and x < m - 1:
        raise PocketMatchError(f"{SCORE_FILE} ends early before pair {names[x]} {names[y]}")


_CRASH_ENV = "POLYPHAR_SIMILARITY_TEST_CRASH_ON"   # test hook: worker dies on this pocket ("*": all)


def _cabbage_worker(pocket: str, pm_bin_dir: str, record: str, failed: str) -> str | None:
    """cabbage-core on one pocket into the cache: `record` on success (returns None); on a
    cabbage-core failure writes `failed` and returns the reason. Both atomically."""
    if os.environ.get(_CRASH_ENV) in (Path(pocket).name, "*"):
        os._exit(1)
    return _task(lambda: _cabbage(Path(pocket), pm_bin_dir, Path(record), Path(failed)))


def _cabbage(pocket: Path, pm_bin_dir, record: Path, failed: Path) -> str | None:
    fd, tmp = tempfile.mkstemp(dir=record.parent, prefix=record.name + ".", suffix=".tmp")
    os.close(fd)
    try:
        try:
            cabbage_record(pocket, pm_bin_dir, tmp)
        except PocketMatchError as exc:
            reason = str(exc)
        else:
            os.chmod(tmp, 0o644)
            os.replace(tmp, record)
            return None
    finally:
        Path(tmp).unlink(missing_ok=True)
    _atomic_write_lines(failed, [reason + "\n"], encoding="utf-8")
    return reason


def _job_worker(paths: list[str], records: list[str], block_of: list[int], job: Job,
                pm_bin_dir: str, jobs_dir: str, fingerprint: str) -> int:
    return _task(lambda: _job(paths, records, block_of, job, pm_bin_dir, jobs_dir, fingerprint))


def _job(paths, records, block_of, job, pm_bin_dir, jobs_dir, fingerprint) -> int:
    t0 = time.monotonic()
    paths = [Path(p) for p in paths]
    names = [p.name for p in paths]
    jobs_dir = Path(jobs_dir)
    with tempfile.TemporaryDirectory(prefix=f"polyphar_pm_{job.tag}_") as tmp:
        score = run_pocketmatch(paths, pm_bin_dir, tmp, records=records)
        kept = 0

        def counted():
            nonlocal kept
            for line in _job_lines(score, names, block_of, job):
                kept += 1
                yield line
        txt = jobs_dir / f"{job.tag}.txt"
        _atomic_write_lines(txt, counted())
    marker = {"fingerprint": fingerprint, "pairs": kept, "bytes": txt.stat().st_size,
              "pockets": len(paths), "seconds": round(time.monotonic() - t0, 3)}
    _atomic_write_lines(jobs_dir / f"{job.tag}.done", [json.dumps(marker) + "\n"])
    return kept


def _reason(exc) -> str:
    return str(exc) if isinstance(exc, PocketMatchError) else f"{type(exc).__name__}: {exc}"


def _terminate(ex) -> None:
    """Cancel pending jobs and stop running ones now (an interrupted run leaves no
    workers and no PocketMatch children behind)."""
    procs = list((getattr(ex, "_processes", None) or {}).values())
    ex.shutdown(wait=False, cancel_futures=True)
    for p in procs:
        if p.is_alive():
            p.terminate()
    deadline = time.monotonic() + TERMINATE_GRACE_S
    for p in procs:
        p.join(timeout=max(0.0, deadline - time.monotonic()))
    for p in procs:   # a worker that did not exit in time is killed
        if p.is_alive():
            p.kill()
            p.join(timeout=10)


@dataclass
class PoolOutcome:
    failures: list[tuple[object, str]]   # (item, reason)
    isolated: list[object]               # items run alone after their pool broke twice


def _pool_pass(tasks, queue, nprocs: int, on_result, failures):
    """Run tasks from `queue` (indices into tasks) on one fresh spawn pool, keeping at most
    nprocs in flight, until the queue is empty or the pool breaks. Returns (indices that
    were in flight when the pool broke, number of tasks completed)."""
    import multiprocessing
    from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
    from concurrent.futures.process import BrokenProcessPool

    ctx = multiprocessing.get_context("spawn")   # workers start clean (no fork of a threaded parent)
    completed = 0
    with ProcessPoolExecutor(max_workers=max(1, min(nprocs, len(queue))), mp_context=ctx,
                             initializer=_worker_init) as ex:
        try:
            inflight: dict = {}
            while queue or inflight:
                while queue and len(inflight) < nprocs:
                    i = queue.popleft()
                    item, fn, args = tasks[i]
                    try:
                        inflight[ex.submit(fn, *args)] = i
                    except BrokenProcessPool:   # never ran: back to the queue
                        queue.appendleft(i)
                        return list(inflight.values()), completed
                finished, _ = wait(inflight, return_when=FIRST_COMPLETED)
                broke = []
                for fut in finished:
                    i = inflight.pop(fut)
                    try:
                        result = fut.result()
                    except BrokenProcessPool:
                        broke.append(i)
                        continue
                    except Exception as exc:
                        failures.append((tasks[i][0], _reason(exc)))
                    else:
                        if on_result is not None:
                            on_result(tasks[i][0], result)
                    completed += 1
                if broke:   # everything still in flight went down with it
                    return broke + list(inflight.values()), completed
            return [], completed
        except BaseException:
            _terminate(ex)
            raise


def _run_pool(tasks, nprocs: int, on_result=None) -> PoolOutcome:
    """Run tasks [(item, fn, args)] on spawn process pools; returns the failed items with
    their reasons. on_result(item, result) runs in this (parent) process. When a worker
    dies, the tasks in flight go back to a fresh full-size pool; a task in flight during
    a second break runs alone, so a crash is attributed to its own task. When pools keep
    breaking without finishing anything, or every solo run dies, the rest are failed
    without retrying. Any BaseException (Ctrl-C, SIGTERM) stops the workers and propagates."""
    from collections import deque

    failures: list = []
    queue = deque(range(len(tasks)))
    breaks: dict[int, int] = {}
    solo: list[int] = []
    empty = total = 0
    while queue:
        broken, completed = _pool_pass(tasks, queue, nprocs, on_result, failures)
        total += completed
        if not broken:
            continue
        empty = 0 if completed else empty + 1
        retry = []
        for i in broken:
            breaks[i] = breaks.get(i, 0) + 1
            (solo if breaks[i] >= 2 else retry).append(i)
        queue.extendleft(reversed(sorted(retry)))
        if total == 0 and empty >= SOLO_RETRY_PROBES:   # pools cannot run anything at all
            failures += [(tasks[i][0], REASON_DIED + "; not retried, every pool died")
                         for i in list(queue) + solo]
            return PoolOutcome(failures, [])
    isolated, tried, died = [], 0, 0
    for i in sorted(solo):
        if tried >= SOLO_RETRY_PROBES and died == tried:
            failures.append((tasks[i][0], REASON_DIED + "; not retried, every retry died"))
            continue
        tried += 1
        isolated.append(tasks[i][0])
        broken, _ = _pool_pass(tasks, deque([i]), 1, on_result, failures)
        if broken:
            failures.append((tasks[i][0], REASON_DIED))
            died += 1
    return PoolOutcome(failures, isolated)


def _bin_hashes(pm: Path) -> dict[str, str]:
    return {b: sha256_file(Path(pm) / b) for b in PM_BINARIES}


def _cache_paths(work: Path, path: Path, digest: str) -> tuple[Path, Path]:
    key = f"{digest}_{path.name}"   # the record embeds the name, so it is part of the key
    d = Path(work) / CABBAGE_DIR
    return d / f"{key}.cabbage", d / f"{key}.failed"


def _cache_valid(work: Path, cabbage_sha: str | None) -> bool:
    stamp = Path(work) / CABBAGE_DIR / CABBAGE_STAMP
    try:
        return cabbage_sha is not None and stamp.read_text().strip() == cabbage_sha
    except OSError:
        return False


def _cached(work: Path, paths, digests, cabbage_sha) -> tuple[dict[int, Path], dict[int, str]]:
    """(records, failures) already in the cabbage cache, by pocket index (read only)."""
    records, failed = {}, {}
    if not _cache_valid(work, cabbage_sha):
        return records, failed
    for i, (p, dg) in enumerate(zip(paths, digests)):
        rec, fl = _cache_paths(work, p, dg)
        if rec.is_file():
            records[i] = rec
        elif fl.is_file():
            failed[i] = fl.read_text(encoding="utf-8", errors="replace").strip()
    return records, failed


def _cabbage_phase(paths, digests, pm: Path, bins: dict, work: Path, nprocs: int):
    """Phase 1: cabbage records for every pocket, from the cache or computed in parallel.
    Returns (records by index, reasons by index for the pockets cabbage-core rejected or
    whose worker died)."""
    d = Path(work) / CABBAGE_DIR
    if not _cache_valid(work, bins[CABBAGE]):   # new binary: old records are stale
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)
        _atomic_write_lines(d / CABBAGE_STAMP, [bins[CABBAGE] + "\n"])
    records, failed = _cached(work, paths, digests, bins[CABBAGE])
    todo = [i for i in range(len(paths)) if i not in records and i not in failed]
    # Largest first, so a giant pocket (cabbage-core is superlinear: 42 s for a 265 kB
    # net2 pocket) does not start last; results do not depend on run order.
    todo.sort(key=lambda i: (-paths[i].stat().st_size, i))
    tasks = []
    for i in todo:
        rec, fl = _cache_paths(work, paths[i], digests[i])
        tasks.append((i, _cabbage_worker, (str(paths[i]), str(pm), str(rec), str(fl))))

    def collect(i, reason):
        if reason is None:
            records[i] = _cache_paths(work, paths[i], digests[i])[0]
        else:
            failed[i] = reason
    outcome = _run_pool(tasks, nprocs, collect)
    for i, reason in outcome.failures:
        failed[i] = reason
    return records, failed


def _run_jobs(plan: _Plan, records: list[Path], pm_bin_dir: Path, jobs_dir: Path,
              nprocs: int, progress) -> None:
    """Phase 2: the block-pair jobs not done yet; any failure is fatal (JobsFailed)."""
    todo = [j for j in plan.jobs if j.tag not in plan.done]
    if not todo:
        return
    block_of = [b for b, (s, e) in enumerate(plan.blocks) for _ in range(s, e)]
    tasks = []
    for job in todo:
        members = _members(job, plan.blocks)
        tasks.append((job.tag, _job_worker,
                      ([str(plan.paths[a]) for a in members], [str(records[a]) for a in members],
                       [block_of[a] for a in members], job, str(pm_bin_dir),
                       str(jobs_dir), plan.fingerprints[job.tag])))

    def advance(item, kept):
        if progress is not None:
            progress.advance(kept)
    outcome = _run_pool(tasks, nprocs, advance)
    if outcome.failures:
        raise JobsFailed(sorted(outcome.failures))


def _merge(plan: _Plan, jobs_dir: Path, out: Path) -> int:
    """Merge the job files into `out` in serial order; returns the line count."""
    index = {p.name: i for i, p in enumerate(plan.paths)}

    def keyed(job):
        with open(jobs_dir / f"{job.tag}.txt", encoding="latin-1", newline="") as fh:
            for line in fh:
                f = line.split(None, 2)
                yield (index[f[0].rstrip("_")], index[f[1].rstrip("_")]), line

    count = 0

    def lines():
        nonlocal count
        last = None
        for key, line in heapq.merge(*(keyed(j) for j in plan.jobs), key=lambda t: t[0]):
            if last is not None and key <= last:
                raise PocketMatchError(f"job outputs overlap or are out of order at {key}")
            last = key
            count += 1
            yield line
    _atomic_write_lines(out, lines())
    if count != plan.total:
        out.unlink(missing_ok=True)
        raise PocketMatchError(f"merged {count} pairs, expected {plan.total}")
    return count


@dataclass
class _Result:
    score: Path
    plan: _Plan
    skipped: int                 # pairs of jobs already done before this run
    dropped: dict[Path, str]     # pockets dropped in phase 1 -> reason


def _similarity(pockets, pm_bin_dir, work_dir, nprocs, progress=None) -> _Result:
    work, pm = Path(work_dir), Path(pm_bin_dir)
    paths = [Path(p) for p in pockets]
    bins = _bin_hashes(pm)
    digests = [sha256_file(p) for p in paths]
    records, failed = _cabbage_phase(paths, digests, pm, bins, work, nprocs)
    keep = [i for i in range(len(paths)) if i not in failed]
    plan = _plan([paths[i] for i in keep], nprocs, work, [digests[i] for i in keep], bins)
    jobs_dir = work / JOBS
    jobs_dir.mkdir(parents=True, exist_ok=True)
    current = {f"{j.tag}{s}" for j in plan.jobs for s in (".txt", ".done")}
    for f in jobs_dir.iterdir():   # job files of another plan (other nprocs or pockets)
        if f.suffix in (".txt", ".done") and f.name not in current:
            f.unlink(missing_ok=True)
    _run_jobs(plan, [records[i] for i in keep], pm, jobs_dir, nprocs, progress)
    out = work / SCORE_FILE
    if plan.jobs:
        _merge(plan, jobs_dir, out)
    else:
        _atomic_write_lines(out, [])
    return _Result(out, plan, plan.total - plan.remaining,
                   {paths[i]: failed[i] for i in sorted(failed)})


def run_pocketmatch_parallel(pockets, pm_bin_dir, work_dir, nprocs, progress=None) -> Path:
    """All-vs-all PocketMatch on `pockets` as block-pair jobs on `nprocs` processes;
    cabbage records come from the cache in work_dir/cabbage and finished jobs (matching
    .done markers in work_dir/jobs) are skipped. Pockets cabbage-core rejects are left
    out. Writes work_dir/PocketMatch_score.txt, byte-identical to a serial run, and
    returns it."""
    return _similarity(pockets, pm_bin_dir, work_dir, nprocs, progress).score


# --- scores --------------------------------------------------------------------------

def _gene_names(protein_info) -> dict[str, str]:
    import pandas as pd

    info = protein_info if isinstance(protein_info, pd.DataFrame) else pd.read_csv(
        protein_info, dtype=str, keep_default_na=False)
    out: dict[str, str] = {}
    for acc, gene in zip(info["UniProt ID"].astype(str), info["Gene Name"].astype(str)):
        if gene.strip() and gene.strip() != "-":
            out.setdefault(acc, gene)
    return out


def parse_scores(score_file, pockets, protein_info):
    """The Pocketsim table from a PocketMatch score file. `pockets` is the DrugScore table
    of the scored pockets (UniProt ID, Gene Name, Protein, Pocket, Pocket Filename);
    every name in the file must be one of its Pocket Filenames. Legacy filters: drop
    NULL____ scores, Pmin or Pmax == 1, and same-UniProt pairs. Rows are sorted by Pmax
    descending (stable); the index is the 0-based line number in the score file.

    Memory: the file is streamed into integer pocket codes and float scores (24 bytes a
    line); the per-pocket text columns are categoricals over the pocket table."""
    from array import array

    import numpy as np
    import pandas as pd

    names = [str(n) for n in pockets["Pocket Filename"]]
    code: dict[str, int] = {}
    for i, name in enumerate(names):
        if name in code:
            raise ValueError(f"pocket file name {name} is listed twice")
        code[name] = i
    a, b, pmin, pmax = array("i"), array("i"), array("d"), array("d")
    unknown: set[str] = set()
    nan = float("nan")
    with open(score_file, encoding="latin-1") as fh:
        for ln, line in enumerate(fh):
            f = line.split(None, 4)
            if len(f) < 4:
                raise ValueError(f"{score_file}:{ln + 1}: malformed line {line[:80]!r}")
            n1, n2 = f[0].rstrip("_"), f[1].rstrip("_")
            x, y = code.get(n1, -1), code.get(n2, -1)
            if x < 0 or y < 0:
                unknown.update(n for n, c in ((n1, x), (n2, y)) if c < 0)
            a.append(x)
            b.append(y)
            pmin.append(nan if f[2] == NULL else float(f[2]))
            pmax.append(nan if f[3] == NULL else float(f[3]))
    if unknown:
        unknown_s = sorted(unknown)
        raise ValueError(f"{len(unknown_s)} pocket name(s) in {score_file} are not in the "
                         f"pocket table: " + ", ".join(unknown_s[:20]))
    a, b = np.frombuffer(a, dtype=np.intc), np.frombuffer(b, dtype=np.intc)
    pmin, pmax = np.frombuffer(pmin, dtype=np.float64), np.frombuffer(pmax, dtype=np.float64)

    accs = np.asarray([str(x) for x in pockets["UniProt ID"]], dtype=object)
    acc_code, _ = pd.factorize(accs)
    keep = (~np.isnan(pmin) & ~np.isnan(pmax) & (pmin != 1) & (pmax != 1)
            & (acc_code[a] != acc_code[b]))
    idx = np.flatnonzero(keep)
    idx = idx[np.argsort(-pmax[idx], kind="stable")]
    ia, ib = a[idx], b[idx]
    del a, b, keep

    genes = _gene_names(protein_info)
    per_pocket = {
        "UniProt ID": list(accs),
        "Gene Name": [genes.get(acc) or str(g) for acc, g in zip(accs, pockets["Gene Name"])],
        "Protein": [str(x) for x in pockets["Protein"]],
        "Pocket": [str(x) for x in pockets["Pocket"]],
        "Pocket_Path": [f"{POCKETS_REL}/{acc}/{n}" for acc, n in zip(accs, names)],
    }
    cols = {}
    for suffix, sel in (("", ia), ("_2", ib)):
        for col, values in per_pocket.items():
            codes, cats = pd.factorize(np.asarray(values, dtype=object))
            cols[col + suffix] = pd.Categorical.from_codes(codes[sel], categories=cats)
    cols["Pmin"] = pmin[idx]
    cols["Pmax"] = pmax[idx]
    return pd.DataFrame(cols, index=pd.Index(idx))[POCKETSIM_COLUMNS]


# --- stage ---------------------------------------------------------------------------

def _inputs(layout):
    table, missing = pocket_table(layout.pockets)
    return table, _paths(layout.pockets, table), missing


def _guard(layout, paths, record: bool, fresh: bool = True) -> None:
    """The name guard; with record, its findings go to failed.tsv first (a fresh one, and
    any manifest removed, unless the caller already started them)."""
    rows, summary = _name_problems(paths)
    if not summary:
        return
    if record:
        d = layout.similarity
        if fresh:
            d.mkdir(parents=True, exist_ok=True)
            (d / "manifest.json").unlink(missing_ok=True)
            (d / "failed.tsv").write_text("item\tstage\treason\n")
        for item, reason in rows:
            append_failure(d, item, STAGE, reason)
    raise ValueError(summary)


def _pm_or_none(args) -> Path | None:
    try:
        return _pm_bin_dir(args)
    except FileNotFoundError:
        return None


def units(layout, args) -> WorkUnits:
    """Pairs still to compute: present pockets minus the ones the cabbage cache knows
    cabbage-core rejects, minus finished jobs. Fails on the name guard (in --dry-run too;
    otherwise the findings are written to failed.tsv). Without the PocketMatch binaries
    no job can be matched, so every pair counts."""
    _, paths, _ = _inputs(layout)
    _guard(layout, paths, record=not getattr(args, "dry_run", False))
    pm = _pm_or_none(args)
    bins = _bin_hashes(pm) if pm is not None else None
    digests = [sha256_file(p) for p in paths]
    _, failed = _cached(layout.similarity, paths, digests, bins and bins[CABBAGE])
    keep = [i for i in range(len(paths)) if i not in failed]
    plan = _plan([paths[i] for i in keep], args.nprocs, layout.similarity,
                 [digests[i] for i in keep], bins)
    return WorkUnits(plan.remaining, "pair")


def _versions(pm: Path) -> dict:
    return {"pocketmatch_commit": os.environ.get("PM_COMMIT") or "unknown",
            "pocketmatch_sha256": _bin_hashes(pm)}


def run(layout, args, progress) -> dict:
    d = layout.similarity
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    table, paths, missing = _inputs(layout)
    for item in missing:
        append_failure(d, item, STAGE, REASON_MISSING)
    _guard(layout, paths, record=True, fresh=False)   # before anything touches PocketMatch
    pm = _pm_bin_dir(args)
    (d / SCORE_FILE).unlink(missing_ok=True)
    (d / "Pocketsim.csv").unlink(missing_ok=True)

    try:
        res = _similarity(paths, pm, d, args.nprocs, progress)
    except JobsFailed as exc:
        for item, reason in exc.failures:
            append_failure(d, item, STAGE, reason)
        raise
    for path, reason in res.dropped.items():
        append_failure(d, _item(path), STAGE, reason)

    df = parse_scores(res.score, table, layout.protein_info)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".Pocketsim.csv.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            df.to_csv(fh)
        os.chmod(tmp, 0o644)
        os.replace(tmp, layout.pocketsim_csv)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    rows = len(df)
    del df
    shutil.rmtree(d / JOBS, ignore_errors=True)   # merged; ~100 bytes a pair

    counts = {"pockets": len(res.plan.paths), "pairs": res.plan.total,
              "jobs": len(res.plan.jobs), "skipped_existing": res.skipped, "rows": rows,
              "missing": len(missing), "cabbage_failed": len(res.dropped),
              "failed": len(missing) + len(res.dropped)}
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, STAGE, args={"nprocs": args.nprocs, "drugscore_cutoff": DRUGSCORE_CUTOFF},
                   inputs=[layout.drugscore_csv, layout.protein_info], versions=_versions(pm),
                   counts=counts, timing={})
    return counts
