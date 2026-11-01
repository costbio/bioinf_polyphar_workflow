"""Stage 05: pocket similarity with upstream PocketMatch v2.0, serial or block-pair parallel.

Inputs: 04_pockets/DrugScore.csv, rows with drugScore >= 0.5; each pocket file is
04_pockets/<UniProt ID>/<Pocket Filename>, found by its exact on-disk name. The pocket
list is ordered by (UniProt ID, Pocket Filename) so every run sees the same order.

PocketMatch (the upstream prebuilt 32-bit binaries in $PM_BIN_DIR) is run as the upstream
pipeline does: `Step0-cabbage-core ./<pocket>` (from the pocket's directory) for each
pocket into outfile.cabbage, then `Step0-END-FILE`, then `Step3-PM_serial outfile.cabbage`,
every exit code checked, in private temp dirs. cabbage-core's output depends only on the
pocket, so it runs once per pocket (in parallel) and each job concatenates the records.
PocketMatch records the file's base name, truncated to 32 characters, so names longer
than 32 or colliding fail before it runs. Its PocketMatch_pairs.txt (hundreds of bytes
per pair, not used) goes to /dev/null.

Parallelism: the ordered list is cut into k contiguous blocks B_0..B_{k-1}, k the smallest
with k(k-1)/2 >= nprocs. Job (i, j), i < j, runs all-vs-all on B_i ++ B_j and keeps the
cross pairs; the within-block pairs of B_i are kept by exactly one job (the ring
neighbour (i, i+1 mod k)). Each job writes 05_similarity/jobs/<i>-<j>.txt and then a
<i>-<j>.done marker holding a fingerprint of its pockets (names + file sha256); a rerun
skips jobs whose marker matches. The job files are merged in the order a serial run
emits pairs ((index(name1), index(name2)) in the ordered list), so the merged
PocketMatch_score.txt is byte-identical for every nprocs.

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
NULL = "NULL____"
POCKETSIM_COLUMNS = ["UniProt ID", "Gene Name", "Protein", "Pocket", "Pocket_Path",
                     "UniProt ID_2", "Gene Name_2", "Protein_2", "Pocket_2", "Pocket_Path_2",
                     "Pmin", "Pmax"]
REASON_DIED = "worker process died (BrokenProcessPool)"
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
    """DrugScore.csv rows with drugScore >= cutoff, one per pocket file, sorted by
    (UniProt ID, Pocket Filename). Raises FileNotFoundError, listing them, when files are
    missing (compared by exact name, so a case mismatch is caught on any filesystem)."""
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
    missing = []
    for acc, name in zip(df["UniProt ID"], df["Pocket Filename"]):
        if acc not in listing:
            d = pockets_dir / acc
            listing[acc] = set(os.listdir(d)) if d.is_dir() else set()
        if name not in listing[acc]:
            missing.append(f"{acc}/{name}")
    if missing:
        raise FileNotFoundError(f"{len(missing)} pocket file(s) listed in {path} are missing: "
                                + ", ".join(missing[:20]))
    return df


def _paths(pockets_dir, table) -> list[Path]:
    return [Path(pockets_dir) / a / f for a, f in zip(table["UniProt ID"], table["Pocket Filename"])]


def select_pockets(pockets_dir, cutoff: float = DRUGSCORE_CUTOFF) -> list[Path]:
    """Pocket files with drugScore >= cutoff, in the stage's deterministic order."""
    return _paths(pockets_dir, pocket_table(pockets_dir, cutoff))


def check_pocket_names(pockets) -> None:
    """PocketMatch keeps the first 32 characters of the base name and pads it with '_'.
    Raise ValueError, listing the names, when a name is longer than 32 characters,
    cannot be read back (whitespace, trailing '_'), or two names collide."""
    names = [Path(p).name for p in pockets]
    long = sorted({n for n in names if len(n) > PM_NAME_LIMIT})
    bad = sorted({n for n in names if n != n.strip() or any(c.isspace() for c in n)
                  or n.endswith("_")})
    seen: dict[str, list[str]] = {}
    for p in pockets:
        seen.setdefault(Path(p).name[:PM_NAME_LIMIT].rstrip("_"), []).append(str(p))
    collide = {k: v for k, v in seen.items() if len(v) > 1}
    problems = []
    if long:
        problems.append(f"{len(long)} pocket name(s) longer than {PM_NAME_LIMIT} characters "
                        f"(PocketMatch truncates them): " + ", ".join(long[:20]))
    if collide:
        problems.append(f"{len(collide)} pocket name(s) collide in PocketMatch: " + "; ".join(
            f"{k}: {', '.join(v)}" for k, v in list(collide.items())[:20]))
    if bad:
        problems.append("pocket name(s) with whitespace or a trailing '_': " + ", ".join(bad[:20]))
    if problems:
        raise ValueError("; ".join(problems))


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
    _check(proc, f"{CABBAGE} {p}")


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


def _plan(paths: list[Path], nprocs: int, work_dir: Path) -> _Plan:
    blocks = plan_blocks(len(paths), nprocs)
    jobs = plan_jobs(blocks) if len(paths) >= 2 else []
    digests = [sha256_file(p) for p in paths]
    jobs_dir = Path(work_dir) / JOBS
    fps, done = {}, {}
    for job in jobs:
        payload = {"blocks": [[[paths[a].name, digests[a]] for a in range(*blocks[b])]
                              for b in job.blocks],
                   "within": list(job.within)}
        fp = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        fps[job.tag] = fp
        marker = _read_done(jobs_dir / f"{job.tag}.done")
        if marker.get("fingerprint") == fp and (jobs_dir / f"{job.tag}.txt").is_file():
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


def _atomic_write_lines(path: Path, lines) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="latin-1", newline="") as fh:
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


def _cabbage_worker(pocket: str, pm_bin_dir: str, out: str) -> None:
    if os.environ.get(_CRASH_ENV) in (Path(pocket).name, "*"):
        os._exit(1)
    _task(lambda: cabbage_record(pocket, pm_bin_dir, out))


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
        _atomic_write_lines(jobs_dir / f"{job.tag}.txt", counted())
    marker = {"fingerprint": fingerprint, "pairs": kept, "pockets": len(paths),
              "seconds": round(time.monotonic() - t0, 3)}
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


def _run_pool(tasks, nprocs: int, on_result=None) -> list[tuple[str, str]]:
    """Run tasks [(item, fn, args)] on a spawn process pool; returns [(item, reason)] for
    the failed ones. on_result(item, result) runs in this (parent) process. Any
    BaseException (Ctrl-C, SIGTERM) stops the workers and propagates."""
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from concurrent.futures.process import BrokenProcessPool

    ctx = multiprocessing.get_context("spawn")   # workers start clean (no fork of a threaded parent)
    failures: list[tuple[str, str]] = []
    broken = []

    def done(task, result):
        if on_result is not None:
            on_result(task[0], result)

    if not tasks:
        return failures
    with ProcessPoolExecutor(max_workers=max(1, min(nprocs, len(tasks))), mp_context=ctx,
                             initializer=_worker_init) as ex:
        try:
            futures = {}
            for task in tasks:
                try:
                    futures[ex.submit(task[1], *task[2])] = task
                except BrokenProcessPool:   # a worker already died: retried below
                    broken.append(task)
            for fut in as_completed(futures):
                task = futures[fut]
                try:
                    result = fut.result()
                except BrokenProcessPool:
                    broken.append(task)
                    continue
                except Exception as exc:
                    failures.append((task[0], _reason(exc)))
                    continue
                done(task, result)
        except BaseException:
            _terminate(ex)
            raise
    # The pool died under the tasks still in flight: redo each alone in a fresh pool,
    # so a crash is attributed to its own task. When the first solo retries all die too,
    # workers cannot run at all: fail the rest without spawning a pool for each.
    tried = died = 0
    for task in sorted(broken, key=tasks.index):
        item, fn, args = task
        if tried >= SOLO_RETRY_PROBES and died == tried:
            failures.append((item, REASON_DIED + "; not retried, every retry died"))
            continue
        tried += 1
        with ProcessPoolExecutor(max_workers=1, mp_context=ctx, initializer=_worker_init) as ex:
            try:
                result = ex.submit(fn, *args).result()
            except BrokenProcessPool:
                failures.append((item, REASON_DIED))
                died += 1
                continue
            except Exception as exc:
                failures.append((item, _reason(exc)))
                continue
            except BaseException:
                _terminate(ex)
                raise
        done(task, result)
    return failures


def _item(path: Path) -> str:
    return f"{path.parent.name}/{path.name}"


def _run_jobs(plan: _Plan, pm_bin_dir: Path, jobs_dir: Path, nprocs: int, progress) -> None:
    """Phase 1: cabbage-core once per pocket the remaining jobs need, in parallel, into a
    private temp dir. Phase 2: the block-pair jobs, which concatenate those records."""
    todo = [j for j in plan.jobs if j.tag not in plan.done]
    if not todo:
        return
    block_of = [b for b, (s, e) in enumerate(plan.blocks) for _ in range(s, e)]
    # Largest first, so a giant pocket (cabbage-core is superlinear: 42 s for a 265 kB
    # net2 pocket) does not start last; result order does not depend on run order.
    size = {a: plan.paths[a].stat().st_size for job in todo for a in _members(job, plan.blocks)}
    needed = sorted(size, key=lambda a: (-size[a], a))
    with tempfile.TemporaryDirectory(prefix="polyphar_cabbage_") as tmp:
        rec = {a: str(Path(tmp) / f"{a}.cabbage") for a in needed}
        failures = _run_pool([(_item(plan.paths[a]), _cabbage_worker,
                               (str(plan.paths[a]), str(pm_bin_dir), rec[a])) for a in needed],
                             nprocs)
        if failures:
            raise JobsFailed(sorted(failures))

        tasks = []
        for job in todo:
            members = _members(job, plan.blocks)
            tasks.append((job.tag, _job_worker,
                          ([str(plan.paths[a]) for a in members], [rec[a] for a in members],
                           [block_of[a] for a in members], job, str(pm_bin_dir),
                           str(jobs_dir), plan.fingerprints[job.tag])))

        def advance(item, kept):
            if progress is not None:
                progress.advance(kept)
        failures = _run_pool(tasks, nprocs, advance)
    if failures:
        raise JobsFailed(sorted(failures))


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


def run_pocketmatch_parallel(pockets, pm_bin_dir, work_dir, nprocs, progress=None) -> Path:
    """All-vs-all PocketMatch on `pockets` as block-pair jobs on `nprocs` processes;
    finished jobs (matching .done markers in work_dir/jobs) are skipped. Writes
    work_dir/PocketMatch_score.txt, byte-identical to a serial run, and returns it."""
    work = Path(work_dir)
    jobs_dir = work / JOBS
    jobs_dir.mkdir(parents=True, exist_ok=True)
    plan = _plan([Path(p) for p in pockets], nprocs, work)
    current = {f"{j.tag}{s}" for j in plan.jobs for s in (".txt", ".done")}
    for f in jobs_dir.iterdir():   # job files of another plan (other nprocs or pockets)
        if f.suffix in (".txt", ".done") and f.name not in current:
            f.unlink(missing_ok=True)
    _run_jobs(plan, Path(pm_bin_dir), jobs_dir, nprocs, progress)
    out = work / SCORE_FILE
    if plan.jobs:
        _merge(plan, jobs_dir, out)
    else:
        _atomic_write_lines(out, [])
    return out


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
    descending (stable); the index is the 0-based line number in the score file."""
    import pandas as pd

    meta: dict[str, dict] = {}
    for row in pockets.to_dict("records"):
        name = str(row["Pocket Filename"])
        if name in meta:
            raise ValueError(f"pocket file name {name} is listed twice")
        meta[name] = row
    genes = _gene_names(protein_info)

    rows = []
    with open(score_file, encoding="latin-1") as fh:
        for i, line in enumerate(fh):
            f = line.split()
            if len(f) < 4:
                raise ValueError(f"{score_file}:{i + 1}: malformed line {line[:80]!r}")
            rows.append((i, f[0].rstrip("_"), f[1].rstrip("_"), f[2], f[3]))
    df = pd.DataFrame(rows, columns=["line", "n1", "n2", "Pmin", "Pmax"]).set_index("line")
    df.index.name = None
    unknown = sorted((set(df["n1"]) | set(df["n2"])) - meta.keys())
    if unknown:
        raise ValueError(f"{len(unknown)} pocket name(s) in {score_file} are not in the pocket "
                         f"table: " + ", ".join(unknown[:20]))

    df = df[(df["Pmin"] != NULL) & (df["Pmax"] != NULL)].copy()
    df["Pmin"] = pd.to_numeric(df["Pmin"])
    df["Pmax"] = pd.to_numeric(df["Pmax"])
    df = df[(df["Pmin"] != 1) & (df["Pmax"] != 1)]

    out = pd.DataFrame(index=df.index)
    for suffix, col in (("", "n1"), ("_2", "n2")):
        m = df[col].map(meta)
        acc = m.map(lambda r: str(r["UniProt ID"]))
        out["UniProt ID" + suffix] = acc
        out["Gene Name" + suffix] = [genes.get(a) or str(r["Gene Name"]) for a, r in zip(acc, m)]
        out["Protein" + suffix] = m.map(lambda r: str(r["Protein"]))
        out["Pocket" + suffix] = m.map(lambda r: str(r["Pocket"]))
        out["Pocket_Path" + suffix] = [f"{POCKETS_REL}/{a}/{n}" for a, n in zip(acc, df[col])]
    out["Pmin"] = df["Pmin"]
    out["Pmax"] = df["Pmax"]
    out = out[out["UniProt ID"] != out["UniProt ID_2"]]
    out = out.sort_values("Pmax", ascending=False, kind="mergesort")
    return out[POCKETSIM_COLUMNS]


# --- stage ---------------------------------------------------------------------------

def _inputs(layout):
    table = pocket_table(layout.pockets)
    return table, _paths(layout.pockets, table)


def units(layout, args) -> WorkUnits:
    _, paths = _inputs(layout)
    plan = _plan(paths, args.nprocs, layout.similarity)
    return WorkUnits(plan.remaining, "pair")


def _versions(pm: Path) -> dict:
    return {"pocketmatch_commit": os.environ.get("PM_COMMIT") or "unknown",
            "pocketmatch_sha256": {b: sha256_file(pm / b) for b in PM_BINARIES}}


def run(layout, args, progress) -> dict:
    d = layout.similarity
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    table, paths = _inputs(layout)
    check_pocket_names(paths)   # before anything touches PocketMatch
    pm = _pm_bin_dir(args)
    (d / SCORE_FILE).unlink(missing_ok=True)
    (d / "Pocketsim.csv").unlink(missing_ok=True)

    plan = _plan(paths, args.nprocs, d)
    skipped = plan.total - plan.remaining
    try:
        score = run_pocketmatch_parallel(paths, pm, d, args.nprocs, progress)
    except JobsFailed as exc:
        for item, reason in exc.failures:
            append_failure(d, item, STAGE, reason)
        raise

    df = parse_scores(score, table, layout.protein_info)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".Pocketsim.csv.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            df.to_csv(fh)
        os.chmod(tmp, 0o644)
        os.replace(tmp, layout.pocketsim_csv)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise

    counts = {"pockets": len(paths), "pairs": plan.total, "jobs": len(plan.jobs),
              "skipped_existing": skipped, "rows": len(df), "failed": 0}
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, STAGE, args={"nprocs": args.nprocs, "drugscore_cutoff": DRUGSCORE_CUTOFF},
                   inputs=[layout.drugscore_csv, layout.protein_info], versions=_versions(pm),
                   counts=counts, timing={})
    return counts
