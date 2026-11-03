"""Stage 05: pocket similarity with upstream PocketMatch v2.0, serial or block-pair parallel.

Inputs: 04_pockets/DrugScore.csv, rows with Selected == True (the pockets stage's rule for its
tool: dogsite3 volume >= --min-volume, dogsite2 drugScore >= 0.5), or, for a legacy table
without that column, drugScore >= 0.5 (so the 2023 net2 regression still applies); each
pocket file is 04_pockets/<UniProt ID>/<Pocket Filename>, found by its exact on-disk name.
The pocket list is ordered by (UniProt ID, Pocket Filename) so every run sees the same order. A
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
tail (and cached as <key>.failed, so a resume does not rerun it) and dropped. A kill from
outside (SIGKILL/OOM, SIGTERM, SIGINT, SIGHUP, SIGQUIT, or exit 128+those) or a failure
while the worker is being stopped is recorded for this run only and retried next run.

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
Pocket_Path relative to the run folder), failed.tsv, manifest.json (written last). The
manifest records a stage fingerprint (surviving pockets and their sha256, PocketMatch
binaries, cutoff, DrugScore.csv, Protein_info.csv; not nprocs) and the output sizes; a
rerun on the same inputs reports 0 units and changes nothing.
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

from polyphar import pool
from polyphar.common import PM_NAME_LIMIT, WorkUnits, append_failure, sha256_file, write_manifest
from polyphar.pool import REASON_DIED, TaskError
from polyphar.pool import run_pool as _run_pool
from polyphar.pool import task as _task

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
REASON_MISSING = "pocket file missing"


def add_arguments(parser) -> None:
    parser.add_argument("--pm-bin-dir", default=None,
                        help="directory holding the PocketMatch binaries (default: $PM_BIN_DIR)")


class PocketMatchError(TaskError):
    """A PocketMatch step failed or wrote unexpected output (a plain message, so it
    always pickles across the process boundary; `returncode` does not cross it)."""

    def __init__(self, message, returncode=None):
        super().__init__(message)
        self.returncode = returncode


# Signals that come from outside (OOM killer, docker stop, Ctrl-C), not from the tool:
# such a failure says nothing about the pocket and is never cached. A tool's own crash
# (SIGABRT / SIGSEGV, e.g. cabbage-core on a one-atom pocket) is deterministic and is.
EXTERNAL_SIGNALS = frozenset({signal.SIGHUP, signal.SIGINT, signal.SIGQUIT,
                              signal.SIGKILL, signal.SIGTERM})


def _external_kill(returncode) -> bool:
    if returncode is None:
        return False
    return (-returncode if returncode < 0 else returncode - 128) in EXTERNAL_SIGNALS


class JobsFailed(RuntimeError):
    def __init__(self, failures: list[tuple[str, str]]):
        self.failures = failures
        first = "; ".join(f"{item}: {reason}" for item, reason in failures[:3])
        super().__init__(f"{len(failures)} PocketMatch job(s) failed: {first}")


# --- inputs --------------------------------------------------------------------------

def _truthy(series):
    """A `Selected` column as read back from CSV (bool, or the strings True/False)."""
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))


def _drugscore_table(pockets_dir):
    import pandas as pd

    path = Path(pockets_dir) / "DrugScore.csv"
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; run `polyphar pockets` first")
    text_cols = ["UniProt ID", "Gene Name", "Protein", "Pocket", "Pocket Filename"]
    df = pd.read_csv(path, index_col=0, dtype={c: str for c in text_cols})
    df[text_cols] = df[text_cols].fillna("")
    return df


def _selection(df, cutoff: float) -> str:
    return "Selected" if "Selected" in df.columns else f"drugScore >= {cutoff:g}"


def selection_used(pockets_dir, cutoff: float = DRUGSCORE_CUTOFF) -> str:
    """How pockets are chosen from this DrugScore.csv: "Selected" (the column the pockets
    stage writes by its tool's rule) or, for a legacy table, "drugScore >= <cutoff>"."""
    return _selection(_drugscore_table(pockets_dir), cutoff)


def pocket_table(pockets_dir, cutoff: float = DRUGSCORE_CUTOFF):
    """(table, missing): the selected DrugScore.csv rows whose pocket file exists, one per
    file, sorted by (UniProt ID, Pocket Filename); and the listed files that do not exist,
    as `<UniProt ID>/<Pocket Filename>` (compared by exact name, so a case mismatch counts
    as missing on any filesystem). Selected = the `Selected` column when present (written
    by the pockets stage per tool: dogsite3 volume, dogsite2 drugScore), else the legacy
    drugScore >= cutoff (DrugScore.csv of 2023 and of earlier polyphar versions)."""
    import pandas as pd

    pockets_dir = Path(pockets_dir)
    df = _drugscore_table(pockets_dir)
    if "Selected" in df.columns:
        df = df[_truthy(df["Selected"])]
    else:
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
    """Existing selected pocket files (pocket_table), in the stage's deterministic order."""
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
    rc = proc.returncode
    if rc == 0:
        return
    tail = (proc.stderr or b"")[-500:].decode(errors="replace").strip()
    if _external_kill(rc):
        sig = signal.Signals(-rc if rc < 0 else rc - 128)
        status = f"killed by signal {sig.value} ({sig.name}, exit {rc}); not cached, retried next run"
    elif rc < 0:
        status = f"crashed with signal {-rc} ({signal.Signals(-rc).name})"
    else:
        status = f"failed (exit {rc})"
    raise PocketMatchError(f"{what} {status}" + (f": {tail}" if tail else ""), rc)


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
            # only the tool's own verdict on this pocket is cached: not a kill from
            # outside, and nothing while this worker is being stopped
            cache = not (pool.stopping() or _external_kill(exc.returncode))
        else:
            os.chmod(tmp, 0o644)
            os.replace(tmp, record)
            return None
    finally:
        Path(tmp).unlink(missing_ok=True)
    if cache:
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
    digests: list[str]           # sha256 of plan.paths


def _similarity(pockets, pm_bin_dir, work_dir, nprocs, progress=None, digests=None) -> _Result:
    work, pm = Path(work_dir), Path(pm_bin_dir)
    paths = [Path(p) for p in pockets]
    bins = _bin_hashes(pm)
    digests = digests if digests is not None else [sha256_file(p) for p in paths]
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
                   {paths[i]: failed[i] for i in sorted(failed)}, [digests[i] for i in keep])


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
    every name in the file must be one of its Pocket Filenames. Filters: drop NULL____
    scores, same-UniProt pairs, and pairs on the same chain stem (the Protein column) --
    a chimeric chain SIFTS maps to two accessions is extracted into both accession
    folders, so the same physical atoms appear twice under different accessions and would
    otherwise self-compare at Pmin == Pmax == 1. Rows are sorted by Pmax descending
    (stable); the index is the 0-based line number in the score file.

    The 2023 script (394a337) also dropped every pair with Pmin == 1 or Pmax == 1,
    alongside a working `UniProt ID != UniProt ID_2` filter (added in 90f56f2); only the
    very first version, c2132f0, lacked that filter. So the == 1 drop was a deliberate
    filter of its own, not a self-comparison workaround. It is gone because, measured over
    every 2023 PocketMatch_score.txt on the drive, it discarded real cross-protein
    containment hits above the selection threshold -- e.g. CDC25C `P30307_3bzi_a` P_4
    against PLK1 `P53350_5j19_a` P_5, 210 of 276 elements matched, Pmin 1.0 but
    Pmax 0.76, present in net1 and net2.

    Caveat for anyone reading the output: because PocketMatch normalises Pmin by the
    smaller pocket and Pmax by the larger, Pmin >= Pmax always, so Pmax == 1 implies
    Pmin == 1 and the old mask was effectively just `Pmin != 1`. No cross-protein
    Pmax == 1 pair exists in any 2023 subtype, and most of the ~49.6k restored rows are
    degenerate containment (a 3-element pocket inside a 682,696-element AlphaFold
    whole-surface blob at Pmax 4e-06). Pmin == 1 on its own is therefore a near-meaningless
    similarity signal; judge pairs on Pmax.

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
    # A chimeric chain is extracted once per accession, so two accessions sharing a chain
    # stem (<pdb>_<chain>, or AF-<acc>-F1, which cannot collide) are the same atoms.
    stems = np.asarray([str(x) for x in pockets["Protein"]], dtype=object)
    stem_code, _ = pd.factorize(stems)
    keep = (~np.isnan(pmin) & ~np.isnan(pmax)
            & (acc_code[a] != acc_code[b]) & (stem_code[a] != stem_code[b]))
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

def eligible_csv(layout):
    return layout.targets / "eligible_pockets.csv"


def _inputs(layout):
    """The pockets to compare: the panel's eligible set, written by `select`.

    This stage used to run over every selected pocket, because the old hub-and-neighbours
    selection was computed from the resulting network. Selection is now panel filters
    alone, and similarity is only used post hoc (are a drug's hit targets more mutually
    similar than a random same-size panel set?), so comparing anything beyond the
    eligible pockets is wasted work: on net1, ~179k pairs instead of 67.2M.
    """
    import csv as _csv

    # pocket_table first, so a run folder missing BOTH inputs is told to run the
    # earlier stage rather than sent to `select` for want of a file `select` cannot
    # produce either.
    table, missing = pocket_table(layout.pockets)
    path = eligible_csv(layout)
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found: run `select` first (similarity now runs on the panel's "
            f"eligible pockets, after selection, not before it)")
    with open(path, newline="", encoding="utf-8") as fh:
        # keyed on the pair: a filename-only key can admit another accession's pocket
        keep = {(r["UniProt_ID"], r["Pocket Filename"]) for r in _csv.DictReader(fh)}
    mask = [(a, f) in keep for a, f in zip(table["UniProt ID"], table["Pocket Filename"])]
    table = table[mask].reset_index(drop=True)
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


OUTPUTS = (SCORE_FILE, "Pocketsim.csv")


def selection_record(layout) -> dict:
    """The pocket selection this stage used, for its manifest and fingerprint: a legacy
    DrugScore.csv -> {"drugscore_cutoff": 0.5}; one with `Selected` -> {"selection":
    "Selected", "selection_rule": the pockets manifest's selection_rule (None when
    unreadable)}."""
    if selection_used(layout.pockets) != "Selected":
        return {"drugscore_cutoff": DRUGSCORE_CUTOFF}
    rule = None
    try:
        rule = json.loads((layout.pockets / "manifest.json").read_text()).get("selection_rule")
    except (OSError, ValueError, AttributeError):
        pass
    return {"selection": "Selected", "selection_rule": rule}


def _fingerprint_payload(layout, survivors, digests, bins) -> dict:
    """What the finished outputs depend on, independent of nprocs: the ordered surviving
    pockets (item, sha256), the PocketMatch binaries, the pocket selection (selection_record;
    for a legacy table exactly the payload of earlier versions), DrugScore.csv and
    Protein_info.csv."""
    return {"pockets": [[_item(p), dg] for p, dg in zip(survivors, digests)],
            "pocketmatch_sha256": bins, **selection_record(layout),
            "drugscore_sha256": sha256_file(layout.drugscore_csv),
            "protein_info_sha256": sha256_file(layout.protein_info)}


def _stage_fingerprint(layout, survivors, digests, bins) -> str:
    payload = _fingerprint_payload(layout, survivors, digests, bins)
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _finished(layout, paths, digests, pm) -> dict | None:
    """The recorded counts when the stage already finished on exactly these inputs:
    manifest and both outputs present (sizes as written) and the stage fingerprint equal
    to the current one. Every pocket must have a cached cabbage verdict (one dropped for
    a worker death or an outside kill has none, so the stage is redone and it is retried)."""
    d = layout.similarity
    if pm is None:
        return None
    try:
        manifest = json.loads((d / "manifest.json").read_text())
        sizes = {name: (d / name).stat().st_size for name in OUTPUTS}
    except (OSError, ValueError):
        return None
    if manifest.get("outputs") != sizes or "fingerprint" not in manifest:
        return None
    bins = _bin_hashes(pm)
    records, failed = _cached(d, paths, digests, bins[CABBAGE])
    if len(records) + len(failed) != len(paths):
        return None
    keep = [i for i in range(len(paths)) if i in records]
    fp = _stage_fingerprint(layout, [paths[i] for i in keep], [digests[i] for i in keep], bins)
    if manifest["fingerprint"] != fp:
        return None
    return dict(manifest["counts"])


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
    if _finished(layout, paths, digests, pm) is not None:
        return WorkUnits(0, "pair")   # finished on these inputs: run() changes nothing
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
    table, paths, missing = _inputs(layout)
    digests = None
    if not _name_problems(paths)[1]:
        digests = [sha256_file(p) for p in paths]
        done = _finished(layout, paths, digests, _pm_or_none(args))
        if done is not None:   # finished on exactly these inputs: touch nothing
            return done
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    for item in missing:
        append_failure(d, item, STAGE, REASON_MISSING)
    _guard(layout, paths, record=True, fresh=False)   # before anything touches PocketMatch
    pm = _pm_bin_dir(args)
    (d / SCORE_FILE).unlink(missing_ok=True)
    (d / "Pocketsim.csv").unlink(missing_ok=True)

    try:
        res = _similarity(paths, pm, d, args.nprocs, progress, digests)
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
    versions = _versions(pm)
    fingerprint = _stage_fingerprint(layout, res.plan.paths, res.digests,
                                     versions["pocketmatch_sha256"])
    write_manifest(d, STAGE, args={"nprocs": args.nprocs, **selection_record(layout)},
                   inputs=[layout.drugscore_csv, layout.protein_info], versions=versions,
                   counts=counts, timing={},
                   extra={"fingerprint": fingerprint,
                          "outputs": {name: (d / name).stat().st_size for name in OUTPUTS}})
    return counts
