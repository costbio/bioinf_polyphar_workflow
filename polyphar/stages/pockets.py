"""Stage 04: pocket detection with DoGSiteScorer, run on the HOST through a file queue.

DoGSiteScorer's license is bound to the host, so it never runs in the container.
This stage (container side) writes one job file per model into 04_pockets/queue/ and
waits while bin/polyphar-dogsite-worker (host side) consumes them.

Inputs (from 02_structures/Protein_info.csv, never globbed):
  AlphaFold == False  model = 03_pdbfixer/<UniProt>/<file name of File>
                      (a missing fixed file is a failure `no pdbfixer output`, not queued)
  AlphaFold == True   model = File itself (02_structures/<UniProt>/AF-<acc>-F1.pdb)
Outputs (DoGSiteScorer `-o <out>` writes <out>_desc.txt, <out>_res_<P>.pdb, ...):
  04_pockets/<UniProt>/<UniProt>_<pdb>_<chain>.pdb   for PDB chains
  04_pockets/<UniProt>/AF-<acc>-F1.pdb               for AlphaFold models
  04_pockets/DrugScore.csv                           the legacy columns (parse_pockets)

Queue protocol (paths in job files are RELATIVE to the run folder, so the container's
/work and the host's path both resolve them); job id = file name of <out> without .pdb:
  <id>.job                       {"model": ..., "out": ..., "size": <bytes>}, written by
                                 this stage through a unique temp name + os.replace
  <id>.running.<host>.<pid>      claimed by a worker (atomic os.rename)
  <id>.done                      {"exit": 0, "seconds": s, "dogsite_sha256": ...}
  <id>.failed                    {"exit": n, "stderr_tail": ...}
  WORKER.<host>.<pid>            worker heartbeat (touched every 30 s)
  CLOSED                         written here once every job is done or failed
A model is done when <out>_desc.txt exists and no job for it is queued or running.
Interrupting the wait is safe: the stage is resumable and queued jobs stay queued.
"""

from __future__ import annotations

import json
import os
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

from polyphar.common import WorkUnits, append_failure, write_manifest

STAGE = "pockets"
QUEUE = "queue"
POLL_S = 5.0
DEFAULT_WORKER_TIMEOUT = 600.0
REASON_NO_FIXED = "no pdbfixer output"
NO_WORKER = ("No DoGSiteScorer worker is consuming the queue — "
             "start bin/polyphar-dogsite-worker on the host")


def add_arguments(parser) -> None:
    parser.add_argument(
        "--worker-timeout", type=float, default=DEFAULT_WORKER_TIMEOUT,
        help="seconds without a host worker heartbeat before one warning notification "
             f"(default: {DEFAULT_WORKER_TIMEOUT:g})")


# --- inputs / outputs ----------------------------------------------------------------

@dataclass(frozen=True)
class Model:
    """One DoGSiteScorer work unit. `model` and `out` are relative to the run folder."""
    job_id: str
    uniprot: str
    gene: str
    protein: str        # DrugScore.csv `Protein`: 2OGQ_A, or AF-O14746-F1
    model: str
    out: str
    size: int


def _info(protein_info):
    import pandas as pd

    if isinstance(protein_info, pd.DataFrame):
        return protein_info
    path = Path(protein_info)
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; run `polyphar structures` first")
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _size(path: Path) -> int:
    try:
        return path.stat().st_size
    except OSError:
        return 0


def _models_from(info, run_dir: Path | None) -> tuple[list[Model], list[str]]:
    """(models, items without pdbfixer output). With run_dir None, existence and size
    of the model files are not checked (parsing only)."""
    models, missing = [], []
    for row in info.to_dict("records"):
        acc = str(row["UniProt ID"])
        name = Path(str(row["File"])).name
        if str(row["AlphaFold"]).strip().lower() == "true":
            model = str(row["File"])
            out_name = name
            protein = Path(name).stem
        else:
            model = f"03_pdbfixer/{acc}/{name}"
            out_name = f"{acc}_{name}"
            protein = f"{row['PDB']}_{row['Chain ID']}".upper()
            if run_dir is not None and not (run_dir / model).is_file():
                missing.append(f"{acc}_{Path(name).stem}")
                continue
        size = _size(run_dir / model) if run_dir is not None else 0
        models.append(Model(job_id=Path(out_name).stem, uniprot=acc,
                            gene=str(row["Gene Name"]), protein=protein, model=model,
                            out=f"04_pockets/{acc}/{out_name}", size=size))
    return models, missing


def select_models(layout) -> tuple[list[Model], list[str]]:
    """Models listed in Protein_info.csv (in listed order), and the PDB rows whose
    pdbfixer output is missing."""
    return _models_from(_info(layout.protein_info), Path(layout.run_dir))


def _queue_dir(layout) -> Path:
    return layout.pockets / QUEUE


def _queue_state(qdir: Path) -> tuple[dict, list[Path]]:
    """job id -> "done" | "failed" | "open", and the heartbeat files."""
    state, heartbeats = {}, []
    rank = {"open": 0, "failed": 1, "done": 2}
    try:
        names = os.listdir(qdir)
    except FileNotFoundError:
        return state, heartbeats
    for name in names:
        if name.startswith("WORKER."):
            heartbeats.append(qdir / name)
            continue
        if name.endswith(".job"):
            jid, s = name[:-4], "open"
        elif ".running." in name:
            jid, s = name.split(".running.", 1)[0], "open"
        elif name.endswith(".done"):
            jid, s = name[:-5], "done"
        elif name.endswith(".failed"):
            jid, s = name[:-7], "failed"
        else:
            continue
        if rank[s] >= rank[state.get(jid, "open")]:
            state[jid] = s
    return state, heartbeats


def _has_output(layout, m: Model) -> bool:
    return (layout.run_dir / f"{m.out}_desc.txt").is_file()


def _pending(layout) -> tuple[list[Model], list[Model], list[str]]:
    models, missing = select_models(layout)
    state, _ = _queue_state(_queue_dir(layout))
    todo = [m for m in models if not _has_output(layout, m) or state.get(m.job_id) == "open"]
    return models, todo, missing


def units(layout, args) -> WorkUnits:
    _, todo, _ = _pending(layout)
    return WorkUnits(len(todo), "bytes", [float(m.size) for m in todo])


# --- queue (container side) ----------------------------------------------------------

def _write_json(path: Path, obj) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(obj, fh)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def enqueue(layout) -> list[str]:
    """Write a job file for every model that is not done yet; returns their job ids.
    Idempotent: queued or running jobs are left alone. A failed job, or a stale .done
    whose output is gone, is queued again. A CLOSED left by an earlier run is removed."""
    qdir = _queue_dir(layout)
    qdir.mkdir(parents=True, exist_ok=True)
    (qdir / "CLOSED").unlink(missing_ok=True)
    _, todo, _ = _pending(layout)
    state, _ = _queue_state(qdir)
    for m in todo:
        s = state.get(m.job_id)
        if s == "open":
            continue
        for suffix in (".done", ".failed"):
            (qdir / f"{m.job_id}{suffix}").unlink(missing_ok=True)
        _write_json(qdir / f"{m.job_id}.job", {"model": m.model, "out": m.out, "size": m.size})
    return [m.job_id for m in todo]


def _read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _failure_reason(info: dict) -> str:
    tail = " ".join(str(info.get("stderr_tail", "")).split())
    return f"dogsite exit {info.get('exit', '?')}" + (f": {tail}" if tail else "")


def _notify(progress, title, message, priority=1) -> None:
    notifier = getattr(progress, "notifier", None)
    if notifier is None:
        return
    try:
        notifier(title, message, priority=priority)
    except Exception:
        pass


def wait_for_queue(layout, progress, worker_timeout: float = DEFAULT_WORKER_TIMEOUT, *,
                   jobs: list[str] | None = None) -> dict:
    """Wait until every job in `jobs` (default: those queued or running now) is done or
    failed, advancing `progress` per job weighted by model size. Failed jobs go to
    failed.tsv. If no worker heartbeat is newer than `worker_timeout` seconds, one
    priority-1 notification is sent and the wait goes on. Writes CLOSED at the end."""
    qdir = _queue_dir(layout)
    qdir.mkdir(parents=True, exist_ok=True)
    if jobs is None:
        state, _ = _queue_state(qdir)
        jobs = sorted(j for j, s in state.items() if s == "open")
    sizes = {m.job_id: m.size for m in select_models(layout)[0]}
    remaining = list(dict.fromkeys(jobs))
    done = failed = 0
    started = time.time()
    warned = False
    while True:
        state, heartbeats = _queue_state(qdir)
        still = []
        for jid in remaining:
            s = state.get(jid)
            if s == "done":
                done += 1
            elif s == "failed":
                failed += 1
                reason = _failure_reason(_read_json(qdir / f"{jid}.failed"))
                append_failure(layout.pockets, jid, STAGE, reason)
            else:
                still.append(jid)
                continue
            progress.advance(1, size=sizes.get(jid))
        remaining = still
        if not remaining:
            break
        last_seen = started
        for hb in heartbeats:
            try:
                last_seen = max(last_seen, hb.stat().st_mtime)
            except OSError:
                pass
        if not warned and time.time() - last_seen > worker_timeout:
            warned = True
            print(f"{STAGE}: {NO_WORKER}", flush=True)
            _notify(progress, f"{STAGE}: no worker", NO_WORKER, priority=1)
        time.sleep(POLL_S)
    (qdir / "CLOSED").write_text("")
    return {"done": done, "failed": failed}


# --- parsing -------------------------------------------------------------------------

LEAD_COLUMNS = ["UniProt ID", "Gene Name", "Protein"]


def parse_pockets(pockets_dir, protein_info):
    """The legacy DrugScore.csv table for every model in Protein_info whose
    <out>_desc.txt exists: UniProt ID, Gene Name, Protein, Pocket, the DoGSiteScorer
    descriptor columns, and `Pocket Filename` = <out name>_res_<Pocket>.pdb exactly as
    on disk. Raises FileNotFoundError, listing them, when referenced files are missing
    (compared by exact name, so a case mismatch is caught on any filesystem)."""
    import pandas as pd

    pockets_dir = Path(pockets_dir)
    models, _ = _models_from(_info(protein_info), None)
    frames, missing = [], []
    listing: dict[Path, set[str]] = {}
    for m in models:
        out_name = Path(m.out).name
        d = pockets_dir / m.uniprot
        desc = d / f"{out_name}_desc.txt"
        if not desc.is_file():
            continue
        df = pd.read_csv(desc, sep="\t").rename(columns={"name": "Pocket"})
        if d not in listing:
            listing[d] = set(os.listdir(d))
        files = [f"{out_name}_res_{p}.pdb" for p in df["Pocket"]]
        missing += [str(d / f) for f in files if f not in listing[d]]
        df.insert(0, "Protein", m.protein)
        df.insert(0, "Gene Name", m.gene)
        df.insert(0, "UniProt ID", m.uniprot)
        df["Pocket Filename"] = files
        frames.append(df)
    if missing:
        raise FileNotFoundError(f"{len(missing)} pocket file(s) listed in _desc.txt are "
                                f"missing: " + ", ".join(missing[:20]))
    if not frames:
        return pd.DataFrame(columns=LEAD_COLUMNS + ["Pocket", "Pocket Filename"])
    return pd.concat(frames, ignore_index=True)


# --- run -----------------------------------------------------------------------------

def _dogsite_sha(layout, models: list[Model]) -> str | None:
    """The DoGSiteScorer binary hash recorded by the worker for the current outputs.
    Raises when outputs come from different binaries."""
    qdir = _queue_dir(layout)
    shas = {}
    for m in models:
        if _has_output(layout, m):
            sha = _read_json(qdir / f"{m.job_id}.done").get("dogsite_sha256")
            if sha:
                shas.setdefault(sha, []).append(m.job_id)
    if len(shas) > 1:
        detail = "; ".join(f"{s[:12]}…: {', '.join(sorted(j)[:5])}" for s, j in shas.items())
        raise RuntimeError(f"pocket outputs come from different DoGSiteScorer binaries "
                           f"({detail}); rerun the stale ones with one binary")
    return next(iter(shas), None)


def run(layout, args, progress) -> dict:
    d = layout.pockets
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    models, missing = select_models(layout)
    for item in missing:
        append_failure(d, item, STAGE, REASON_NO_FIXED)

    jobs = enqueue(layout)
    res = wait_for_queue(layout, progress,
                         getattr(args, "worker_timeout", DEFAULT_WORKER_TIMEOUT), jobs=jobs)
    sha = _dogsite_sha(layout, models)

    df = parse_pockets(d, layout.protein_info)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".DrugScore.csv.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            df.to_csv(fh)
        os.replace(tmp, layout.drugscore_csv)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise

    counts = {"models": len(models) + len(missing), "queued": len(jobs),
              "done": res["done"], "skipped_existing": len(models) - len(jobs),
              "pockets": len(df), "failed": res["failed"] + len(missing)}
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, STAGE, args={"worker_timeout": getattr(args, "worker_timeout", None)},
                   inputs=[layout.protein_info], versions={"dogsite_sha256": sha},
                   counts=counts, timing={})
    return counts
