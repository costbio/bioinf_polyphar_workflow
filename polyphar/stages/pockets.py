"""Stage 04: pocket detection with DoGSiteScorer 2 (default) or DoGSite3, run on the HOST
through a file queue.

Both tools are licensed for the host, so they never run in the container. This stage
(container side) writes one job file per model into 04_pockets/queue/ and waits while
bin/polyphar-dogsite-worker (host side) consumes them.

Tool (--pocket-tool, default $POLYPHAR_POCKET_TOOL or dogsite2): one tool per run folder.
It is recorded in 04_pockets/TOOL when a run starts, in every job file, .done and the
manifest; a run asking for the other tool than the run folder's existing pocket outputs
(TOOL, manifest, or queue records; records without a tool are DoGSiteScorer 2, the only
tool before) fails before anything is queued. Pocket selection (the `Selected` column):
  dogsite3   volume >= --min-volume (default 120 A^3, calibrated; DoGSite3 has no drugScore)
  dogsite2   drugScore >= --drugscore-cutoff (default 0.5; the 2023 rule)

Inputs (from 02_structures/Protein_info.csv, never globbed):
  AlphaFold == False  model = 03_pdbfixer/<UniProt>/<file name of File>
                      (a missing fixed file is a failure `no pdbfixer output`, not queued)
  AlphaFold == True   model = File itself (02_structures/<UniProt>/AF-<acc>-F1.pdb)
Outputs, `-o <out>` (the same base name for both tools):
  <out> = 04_pockets/<UniProt>/<pdb>_<chain>   (PDB chain; <UniProt>_<pdb>_<chain> when the
          same <pdb>_<chain> is a model of several accessions, a chimeric chain: PocketMatch
          keeps only base names, so they must stay unique)
          04_pockets/<UniProt>/AF-<acc>-F1     (AlphaFold model)
  dogsite3  -> <out>_desc.txt, <out>_P_<n>[_<m>]_res.pdb
  dogsite2  -> <out>_desc.txt, <out>_res_P_<n>[_<m>].pdb
  2023 wrote the dogsite2 base name as <UniProt>_<pdb>_<chain>.pdb, which with a current
  PDB entry's multi-character chain ID (e.g. P53350 9r1w chain AAA) exceeds PocketMatch's
  32-character limit; the prefix and the .pdb infix were dropped for that reason (decided
  2026-11-02), so dogsite2 pocket file names differ from 2023.
  04_pockets/DrugScore.csv   UniProt ID, Gene Name, Protein, Pocket, the tool's descriptor
                             columns, Pocket Filename (the file on disk), Tool, Selected;
                             minus pockets outside the target's SIFTS segment on a fusion
                             chain (filter_segments). The file name is kept from 2023.

Queue protocol (paths in job files are RELATIVE to the run folder, so the container's
/work and the host's path both resolve them); job id = <UniProt>_<pdb>_<chain> or
AF-<acc>-F1 for both tools:
  <id>.job                       {"model": ..., "out": ..., "size": <bytes>, "tool": ...},
                                 written by this stage through a unique temp name + os.replace
  <id>.running.<host>.<pid>      claimed by a worker (atomic os.rename)
  <id>.done                      {"exit": 0, "seconds": s, "dogsite_sha256": ..., "tool": ...,
                                  "tool_version": ...}
  <id>.failed                    {"exit": n, "stderr_tail": ...}
  WORKER.<host>.<pid>            worker heartbeat (touched every 30 s; a priority-1
                                 warning repeats hourly while none is fresh, and the stage
                                 fails after --worker-dead-after without any)
  CLOSED                         written here once every job is done or failed
A model is done when <out>_desc.txt exists and no job for it is queued or running.
Interrupting the wait is safe: the stage is resumable and queued jobs stay queued.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

from polyphar import modelling
from polyphar.common import WorkUnits, append_failure, write_manifest

STAGE = "pockets"
QUEUE = "queue"
POLL_S = 5.0
DEFAULT_WORKER_TIMEOUT = 600.0
DEFAULT_WORKER_DEAD_AFTER = 6 * 3600.0
NO_WORKER_REPEAT_S = 3600.0
REASON_NO_FIXED = "no pdbfixer output"
NO_WORKER = ("No DoGSiteScorer worker is consuming the queue — "
             "start bin/polyphar-dogsite-worker on the host")
TOOLS = ("dogsite3", "dogsite2")
DEFAULT_TOOL = "dogsite2"
LEGACY_TOOL = "dogsite2"        # queue records / manifests without a tool field
TOOL_ENV = "POLYPHAR_POCKET_TOOL"
TOOL_FILE = "TOOL"
DEFAULT_MIN_VOLUME = 120.0      # A^3, dogsite3 selection (scripts/calibrate_dogsite3_volume.py)
DEFAULT_DRUGSCORE_CUTOFF = 0.5  # dogsite2 selection (2023)
_DS3_POCKET = re.compile(r"^(.*)_(P_\d+(?:_\d+)*)_res\.pdb$")


def default_tool() -> str:
    tool = os.environ.get(TOOL_ENV) or DEFAULT_TOOL
    if tool not in TOOLS:
        raise SystemExit(f"polyphar: {TOOL_ENV}={tool!r} is not one of {', '.join(TOOLS)}")
    return tool


def _tool(args) -> str:
    return getattr(args, "pocket_tool", None) or DEFAULT_TOOL


def add_arguments(parser) -> None:
    parser.add_argument(
        "--pocket-tool", choices=TOOLS, default=default_tool(),
        help="pocket detection tool, one per run folder (default: $POLYPHAR_POCKET_TOOL or "
             f"{DEFAULT_TOOL}); the host worker must run the same tool")
    parser.add_argument(
        "--min-volume", type=float, default=DEFAULT_MIN_VOLUME, metavar="A3",
        help="dogsite3: select pockets (and subpockets) with volume >= this many A^3 "
             f"(default: {DEFAULT_MIN_VOLUME:g})")
    parser.add_argument(
        "--drugscore-cutoff", type=float, default=DEFAULT_DRUGSCORE_CUTOFF, metavar="S",
        help="dogsite2: select pockets with drugScore >= this "
             f"(default: {DEFAULT_DRUGSCORE_CUTOFF:g})")
    parser.add_argument(
        "--worker-timeout", type=float, default=DEFAULT_WORKER_TIMEOUT,
        help="seconds without a host worker heartbeat before a priority-1 warning, repeated "
             f"hourly (default: {DEFAULT_WORKER_TIMEOUT:g})")
    parser.add_argument(
        "--worker-dead-after", type=float, default=DEFAULT_WORKER_DEAD_AFTER,
        help="seconds without any host worker heartbeat before the stage fails (queued "
             f"jobs stay queued; default: {DEFAULT_WORKER_DEAD_AFTER:g})")


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


def _is_af(row) -> bool:
    return str(row["AlphaFold"]).strip().lower() == "true"


def out_names(info, tool: str) -> list[str]:
    """The `-o` base name of every Protein_info row, in order (see the module docstring).

    The same for both tools: `<pdb>_<chain>`, or `<UniProt>_<pdb>_<chain>` when that chain
    is a model of several accessions (a chimeric chain), because PocketMatch keeps only base
    names and they must stay unique. `tool` is kept in the signature because the pocket file
    names built from this base differ (see pocket_filename).
    """
    rows = info.to_dict("records")
    stems = [Path(str(r["File"])).stem for r in rows]
    owners: dict[str, set] = {}
    for r, stem in zip(rows, stems):
        if not _is_af(r):
            owners.setdefault(stem, set()).add(str(r["UniProt ID"]))
    shared = {stem for stem, accs in owners.items() if len(accs) > 1}
    return [f"{r['UniProt ID']}_{stem}" if not _is_af(r) and stem in shared else stem
            for r, stem in zip(rows, stems)]


def pocket_filename(tool: str, out_name: str, pocket: str) -> str:
    """The residue file of `pocket` (P_1, P_1_2, ...) for the `-o` base name `out_name`."""
    if tool == "dogsite2":
        return f"{out_name}_res_{pocket}.pdb"
    return f"{out_name}_{pocket}_res.pdb"


def _model_of(tool: str, filename: str) -> str:
    """The `-o` base name a pocket file belongs to (inverse of pocket_filename)."""
    if tool == "dogsite2":
        return filename.split("_res_", 1)[0]
    m = _DS3_POCKET.match(filename)
    return m.group(1) if m else filename


def _models_from(info, run_dir: Path | None, tool: str) -> tuple[list[Model], list[str]]:
    """(models, items without pdbfixer output). With run_dir None, existence and size
    of the model files are not checked (parsing only)."""
    models, missing = [], []
    for row, out_name in zip(info.to_dict("records"), out_names(info, tool)):
        acc = str(row["UniProt ID"])
        name = Path(str(row["File"])).name
        if _is_af(row):
            model = str(row["File"])
            job_id = Path(name).stem
            protein = Path(name).stem
        else:
            model = f"03_pdbfixer/{acc}/{name}"
            job_id = f"{acc}_{Path(name).stem}"
            protein = f"{row['PDB']}_{row['Chain ID']}".upper()
            if run_dir is not None and not (run_dir / model).is_file():
                missing.append(job_id)
                continue
        size = _size(run_dir / model) if run_dir is not None else 0
        models.append(Model(job_id=job_id, uniprot=acc,
                            gene=str(row["Gene Name"]), protein=protein, model=model,
                            out=f"04_pockets/{acc}/{out_name}", size=size))
    return models, missing


def select_models(layout, tool: str = DEFAULT_TOOL) -> tuple[list[Model], list[str]]:
    """Models listed in Protein_info.csv (in listed order), and the PDB rows whose
    pdbfixer output is missing."""
    return _models_from(_info(layout.protein_info), Path(layout.run_dir), tool)


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


def _pending(layout, tool: str) -> tuple[list[Model], list[Model], list[str]]:
    models, missing = select_models(layout, tool)
    state, _ = _queue_state(_queue_dir(layout))
    todo = [m for m in models if not _has_output(layout, m) or state.get(m.job_id) == "open"]
    return models, todo, missing


def _read_json(path: Path) -> dict:
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
        return obj if isinstance(obj, dict) else {}
    except (OSError, ValueError):
        return {}


def recorded_tool(layout) -> str | None:
    """The pocket tool of this run folder's existing pocket outputs: 04_pockets/TOOL, else
    the manifest and the queue records (no tool field = dogsite2), else None. Raises when
    the records disagree."""
    marker = layout.pockets / TOOL_FILE
    try:
        tool = marker.read_text(encoding="utf-8").strip()
        if tool:
            return tool
    except OSError:
        pass
    found: dict[str, str] = {}
    man = layout.pockets / "manifest.json"
    if man.is_file():
        tool = (_read_json(man).get("versions") or {}).get("pocket_tool") or LEGACY_TOOL
        found.setdefault(tool, str(man))
    qdir = _queue_dir(layout)
    try:
        names = sorted(os.listdir(qdir))
    except FileNotFoundError:
        names = []
    for name in names:
        if name.startswith(".") or not (name.endswith((".job", ".done", ".failed"))
                                        or ".running." in name):
            continue
        rec = _read_json(qdir / name)
        if rec or (qdir / name).is_file():
            found.setdefault(rec.get("tool") or LEGACY_TOOL, str(qdir / name))
    if len(found) > 1:
        detail = "; ".join(f"{t} ({p})" for t, p in found.items())
        raise RuntimeError(f"{layout.pockets} holds outputs of several pocket tools: {detail}; "
                           "one tool per run folder: delete 04_pockets/ and rerun")
    return next(iter(found), None)


def check_tool(layout, tool: str) -> None:
    """Raise when the run folder's pocket outputs were made with another tool."""
    rec = recorded_tool(layout)
    if rec is not None and rec != tool:
        raise RuntimeError(
            f"{layout.pockets} was made with --pocket-tool {rec}, this run asks for {tool}; "
            f"one pocket tool per run folder: rerun with --pocket-tool {rec}, or use a new "
            f"run folder (or delete 04_pockets/ and every later stage) for {tool}")


def calibration_tool(args) -> str:
    """The calibration variant of pockets (eta.Calibration `tool`)."""
    return _tool(args)


def selection_rule(args) -> dict:
    tool = _tool(args)
    if tool == "dogsite3":
        v = float(getattr(args, "min_volume", DEFAULT_MIN_VOLUME))
        return {"tool": tool, "column": "volume", "min": v, "text": f"volume >= {v:g} A^3"}
    c = float(getattr(args, "drugscore_cutoff", DEFAULT_DRUGSCORE_CUTOFF))
    return {"tool": tool, "column": "drugScore", "min": c, "text": f"drugScore >= {c:g}"}


def needs_rerun(layout, args) -> str | None:
    """Why a finished pockets stage must run again although no model is left (its
    selection rule changed), or None."""
    man = layout.pockets / "manifest.json"
    if not man.is_file():
        return None
    old = _read_json(man).get("selection_rule")
    new = selection_rule(args)
    if old != new:
        text = old.get("text") if isinstance(old, dict) else None
        return f"pocket selection changed ({text or 'none recorded'} -> {new['text']})"
    return None


def units(layout, args) -> WorkUnits:
    tool = _tool(args)
    check_tool(layout, tool)
    _, todo, _ = _pending(layout, tool)
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


def enqueue(layout, tool: str = DEFAULT_TOOL) -> list[str]:
    """Write a job file for every model that is not done yet; returns their job ids.
    Idempotent: queued or running jobs are left alone. A failed job, or a stale .done
    whose output is gone, is queued again. A CLOSED left by an earlier run is removed."""
    qdir = _queue_dir(layout)
    qdir.mkdir(parents=True, exist_ok=True)
    (qdir / "CLOSED").unlink(missing_ok=True)
    _, todo, _ = _pending(layout, tool)
    state, _ = _queue_state(qdir)
    for m in todo:
        s = state.get(m.job_id)
        if s == "open":
            continue
        for suffix in (".done", ".failed"):
            (qdir / f"{m.job_id}{suffix}").unlink(missing_ok=True)
        _write_json(qdir / f"{m.job_id}.job",
                    {"model": m.model, "out": m.out, "size": m.size, "tool": tool})
    return [m.job_id for m in todo]


def _failure_reason(info: dict) -> str:
    tail = " ".join(str(info.get("stderr_tail", "")).split())
    if info.get("reason"):   # e.g. the worker's `timeout after N s`
        return str(info["reason"]) + (f": {tail}" if tail else "")
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
                   jobs: list[str] | None = None,
                   worker_dead_after: float = DEFAULT_WORKER_DEAD_AFTER,
                   tool: str = DEFAULT_TOOL) -> dict:
    """Wait until every job in `jobs` (default: those queued or running now) is done or
    failed, advancing `progress` per job weighted by model size. Failed jobs go to
    failed.tsv. While no worker heartbeat is newer than `worker_timeout` seconds, a
    priority-1 notification is sent every NO_WORKER_REPEAT_S (hourly); after
    `worker_dead_after` seconds without any heartbeat the wait fails (RuntimeError; the
    jobs stay queued). Writes CLOSED at the end of a successful wait."""
    qdir = _queue_dir(layout)
    qdir.mkdir(parents=True, exist_ok=True)
    if jobs is None:
        state, _ = _queue_state(qdir)
        jobs = sorted(j for j, s in state.items() if s == "open")
    sizes = {m.job_id: m.size for m in select_models(layout, tool)[0]}
    remaining = list(dict.fromkeys(jobs))
    done = failed = 0
    started = time.time()
    warned_at = None
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
        now = time.time()   # the clock the worker stamps its heartbeat with (os.utime)
        silent = now - last_seen
        if silent > worker_dead_after:
            raise RuntimeError(
                f"no DoGSiteScorer worker heartbeat for {silent:.0f} s (--worker-dead-after "
                f"{worker_dead_after:g}); {len(remaining)} job(s) stay queued in "
                f"{qdir} for the next run")
        if silent > worker_timeout:
            if warned_at is None or now - warned_at >= NO_WORKER_REPEAT_S:
                warned_at = now
                print(f"{STAGE}: {NO_WORKER} (no heartbeat for {silent:.0f} s)", flush=True)
                _notify(progress, f"{STAGE}: no worker", NO_WORKER, priority=1)
        else:
            warned_at = None
        time.sleep(POLL_S)
    (qdir / "CLOSED").write_text("")
    return {"done": done, "failed": failed}


# --- parsing -------------------------------------------------------------------------

LEAD_COLUMNS = ["UniProt ID", "Gene Name", "Protein"]


def parse_pockets(pockets_dir, protein_info, tool: str = DEFAULT_TOOL):
    """The DrugScore.csv table for every model in Protein_info whose <out>_desc.txt
    exists: UniProt ID, Gene Name, Protein, Pocket, the tool's descriptor columns, and
    `Pocket Filename` = the pocket's residue file exactly as on disk (pocket_filename).
    Raises FileNotFoundError, listing them, when referenced files are missing (compared
    by exact name, so a case mismatch is caught on any filesystem)."""
    import pandas as pd

    pockets_dir = Path(pockets_dir)
    models, _ = _models_from(_info(protein_info), None, tool)
    frames, missing = [], []
    listing: dict[Path, set[str]] = {}
    for m in models:
        out_name = Path(m.out).name
        d = pockets_dir / m.uniprot
        desc = d / f"{out_name}_desc.txt"
        if not desc.is_file():
            continue
        df = pd.read_csv(desc, sep="\t", dtype={"name": str}).rename(columns={"name": "Pocket"})
        if d not in listing:
            listing[d] = set(os.listdir(d))
        files = [pocket_filename(tool, out_name, p) for p in df["Pocket"]]
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


# --- fusion / chimeric chains ------------------------------------------------------------
# A chain SIFTS maps to several accessions (a fusion construct, e.g. 2VE7) is kept whole for
# pdbfixer and DoGSiteScorer, but a pocket whose residues lie mostly (> 50 %) outside this
# accession's SIFTS segment(s) of that chain belongs to the partner protein: it is dropped
# (failed.tsv, counts["outside_segment"]) and never reaches DrugScore.csv. Pocket residues
# are numbered as in the pdbfixer output (renumbered 1..N), so each one is mapped back to
# its AUTHOR number through the coordinates of its heavy atoms in the extracted chain
# (02_structures/..., first model), which pdbfixer leaves unchanged. Residues no atom maps
# for count as not outside. The structures stage recovers a bound SIFTS leaves empty or
# None (an unobserved end) from the SEQRES range (SIFTS segments source = seqres); a pocket
# without any mapped residue, or a chain whose segment bounds still cannot be parsed, is
# kept and counted as segment_unverified. Single-accession chains and AlphaFold models are not checked.

SEGMENT_COLUMNS = ("SIFTS segments", "SIFTS other accessions")
MAX_OUTSIDE = 0.5
REASON_OUTSIDE = "pocket outside target segment"

# --- unmodelled pockets ------------------------------------------------------------------
# A pocket whose lining residues are traced main chain with no side chains (UNK, or a
# poly-ALA stub) is not a binding site: the chemistry that would line it was never
# modelled. Docking there scores a backbone groove, and PocketMatch SPINS FOREVER on one
# (a single UNK residue is enough). Dropped here, where DrugScore.csv is written, so it
# reaches neither `select`, nor `similarity`, nor `dock`.
#
# Measured over the 6,124 eligible pockets of AD Cluster1 and HCC net1/net2/net3
# (2026-11-05): the split is total. 6,118 have NO unmodelled residue at all, and 6 are
# 100 % unmodelled -- all of them 6BQ1 chain A, the entry that hung similarity. Nothing
# lies in between, so any threshold in (0, 1] selects the same six; a majority is used
# because a site whose lining is mostly unmodelled is already not a site.
MAX_UNMODELLED = 0.5
REASON_UNMODELLED = "pocket lining is unmodelled backbone (no side chains)"
_BOUND = re.compile(r"^(-?\d+)([A-Za-z]?)$")


def _bound(text: str):
    m = _BOUND.match(text.strip())
    return (int(m.group(1)), m.group(2).upper()) if m else None


def parse_segments(text: str):
    """[(begin, end)] from `PDB_BEG:PDB_END;...` (author numbers, optional insertion
    code), or None when any bound cannot be parsed (empty, None, ...) or there is none."""
    out = []
    for part in str(text).split(";"):
        beg, sep, end = part.partition(":")
        b, e = _bound(beg), _bound(end)
        if not sep or b is None or e is None:
            return None
        out.append((b, e))
    return out or None


def _le(a, b) -> bool:
    """a <= b for (number, insertion code); the same number with different insertion codes
    counts as <= both ways (their order is not reliable), so a bound is lenient."""
    return a[0] < b[0] if a[0] != b[0] else True


def in_segments(res, segments) -> bool:
    return any(_le(b, res) and _le(res, e) for b, e in segments)


def _coord_key(line: str):
    try:
        return tuple(int(round(float(line[i:i + 8]) * 1000)) for i in (30, 38, 46))
    except ValueError:
        return None


def _heavy(line: str) -> bool:
    elem = line[76:78].strip().upper()
    if elem:
        return elem not in ("H", "D")
    return not line[12:16].strip().upper().startswith(("H", "D"))


def _author_atoms(raw_pdb: Path) -> dict:
    """heavy-atom coordinate key -> (author residue number, insertion code), first model."""
    out = {}
    with open(raw_pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")) and _heavy(line):
                key = _coord_key(line)
                try:
                    res = (int(line[22:26]), line[26].strip().upper())
                except (ValueError, IndexError):
                    continue
                if key is not None:
                    out.setdefault(key, res)
    return out


_NEAR = [(dx, dy, dz) for dx in (0, -1, 1) for dy in (0, -1, 1) for dz in (0, -1, 1)]


def _lookup(atoms: dict, key):
    for dx, dy, dz in _NEAR:   # +-0.001 A: a coordinate written twice may round differently
        hit = atoms.get((key[0] + dx, key[1] + dy, key[2] + dz))
        if hit is not None:
            return hit
    return None


def outside_fraction(pocket_pdb: Path, atoms: dict, segments):
    """Fraction of the pocket's residues whose author number lies outside `segments`, or
    None when no residue of the pocket maps to the extracted chain."""
    votes: dict = {}
    with open(pocket_pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            rid = (line[21], line[22:27])
            v = votes.setdefault(rid, {})
            if not _heavy(line):
                continue
            key = _coord_key(line)
            hit = _lookup(atoms, key) if key is not None else None
            if hit is not None:
                v[hit] = v.get(hit, 0) + 1
    mapped = [max(v.items(), key=lambda kv: (kv[1], kv[0]))[0] for v in votes.values() if v]
    if not mapped:
        return None
    outside = sum(not in_segments(r, segments) for r in mapped)
    return outside / len(votes)


def filter_segments(layout, df, info, tool: str = DEFAULT_TOOL):
    """(kept DrugScore rows, [(item, reason)] dropped, n segment_unverified)."""
    missing = [c for c in SEGMENT_COLUMNS if c not in info.columns]
    if missing:
        raise ValueError(f"{layout.protein_info} has no {', '.join(missing)} column(s) (written "
                         f"by an older polyphar); rerun `polyphar structures` first")
    checked = {}
    for row, out_name in zip(info.to_dict("records"), out_names(info, tool)):
        others = str(row["SIFTS other accessions"]).strip()
        if _is_af(row) or not others:
            continue
        checked[(str(row["UniProt ID"]), out_name)] = (row, others)
    if not checked or df.empty:
        return df, [], 0
    keep, dropped, unverified = [], [], 0
    atoms_cache: dict = {}
    for acc, fname in zip(df["UniProt ID"], df["Pocket Filename"]):
        model = _model_of(tool, fname)
        hit = checked.get((acc, model))
        if hit is None:
            keep.append(True)
            continue
        row, others = hit
        segments = parse_segments(row["SIFTS segments"])
        frac = None
        if segments is not None:
            if model not in atoms_cache:
                raw = layout.run_dir / str(row["File"])
                atoms_cache[model] = _author_atoms(raw) if raw.is_file() else {}
            frac = outside_fraction(layout.pockets / acc / fname, atoms_cache[model], segments)
        if frac is None:
            unverified += 1
            keep.append(True)
        elif frac > MAX_OUTSIDE:
            keep.append(False)
            partners = ", ".join(others.split(";"))
            dropped.append((f"{acc}/{fname}", f"{REASON_OUTSIDE}: {frac * 100:.0f}% outside; "
                                              f"partner(s) {partners}"))
        else:
            keep.append(True)
    return df[keep].reset_index(drop=True), dropped, unverified


# --- run -----------------------------------------------------------------------------

def filter_unmodelled(layout, df):
    """(df without pockets whose lining is unmodelled, [(item, reason)]).

    Applies to every pocket, AlphaFold models included -- a predicted model has side
    chains everywhere, so it simply never matches."""
    if df.empty:
        return df, []
    keep, dropped = [], []
    for acc, fname in zip(df["UniProt ID"], df["Pocket Filename"]):
        path = layout.pockets / str(acc) / str(fname)
        try:
            frac, bad, total = modelling.unmodelled_fraction(path)
        except OSError:
            keep.append(True)      # unreadable here is the parse step's business
            continue
        if total and frac > MAX_UNMODELLED:
            keep.append(False)
            dropped.append((f"{acc}/{fname}",
                            f"{REASON_UNMODELLED}: {bad} of {total} residues"))
        else:
            keep.append(True)
    return df[keep].reset_index(drop=True), dropped


def _dogsite_versions(layout, models: list[Model], tool: str) -> tuple[str | None, str | None]:
    """(binary sha256, tool version) recorded by the worker for the current outputs.
    Raises when outputs come from different binaries or from another tool."""
    qdir = _queue_dir(layout)
    shas, versions = {}, set()
    for m in models:
        if _has_output(layout, m):
            done = _read_json(qdir / f"{m.job_id}.done")
            if done and (done.get("tool") or LEGACY_TOOL) != tool:
                raise RuntimeError(f"{m.job_id}: output made with {done.get('tool') or LEGACY_TOOL}"
                                   f", this run uses {tool}; one pocket tool per run folder")
            sha = done.get("dogsite_sha256")
            if sha:
                shas.setdefault(sha, []).append(m.job_id)
            if done.get("tool_version"):
                versions.add(str(done["tool_version"]))
    if len(shas) > 1:
        detail = "; ".join(f"{s[:12]}…: {', '.join(sorted(j)[:5])}" for s, j in shas.items())
        raise RuntimeError(f"pocket outputs come from different DoGSiteScorer binaries "
                           f"({detail}); rerun the stale ones with one binary")
    return next(iter(shas), None), (", ".join(sorted(versions)) or None)


def run(layout, args, progress) -> dict:
    tool = _tool(args)
    rule = selection_rule(args)
    d = layout.pockets
    check_tool(layout, tool)   # before anything in 04_pockets/ changes
    d.mkdir(parents=True, exist_ok=True)
    (d / TOOL_FILE).write_text(tool + "\n")
    (d / "manifest.json").unlink(missing_ok=True)   # outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    models, missing = select_models(layout, tool)
    for item in missing:
        append_failure(d, item, STAGE, REASON_NO_FIXED)

    jobs = enqueue(layout, tool)
    res = wait_for_queue(layout, progress,
                         getattr(args, "worker_timeout", DEFAULT_WORKER_TIMEOUT), jobs=jobs,
                         worker_dead_after=getattr(args, "worker_dead_after",
                                                   DEFAULT_WORKER_DEAD_AFTER), tool=tool)
    sha, version = _dogsite_versions(layout, models, tool)

    info = _info(layout.protein_info)
    df = parse_pockets(d, info, tool)
    df, outside, unverified = filter_segments(layout, df, info, tool)
    for item, reason in outside:
        append_failure(d, item, STAGE, reason)
    df, unmodelled = filter_unmodelled(layout, df)
    for item, reason in unmodelled:
        append_failure(d, item, STAGE, reason)
    df["Tool"] = tool
    if rule["column"] in df.columns:
        df["Selected"] = pd_numeric(df[rule["column"]]) >= rule["min"]
    elif df.empty:
        df["Selected"] = []
    else:
        raise ValueError(f"{tool} _desc.txt files have no {rule['column']} column; cannot "
                         f"apply {rule['text']}")
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".DrugScore.csv.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            df.to_csv(fh)
        os.replace(tmp, layout.drugscore_csv)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise

    _guard_names(d, df)

    counts = {"models": len(models) + len(missing), "queued": len(jobs),
              "done": res["done"], "skipped_existing": len(models) - len(jobs),
              "pockets": len(df), "selected": int(df["Selected"].sum()),
              "outside_segment": len(outside),
              "unmodelled_pockets": len(unmodelled),
              "segment_unverified": unverified, "failed": res["failed"] + len(missing)}
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, STAGE, args={"worker_timeout": getattr(args, "worker_timeout", None),
                                   "worker_dead_after": getattr(args, "worker_dead_after",
                                                                None),
                                   "pocket_tool": tool,
                                   "min_volume": getattr(args, "min_volume", None),
                                   "drugscore_cutoff": getattr(args, "drugscore_cutoff", None)},
                   inputs=[layout.protein_info],
                   versions={"pocket_tool": tool, "dogsite_sha256": sha,
                             "dogsite_version": version},
                   counts=counts, timing={}, extra={"selection_rule": rule})
    return counts


def _guard_names(d: Path, df) -> None:
    """similarity's PocketMatch name guard on the selected pockets (<= 32 characters, no
    collisions, no whitespace / trailing '_'), so a bad name fails here, before similarity:
    each finding goes to failed.tsv, then ValueError (no manifest; DrugScore.csv stays for
    inspection). Only a chimeric chain of a 10-character accession can realistically exceed
    32 characters with DoGSite3 names."""
    from polyphar.stages.similarity import _name_problems

    sel = df[df["Selected"].astype(bool)] if len(df) else df
    paths = [d / a / f for a, f in zip(sel.get("UniProt ID", []), sel.get("Pocket Filename", []))]
    rows, summary = _name_problems(paths)
    if summary:
        for item, reason in rows:
            append_failure(d, item, STAGE, reason)
        raise ValueError(f"selected pocket names PocketMatch cannot keep apart: {summary}")


def pd_numeric(series):
    import pandas as pd

    return pd.to_numeric(series, errors="coerce")
