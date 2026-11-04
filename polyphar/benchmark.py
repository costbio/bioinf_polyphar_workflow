"""Per-machine calibration (`polyphar benchmark`) and whole-run ETA table (`polyphar estimate`).

benchmark times one small fixture workload per stage on ONE process and stores
seconds-per-unit in <POLYPHAR_HOME>/calibration.json, in the same unit basis the stages'
units() use (fix: input bytes; similarity: pairs; dock: ligands; uniprot / structures:
requests). `pockets` is never run here: DoGSite3 / DoGSiteScorer is host-only, so
bin/polyphar-dogsite-worker --benchmark calibrates it.

estimate prints the stages after the last completed one: the next stage exactly (from its
units()), later ones projected from PRIORS, with observed values from the run folder
replacing the priors as soon as the stage that produces them has finished.

Heavy imports (pandas, openmm, pdbfixer, openbabel) stay inside functions: cli.build_parser
imports every stage module at parse time.
"""

from __future__ import annotations

import math
import os
import shutil
import tempfile
import time
from datetime import datetime
from types import SimpleNamespace
from pathlib import Path

from polyphar.common import WorkUnits

BENCH_DIR = Path(__file__).resolve().parent / "data" / "bench"

# Bench chains for `fix` (raw extracted PDB chains of clearly different sizes: 59 kB,
# 116 kB, 280 kB). The worker's BENCH_CHAINS is the same list.
FIX_CHAINS = ("6l49_A.pdb", "2mjw_B.pdb", "5nmm_A.pdb")
DOCK_RECEPTOR = "2mjw_A_fixed.pdb"
DOCK_POCKET = "2mjw_A_pocket.pdb"
DOCK_LIGANDS = "edrug3d_bench3.sdf"
DOCK_SEED = 42
UNIPROT_REQUESTS = 5
UNIPROT_BATCH = 100   # = uniprot.BATCH (tested): one bench request is one stage request
# (accession, PDB, author chain): five live mmCIF downloads + chain extraction; each pair
# was verified against the entry's own UniProt annotation.
STRUCTURE_JOBS = (("P49450", "6l49", "A"), ("P53350", "5nmm", "A"), ("P53350", "4lkm", "A"),
                  ("P53350", "1q4k", "B"), ("P43220", "3iol", "A"))

# Dock library-mix correction. The bench ligands are typical small molecules (24 heavy atoms, the median
# of the library), but the real library costs more per ligand: a 100-ligand
# sample docked on this workstation (e-Drug3D_2083, exhaustiveness 12) had
# mean 71.5 s and median 45 s per ligand. The measured cost is scaled by mean/median; the
# in-run 5 % revised ETA corrects the rest.
DOCK_MIX_FACTOR = 71.5 / 45.0
DOCK_MIX_SOURCE = ("100-ligand sample of e-Drug3D_2083 on the dev workstation, 2026-11-01: "
                   "mean 71.5 s / median 45 s per ligand")

# --- projection priors -------------------------------------------------------------------
# Source: the 2023 HCC net2 run in $ADEP_DATA/disease_subtypes/Liver_cancer/net2 and
# Gene_Disease_Association/target_selection_results/liver_cancer_net2_final_targets.csv,
# computed once on 2026-11-01. Counts are on the legacy extraction (first chain per PDB
# entry, no SIFTS), so they are only priors for the v2 pipeline.
#   156 accessions; 1145 extracted PDB chains (mean 336,561 bytes), 156 AlphaFold models
#   (mean 333,980 bytes), 1145 pdbfixer outputs (mean 333,660 bytes); DrugScore.csv has 9,479
#   pockets of which 1,191 have drugScore >= 0.5, from only 495 models (2023 got pockets for
#   68 of 156 accessions), so kept pockets per submitted model = 1191 / (1145 + 156);
#   8 rows in final_targets.csv.
PRIORS = {
    "accessions": 156,
    "pdb_chains_per_accession": 1145 / 156,
    "af_models_per_accession": 156 / 156,
    "pdb_chain_bytes": 336561,
    "af_model_bytes": 333980,
    "fixed_chain_bytes": 333660,
    "kept_pockets_per_model": 1191 / (1145 + 156),
    "targets": 316,
    "ligands": 2083,         # e-Drug3D_2083 source records (what `ligands` processes)
    "ligands_kept": 1960,    # what the dockability filter keeps (what `dock` docks)
}
PRIORS_SOURCE = "2023 HCC net2 run"
# kept (selected) pockets per submitted model, by pocket tool. dogsite2: the PRIORS value
# (drugScore >= 0.5, 2023 net2). dogsite3: volume >= 300 A^3 on DoGSite3 1.3.0 output for 50
# random 2023 net2 models (40 pdbfixer chains + 10 AlphaFold models), 2026-11-02: 11 of 523
# pockets in 50 models.
KEPT_POCKETS_PRIOR = {"dogsite2": PRIORS["kept_pockets_per_model"], "dogsite3": 11 / 50}

BENCH_STAGES = ("uniprot", "structures", "fix", "pockets", "similarity", "dock")
# share of selected pockets that survive the panel filters, for the pre-run estimate
# only (measured on net1: 1,367 of 11,599). Pairs scale as n^2, so a stale value
# makes the similarity ETA wrong quadratically. The real count comes from `select`.
ELIGIBLE_FRACTION = 0.118
HOST_CALIBRATED = ("host-calibrated: DoGSite3 / DoGSiteScorer runs on the host only; calibrate "
                   "with bin/polyphar-dogsite-worker --benchmark --tool <tool> "
                   "(bin/polyphar-docker does it)")


class Skip(Exception):
    """A stage that cannot be benchmarked here; the message is the reason."""


# --- benchmark workloads (each returns (seconds_per_unit, extra calibration fields)) -------

def _bench_fix(tmp: Path):
    from polyphar.stages import fix

    try:
        import openmm.app  # noqa: F401  (imports stay out of the timed region)
        import pdbfixer  # noqa: F401
    except ImportError as exc:
        raise Skip(f"pdbfixer/OpenMM not importable: {exc}") from None
    fix.reference_platform()
    src_dir, out_dir = tmp / "structures", tmp / "pdbfixer"
    total_s = total_b = 0.0
    for name in FIX_CHAINS:
        src = src_dir / "X" / name
        src.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(BENCH_DIR / name, src)
        t0 = time.monotonic()
        fix.fix_one(src, src_dir, out_dir)
        total_s += time.monotonic() - t0
        total_b += src.stat().st_size
    return total_s / total_b, {"basis": "bytes", "bench_seconds": round(total_s, 2),
                               "bench_bytes": int(total_b)}


def _bench_similarity(tmp: Path):
    from polyphar.stages import similarity

    try:
        pm = similarity._pm_bin_dir(None)
    except FileNotFoundError as exc:
        raise Skip(str(exc)) from None
    pockets = []
    for p in sorted((BENCH_DIR / "sample_pockets").glob("*.pdb")):
        shutil.copyfile(p, tmp / p.name)
        pockets.append(tmp / p.name)
    n = len(pockets)
    pairs = n * (n - 1) // 2
    t0 = time.monotonic()
    similarity.run_pocketmatch(pockets, pm, tmp / "pm")
    seconds = time.monotonic() - t0
    return seconds / pairs, {"basis": "pairs", "bench_seconds": round(seconds, 2),
                             "bench_pairs": pairs}


# Uni-Dock seconds per ligand, by GPU device class. Measured 0.176 s/ligand on an L4 and
# -- the real simplification over CPU smina -- FLAT across box volumes from 3,679 to
# 20,112 A^3, where smina's cost was dominated by box volume. So n_ligands x rate(device)
# is an adequate estimate and no per-box term is needed. The other entries are scaled from
# the L4 figure by relative throughput and are priors, not measurements; a device not
# listed uses UNIDOCK_RATE_DEFAULT and says so.
UNIDOCK_RATE_S_PER_LIGAND = {"L4": 0.176, "T4": 0.42, "A100": 0.075, "H100": 0.05,
                             "V100": 0.30, "L40S": 0.11, "A10G": 0.20, "RTX 4090": 0.09}
UNIDOCK_RATE_DEFAULT = 0.176   # the measured L4 rate, used for an unlisted device
UNIDOCK_RATE_DEFAULT_CLASS = "unlisted-at-L4-rate"
UNIDOCK_RATE_SOURCE = ("Uni-Dock 1.2.0 on an L4, 2026-11-03: 0.176 s/ligand, flat across "
                       "box volumes 3,679-20,112 A^3")


def unidock_rate(model: str | None) -> tuple[float, str]:
    """(seconds per ligand, which device class it came from) for a GPU model name."""
    name = (model or "").strip()
    for key, rate in sorted(UNIDOCK_RATE_S_PER_LIGAND.items(), key=lambda kv: -len(kv[0])):
        if key.lower() in name.lower():
            return rate, key
    # No spaces or parens: this becomes a calibration.json key ("unidock:<device class>").
    return UNIDOCK_RATE_DEFAULT, UNIDOCK_RATE_DEFAULT_CLASS


UNIDOCK_TOOL_PREFIX = "unidock:"


def unidock_tool_prefix() -> str:
    return UNIDOCK_TOOL_PREFIX


def _bench_unidock():
    """Calibrate Uni-Dock from the per-device table rather than by running it.

    There is nothing a local run would add: the rate is flat across box volumes, so one
    number per device class is the whole model, and a benchmark would need a GPU and the
    real engine on the machine doing the estimating."""
    from polyphar.stages import dock_unidock as unidock

    model = unidock.gpu_model()
    if model is None:
        raise Skip("no CUDA GPU on this host")
    rate, device = unidock_rate(model)
    return rate, {"basis": "ligands", "tool": f"{unidock_tool_prefix()}{device}",
                  "gpu_model": model, "device_class": device,
                  "measured": False, "rate_source": UNIDOCK_RATE_SOURCE}


def _bench_dock(tmp: Path):
    from polyphar.stages import dock

    if bench_engine() == dock.ENGINE_UNIDOCK:
        return _bench_unidock()
    smina = dock.smina_bin()
    if shutil.which(smina) is None and not os.access(smina, os.X_OK):
        raise Skip(f"smina not found ({smina})")
    pdbqt = tmp / "receptor.pdbqt"
    dock.receptor_to_pdbqt(BENCH_DIR / DOCK_RECEPTOR, pdbqt)
    lig = tmp / "ligands.sdf"
    dock.write_titled_library(BENCH_DIR / DOCK_LIGANDS, lig)
    n = len(dock._records((BENCH_DIR / DOCK_LIGANDS).read_text()))
    center, size = dock.pocket_box(BENCH_DIR / DOCK_POCKET)
    t0 = time.monotonic()
    dock.run_smina(pdbqt, lig, center, size, tmp / "out.sdf", tmp / "out.log", seed=DOCK_SEED,
                   cpu=1, exhaustiveness=dock.EXHAUSTIVENESS, num_modes=dock.NUM_MODES)
    seconds = time.monotonic() - t0
    measured = seconds / n
    return measured * DOCK_MIX_FACTOR, {
        "basis": "ligands", "tool": dock.ENGINE_SMINA,
        "bench_seconds": round(seconds, 2), "bench_ligands": n,
        "measured_seconds_per_ligand": measured, "library_mix_factor": DOCK_MIX_FACTOR,
        "library_mix_source": DOCK_MIX_SOURCE}


def bench_engine() -> str:
    """Which dock engine to benchmark. $POLYPHAR_DOCK_ENGINE wins; otherwise Uni-Dock only
    when this host actually has BOTH a GPU and a runnable Uni-Dock, so a GPU machine doing
    a CPU run is never calibrated with the wrong engine."""
    from polyphar.stages import dock
    from polyphar.stages import dock_unidock as unidock

    chosen = os.environ.get("POLYPHAR_DOCK_ENGINE")
    if chosen in dock.ENGINES:
        return chosen
    if unidock.gpu_model() is not None and unidock.unidock_version() != "unavailable":
        return dock.ENGINE_UNIDOCK
    return dock.ENGINE_SMINA


def _live_session():
    from polyphar.common import http_session

    return http_session(retries=0)   # attempts and their backoff are handled by _timed_request


REQUEST_ATTEMPTS = 3
RETRY_PAUSE_S = 0.5


def _timed_request(fn, attempts: int = REQUEST_ATTEMPTS, pause: float | None = None):
    """Run fn() until it succeeds (at most `attempts` times, pausing between failures).
    Returns (seconds of the successful attempt only, number of failed attempts before it):
    a transient 5xx must not skip the stage, and neither its time nor the pause may inflate
    the calibration. Raises the last error when every attempt failed."""
    failed = 0
    while True:
        t0 = time.monotonic()
        try:
            fn()
        except Exception:
            failed += 1
            if failed >= attempts:
                raise
            time.sleep(RETRY_PAUSE_S if pause is None else pause)
        else:
            return time.monotonic() - t0, failed


def _bench_uniprot(tmp: Path):
    from polyphar.stages import uniprot

    ids = (BENCH_DIR / "accessions.txt").read_text().split()
    # fill UNIPROT_REQUESTS full batches; repeated accessions are fine for timing
    ids = (ids * UNIPROT_REQUESTS * UNIPROT_BATCH)[:UNIPROT_REQUESTS * UNIPROT_BATCH]
    session = _live_session()
    total, retries = 0.0, 0
    for i in range(UNIPROT_REQUESTS):
        chunk = ids[i * UNIPROT_BATCH:(i + 1) * UNIPROT_BATCH]
        try:
            seconds, failed = _timed_request(
                lambda: uniprot.fetch_uniprot(chunk, session, batch=UNIPROT_BATCH))
        except Exception as exc:
            raise Skip(f"UniProt not reachable: {str(exc)[:200]}") from None
        total += seconds
        retries += failed
    return total / UNIPROT_REQUESTS, {"basis": "requests", "bench_seconds": round(total, 2),
                                      "bench_requests": UNIPROT_REQUESTS, "retries": retries,
                                      "accessions_per_request": UNIPROT_BATCH}


def _bench_structures(tmp: Path):
    from polyphar.stages import structures

    session = _live_session()
    total, retries = 0.0, 0
    for acc, pdb, chain in STRUCTURE_JOBS:
        try:
            seconds, failed = _timed_request(
                lambda: structures._pdb_job(tmp, acc, pdb, chain, session))
        except Exception as exc:
            raise Skip(f"RCSB not reachable: {str(exc)[:200]}") from None
        total += seconds
        retries += failed
    n = len(STRUCTURE_JOBS)
    return total / n, {"basis": "downloads", "bench_seconds": round(total, 2),
                       "bench_requests": n, "retries": retries}


WORKLOADS = {"uniprot": _bench_uniprot, "structures": _bench_structures, "fix": _bench_fix,
             "similarity": _bench_similarity, "dock": _bench_dock}


def run_benchmark(stages, calib, log=print) -> dict:
    """Run the workloads for `stages`, fold results into `calib` (replacing any entry:
    a benchmark is a fresh measurement, not one more observation) and save it. Returns
    {stage: {"status": "calibrated" | "skipped" | "host-calibrated" | "failed", ...}}."""
    results = {}
    for stage in stages:
        if stage == "pockets":
            results[stage] = {"status": "host-calibrated", "reason": HOST_CALIBRATED}
            log(f"benchmark {stage}: {HOST_CALIBRATED}")
            continue
        workload = WORKLOADS.get(stage)
        if workload is None:
            results[stage] = {"status": "skipped", "reason": "trivial stage, nothing to calibrate"}
            continue
        t0 = time.monotonic()
        with tempfile.TemporaryDirectory(prefix="polyphar-bench-") as tmp:
            try:
                spc, extra = workload(Path(tmp))
            except Skip as exc:
                results[stage] = {"status": "skipped", "reason": str(exc)}
                log(f"benchmark {stage}: skipped: {exc}")
                continue
            except Exception as exc:   # one broken workload never stops the others
                results[stage] = {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
                log(f"benchmark {stage}: FAILED: {type(exc).__name__}: {exc}")
                continue
        elapsed = time.monotonic() - t0
        calib.data.pop(stage, None)
        # dock is calibrated PER ENGINE (and, for Uni-Dock, per GPU device class): an L4
        # rate must never be reused for an A100, nor a CPU smina rate for either.
        calib.update(stage, spc, tool=extra.pop("tool", None))
        if stage in calib.data:
            calib.data[stage].update({**extra, "source": "polyphar benchmark"})
        results[stage] = {"status": "calibrated", "seconds_per_unit": spc,
                          "elapsed_s": elapsed, **extra}
        log(f"benchmark {stage}: {spc:.4g} s per {extra.get('basis', 'unit')} "
            f"(one process, {elapsed:.0f} s)")
    calib.save()
    return results


# --- run state: observed values for the projection ---------------------------------------

# `similarity` runs AFTER `select`: it no longer feeds target selection (which is now
# panel filters alone) and is only needed for the post hoc similarity test, so it runs on
# the panel's eligible pockets. On net1 that is ~179k pairs against 67.2M, three orders
# of magnitude cheaper. The stage directory numbers are left as they are — renaming them
# would invalidate every existing run folder for no gain.
ORDER = ("ids", "uniprot", "structures", "fix", "pockets", "annotate", "select",
         "similarity", "ligands", "dock")
STAGE_DIRS = {"ids": "00_ids", "uniprot": "01_uniprot", "structures": "02_structures",
              "fix": "03_pdbfixer", "pockets": "04_pockets",
              "annotate": "04_pockets/annotation", "select": "06_targets",
              "similarity": "05_similarity", "ligands": "07_docking/library",
              "dock": "07_docking"}


def completed(layout) -> list[str]:
    """Stages whose directory has a manifest.json (written last by every stage)."""
    return [s for s in ORDER if (layout.run_dir / STAGE_DIRS[s] / "manifest.json").is_file()]


def _manifest_counts(layout, stage) -> dict:
    import json

    try:
        data = json.loads((layout.run_dir / STAGE_DIRS[stage] / "manifest.json").read_text())
        return data.get("counts") or {}
    except (OSError, ValueError, AttributeError):
        return {}


def _csv_rows(path) -> list[dict]:
    import csv

    try:
        with open(path, newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh))
    except OSError:
        return []


def _mean_size(paths) -> float | None:
    sizes = []
    for p in paths:
        try:
            sizes.append(p.stat().st_size)
        except OSError:
            pass
    return sum(sizes) / len(sizes) if sizes else None


def _count_records(path: Path) -> int | None:
    try:
        with open(path, "rb") as fh:
            n = sum(1 for line in fh if line.rstrip(b"\r\n") == b"$$$$")
        return n or None
    except OSError:
        return None


def run_state(layout, done: list[str], tool: str | None = None) -> tuple[dict, list[str]]:
    """Projection inputs: PRIORS (kept pockets per model by pocket `tool`, when given) with
    every value the finished stages have observed replacing its prior. Returns (state,
    names of the observed values)."""
    s = dict(PRIORS)
    if tool in KEPT_POCKETS_PRIOR:
        s["kept_pockets_per_model"] = KEPT_POCKETS_PRIOR[tool]
    observed = []

    def put(key, value):
        if value is not None and value == value and value > 0:
            s[key] = value
            observed.append(key)

    # accessions: ids.txt is exact once `ids` ran; uniprot.csv (kept records) once `uniprot` did
    n_ids = None
    try:
        n_ids = sum(1 for line in layout.ids_txt.read_text().split() if line) or None
    except OSError:
        pass
    put("accessions", n_ids)
    if "uniprot" in done:
        put("accessions", len(_csv_rows(layout.uniprot_csv)) or None)
    n_acc = s["accessions"]
    if "structures" in done:
        rows = _csv_rows(layout.protein_info)
        is_af = [str(r.get("AlphaFold", "")).strip().lower() == "true" for r in rows]
        pdb = [layout.run_dir / r["File"] for r, a in zip(rows, is_af) if not a and r.get("File")]
        af = [layout.run_dir / r["File"] for r, a in zip(rows, is_af) if a and r.get("File")]
        if rows:
            put("pdb_chains_per_accession", len(pdb) / n_acc if pdb else None)
            put("af_models_per_accession", len(af) / n_acc if af else None)
            put("pdb_chain_bytes", _mean_size(pdb))
            put("af_model_bytes", _mean_size(af))
            if not pdb:
                s["pdb_chains_per_accession"] = 0.0
                if "pdb_chains_per_accession" not in observed:
                    observed.append("pdb_chains_per_accession")
    if "fix" in done:
        # Chains per accession are NOT re-counted here: Protein_info.csv (structures) is the
        # authority on how many PDB chains exist; fix may drop some, which only shows up in
        # the pockets-stage size, so only the mean size of the fixed chains is observed.
        put("fixed_chain_bytes", _mean_size(layout.pdbfixer.glob("*/*.pdb")))
    if "pockets" in done:
        models = _manifest_counts(layout, "pockets").get("models")
        rows = _csv_rows(layout.drugscore_csv)
        if rows and "Selected" in rows[0]:   # the pockets stage's per-tool selection
            kept = sum(1 for r in rows if str(r.get("Selected")).strip().lower() == "true")
        else:                                # legacy DrugScore.csv
            kept = sum(1 for r in rows if _float(r.get("drugScore")) >= 0.5)
        if isinstance(models, int) and models > 0:
            put("kept_pockets_per_model", kept / models if kept else None)
    if "select" in done:
        put("targets", len(_csv_rows(layout.final_targets_csv)) or None)
    put("ligands", _count_records(layout.inputs / "e-Drug3D_2083.sdf"))
    if "ligands" in done:   # exact: the records the stage kept are the ones dock docks
        put("ligands_kept", _manifest_counts(layout, "ligands").get("kept"))
    return s, observed


def _float(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("-inf")


def projected_units(stage: str, s: dict) -> WorkUnits | None:
    """Work a later stage will have, from the (observed or prior) state. None for a
    trivial stage."""
    n_acc = s["accessions"]
    n_pdb = max(0, round(n_acc * s["pdb_chains_per_accession"]))
    n_af = max(0, round(n_acc * s["af_models_per_accession"]))
    if stage in ("ids", "select"):
        return None
    if stage == "uniprot":
        from polyphar.stages.uniprot import BATCH

        return WorkUnits(math.ceil(n_acc / BATCH), "requests")
    if stage == "structures":
        return WorkUnits(n_pdb + n_af, "downloads")
    if stage == "fix":
        return WorkUnits(n_pdb, "bytes", [float(s["pdb_chain_bytes"])] * n_pdb)
    if stage == "pockets":
        sizes = [float(s["fixed_chain_bytes"])] * n_pdb + [float(s["af_model_bytes"])] * n_af
        return WorkUnits(len(sizes), "bytes", sizes)
    if stage == "annotate":
        return WorkUnits(round((n_pdb + n_af) * s["kept_pockets_per_model"]), "pockets")
    if stage == "similarity":
        # pairs among the PANEL's eligible pockets, not every selected pocket: the stage
        # now runs after `select`. On net1 that is ~179k pairs rather than 67.2M.
        pockets = round((n_pdb + n_af) * s["kept_pockets_per_model"]
                        * s.get("eligible_fraction", ELIGIBLE_FRACTION))
        return WorkUnits(pockets * (pockets - 1) // 2, "pairs")
    if stage == "ligands":
        # the source library's records; s["ligands"] is the record count of the library in
        # <run>/inputs when it is there, else the prior (the kept count, close enough for
        # a stage that is seconds of work).
        return WorkUnits(round(s["ligands"]), "ligands")
    if stage == "dock":
        # the ligands actually docked are the ones `ligands` KEPT, not the source
        # library's records: the dockability filter drops ~6% of e-Drug3D_2083, so using
        # the raw count over-projects dock by that much.
        return WorkUnits(round(s["targets"] * s.get("ligands_kept", s["ligands"])),
                         "ligands")
    raise KeyError(stage)


# --- estimate ------------------------------------------------------------------------------

def _stage_args(stage: str, layout, nprocs: int, extra: list[str] | None = None):
    """Default arguments of `stage` (plus `extra`). --dry-run is set so units() never records
    anything (e.g. similarity's name guard would otherwise write failed.tsv): estimate is
    read-only."""
    from polyphar.cli import build_parser

    return build_parser().parse_args([stage, "--run-dir", str(layout.run_dir),
                                      "--nprocs", str(nprocs), "--dry-run", *(extra or [])])


def _pocket_tool(layout) -> str:
    """The run folder's recorded pocket tool, else the default ($POLYPHAR_POCKET_TOOL)."""
    from polyphar.stages import pockets

    try:
        return pockets.recorded_tool(layout) or pockets.default_tool()
    except (RuntimeError, SystemExit):
        return pockets.DEFAULT_TOOL


def estimate_rows(layout, nprocs: int, calib) -> tuple[list[dict], list[str]]:
    """One row per stage after the last completed one: {stage, kind, units, basis, seconds,
    basis_note}. kind is "exact" (the next stage, from its units()), "projected" or
    "trivial". Returns (rows, observed value names)."""
    import importlib

    from polyphar.cli import STAGES
    from polyphar.eta import PARALLEL_EFFICIENCY, estimate_stage

    done = completed(layout)
    last = max((ORDER.index(s) for s in done), default=-1)
    tool = _pocket_tool(layout)
    state, observed = run_state(layout, done, tool)
    rows = []
    for i, stage in enumerate(ORDER[last + 1:], start=last + 1):
        units, kind, note = None, "projected", ""
        trivial = projected_units(stage, state) is None
        if trivial:
            kind = "trivial"
        elif i == last + 1:
            try:
                module = importlib.import_module(STAGES[stage][0])
                extra = ["--pocket-tool", tool] if stage == "pockets" else None
                units = module.units(layout, _stage_args(stage, layout, nprocs, extra))
                kind = "exact"
            except Exception as exc:   # inputs not usable yet: fall back to the projection
                note = f"exact units unavailable ({type(exc).__name__}: {exc})"
        if units is None and not trivial:
            units = projected_units(stage, state)
        if trivial:
            rows.append({"stage": stage, "kind": kind, "units": "-", "seconds": 0.0, "note": ""})
            continue
        # dock's calibration is per engine, and per GPU device class for Uni-Dock: the
        # estimate must not quote a CPU smina rate for a GPU run (or the reverse). The
        # engine here is this host's (bench_engine); `polyphar dock --engine ...` uses the
        # flag, through cli.gate -> dock.calibration_tool.
        stage_tool, tpj, slots = None, 1, None
        if stage == "pockets":
            stage_tool = tool
        elif stage == "dock":
            from polyphar.stages import dock as _dock

            engine = bench_engine()
            eargs = SimpleNamespace(engine=engine, nprocs=nprocs)
            stage_tool = _dock.calibration_tool(eargs)
            # the stage's own slot count (1 for Uni-Dock: it saturates the one device), the
            # same source cli.gate uses -- not the engine-blind PARALLEL_EFFICIENCY table.
            slots = _dock.slots(eargs, tpj)
        est = estimate_stage(stage, units, nprocs, calib, PARALLEL_EFFICIENCY,
                             threads_per_job=tpj, tool=stage_tool, slots=slots)
        rows.append({"stage": stage, "kind": kind, "units": f"{units.count} ({units.cost_basis})",
                     "seconds": est.seconds, "note": note or est.basis})
    return rows, observed


def format_table(rows, observed, now: datetime | None = None) -> str:
    from polyphar.eta import UNCALIBRATED, Estimate, format_duration, format_eta

    now = now or datetime.now()
    if not rows:
        return "all stages are complete: nothing to estimate"
    out = [f"{'stage':<12}{'basis':<11}{'units':<22}{'estimate':<14}note"]
    known, missing = 0.0, []
    for r in rows:
        if r["seconds"] is None:
            eta = "uncalibrated"
            missing.append(r["stage"])
        else:
            known += r["seconds"]
            eta = "-" if r["kind"] == "trivial" else (
                "instant" if r["seconds"] == 0 else format_duration(r["seconds"]))
        out.append(f"{r['stage']:<12}{r['kind']:<11}{r['units']:<22}{eta:<14}{r['note']}")
    total = format_eta(Estimate(known, ""), now) if known else "instant"
    if missing:
        total = (f"at least {total}; uncalibrated, not counted: {', '.join(missing)} "
                 f"({UNCALIBRATED})")
    out.append(f"total: {total}")
    proj = [r["stage"] for r in rows if r["kind"] == "projected"]
    if proj:
        used = (f"observed in this run: {', '.join(observed)}; " if observed else "") + \
            f"the rest from priors ({PRIORS_SOURCE})"
        out.append(f"projected = {', '.join(proj)}; {used}")
    return "\n".join(out)
