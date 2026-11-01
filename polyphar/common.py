"""Shared infrastructure: run-folder layout, work units, manifests, HTTP session."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

USER_AGENT = "polyphar (+https://github.com/costbio/bioinf_polyphar_workflow)"

# PocketMatch limits pocket names to 32 characters.
PM_NAME_LIMIT = 32


class RunLayout:
    """Every path in a run folder (one run folder = one subtype)."""

    def __init__(self, run_dir):
        self.run_dir = Path(run_dir)

    def stage_dir(self, stage: str) -> Path:
        return self.run_dir / stage

    @property
    def inputs(self) -> Path:
        return self.run_dir / "inputs"

    @property
    def ids_txt(self) -> Path:
        return self.run_dir / "00_ids" / "ids.txt"

    @property
    def uniprot_csv(self) -> Path:
        return self.run_dir / "01_uniprot" / "uniprot.csv"

    @property
    def structures(self) -> Path:
        return self.run_dir / "02_structures"

    @property
    def protein_info(self) -> Path:
        return self.structures / "Protein_info.csv"

    @property
    def pdbfixer(self) -> Path:
        return self.run_dir / "03_pdbfixer"

    @property
    def pockets(self) -> Path:
        return self.run_dir / "04_pockets"

    @property
    def drugscore_csv(self) -> Path:
        return self.pockets / "DrugScore.csv"

    @property
    def similarity(self) -> Path:
        return self.run_dir / "05_similarity"

    @property
    def pocketsim_csv(self) -> Path:
        return self.similarity / "Pocketsim.csv"

    @property
    def targets(self) -> Path:
        return self.run_dir / "06_targets"

    @property
    def final_targets_csv(self) -> Path:
        return self.targets / "final_targets.csv"

    @property
    def docking(self) -> Path:
        return self.run_dir / "07_docking"

    @property
    def logs(self) -> Path:
        return self.run_dir / "logs"


@dataclass(frozen=True)
class WorkUnits:
    """Exact work a stage will do. `sizes` are per-unit weights (e.g. file bytes)
    used for ETA and the longest-job bound."""

    count: int
    cost_basis: str
    sizes: list[float] | None = None


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def http_session(retries: int = 5, backoff: float = 2.0) -> requests.Session:
    """Session that retries on 429/5xx and identifies itself."""
    retry = Retry(
        total=retries,
        backoff_factor=backoff,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=None,  # retry all methods, including POST
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry)
    s = requests.Session()
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers["User-Agent"] = USER_AGENT
    return s


def append_failure(stage_dir, item, stage, reason) -> None:
    """Append one row to <stage_dir>/failed.tsv, writing the header once."""
    stage_dir = Path(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    path = stage_dir / "failed.tsv"
    clean = lambda x: str(x).replace("\t", " ").replace("\r", " ").replace("\n", " ")
    new = not path.exists() or path.stat().st_size == 0
    with open(path, "a", encoding="utf-8") as fh:
        if new:
            fh.write("item\tstage\treason\n")
        fh.write(f"{clean(item)}\t{clean(stage)}\t{clean(reason)}\n")


def _git_commit() -> str:
    """Code revision: $POLYPHAR_GIT_COMMIT (baked into the docker image, which
    has no .git) when set to a real value, else `git rev-parse HEAD`."""
    env = os.environ.get("POLYPHAR_GIT_COMMIT", "").strip()
    if env and env != "unknown":
        return env
    repo = Path(__file__).resolve().parent
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True, timeout=10
        )
        sha = out.stdout.strip()
        return sha if out.returncode == 0 and sha else "unknown"
    except Exception:
        return "unknown"


def write_manifest(stage_dir, stage, args: dict, inputs: list, versions: dict,
                   counts: dict, timing: dict, extra: dict | None = None) -> Path:
    """Write <stage_dir>/manifest.json recording inputs (with sha256), tool
    versions, counts, timing (estimated_s / actual_s) and the code's git commit.
    `extra` adds stage-specific top-level keys (it cannot replace the standard ones)."""
    stage_dir = Path(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "stage": stage,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "git_commit": _git_commit(),
        "python": sys.version.split()[0],
        "args": args,
        "inputs": [{"path": str(p), "sha256": sha256_file(p)} for p in inputs],
        "versions": versions,
        "counts": counts,
        "timing": timing,
    }
    for key, value in (extra or {}).items():
        if key in manifest:
            raise ValueError(f"manifest key {key!r} is reserved")
        manifest[key] = value
    path = stage_dir / "manifest.json"
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    tmp.replace(path)
    return path
