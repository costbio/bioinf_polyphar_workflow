"""ETA engine: per-machine calibration, stage estimates, milestone progress."""

from __future__ import annotations

import json
import math
import os
import platform
import time
from collections import deque
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from polyphar.common import WorkUnits

MILESTONES = (0.05, 0.25, 0.5, 0.75)
LOG_INTERVAL_S = 60.0


# --- slots per stage -------------------------------------------------------
# stage -> callable(nprocs, threads_per_job) -> parallel slots. None = trivial
# stage (seconds of work at most; never worth calibrating or waiting for).
# This is the only place the per-stage parallelism model lives.
PARALLEL_EFFICIENCY = {
    "ids": None,
    "select": None,
    "uniprot": lambda nprocs, tpj=1: min(nprocs, 8),       # network bound
    "structures": lambda nprocs, tpj=1: min(nprocs, 8),    # network bound
    "fix": lambda nprocs, tpj=1: nprocs,
    "pockets": lambda nprocs, tpj=1: nprocs,
    "similarity": lambda nprocs, tpj=1: max(1, nprocs / 2),  # block-pair overhead
    # CPU smina: one unit per nprocs//threads-per-job. This entry is ENGINE-BLIND and is
    # therefore NOT what the dock stage's own ETA uses: Uni-Dock runs ONE unit at a time
    # (it batches the whole library onto the one device), so `dock.slots()` -- the stage's
    # own slot count, which `cli.gate` asks for and passes to both estimate_stage and
    # Progress -- returns 1 for it whatever --nprocs says. Deriving the slots from the
    # table here would divide a GPU ETA by --nprocs AND write back an --nprocs-inflated
    # seconds-per-ligand into calibration.json. This entry is the fallback for callers that
    # have no stage object (projected estimates for a stage not yet selected).
    "dock": lambda nprocs, tpj=1: max(1, nprocs // max(1, tpj)),
    "dock-worker": lambda nprocs, tpj=1: max(1, nprocs // max(1, tpj)),
    "dock-export": lambda nprocs, tpj=1: 1,
    "dock-merge": None,
    # reading the panel's result files is I/O bound; the permutation itself is one seeded
    # stream in this process and does not parallelise (it must not: --nprocs cannot change
    # a result).
    "stats": lambda nprocs, tpj=1: min(nprocs, 8),
}


def slots_for(stage: str, nprocs: int, threads_per_job: int = 1,
              table: dict | None = None):
    """Parallel slots for a stage, or None when the stage is trivial."""
    table = PARALLEL_EFFICIENCY if table is None else table
    fn = table.get(stage, lambda n, t=1: n)
    if fn is None:
        return None
    return max(1.0, float(fn(nprocs, threads_per_job)))


def stage_slots(name: str, stage, args, threads_per_job: int = 1):
    """The slot count to estimate AND calibrate `stage` with: its own `slots()` when it has
    one, else the PARALLEL_EFFICIENCY table. None for a trivial stage.

    A stage's parallelism can depend on more than nprocs (dock: the engine -- Uni-Dock runs
    one unit at a time whatever --nprocs says), and the table cannot see that."""
    fn = getattr(stage, "slots", None)
    if fn is None:
        return slots_for(name, getattr(args, "nprocs", 1), threads_per_job)
    return max(1.0, float(fn(args, threads_per_job)))


# --- calibration -----------------------------------------------------------

def polyphar_home() -> Path:
    return Path(os.environ.get("POLYPHAR_HOME") or (Path.home() / ".polyphar"))


def _cpu_model() -> str:
    try:
        with open("/proc/cpuinfo", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine() or "unknown"


# Stages whose calibration depends on an external tool, and the tool an entry without a
# `tool` field was measured with (pockets: DoGSiteScorer 2, the only tool before DoGSite3).
# dock: an entry written before the engine was selectable was measured with smina, the
# only engine then.
LEGACY_TOOLS = {"pockets": "dogsite2", "dock": "smina", "dock-worker": "smina"}

# How to re-measure a stage for a tool it has no calibration for. pockets is the exception:
# DoGSiteScorer runs on the host, not through `polyphar benchmark`.
RECALIBRATE_HINT = {"pockets": "bin/polyphar-dogsite-worker --benchmark --tool {tool}"}

# Stages that are calibrated under another stage's name. `dock-worker` runs exactly the
# engine calls `dock` does, on the same unit basis (ligands), so it uses dock's measured
# rate instead of printing "uncalibrated" on the very host where the engine runs.
# `dock-export` is deliberately NOT here: it copies files, so dock's per-ligand rate would
# be nonsense for it.
CALIBRATION_ALIAS = {"dock-worker": "dock"}


class Calibration:
    """Seconds per unit-basis on ONE process, per stage, in
    <home>/calibration.json. A missing or unreadable file is an empty
    calibration."""

    def __init__(self, path, data: dict | None = None):
        self.path = Path(path)
        self.data = data or {}

    @classmethod
    def load(cls, home) -> "Calibration":
        path = Path(home) / "calibration.json"
        data = {}
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
            for stage, entry in raw.items():
                c = entry.get("seconds_per_unit") if isinstance(entry, dict) else None
                if isinstance(c, (int, float)) and math.isfinite(c) and c > 0:
                    n = entry.get("n_observations")
                    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
                        entry = {**entry, "n_observations": 1}  # hand-edited / older file
                    data[stage] = entry
        except (OSError, ValueError, AttributeError):
            pass
        return cls(path, data)

    def tool_of(self, stage: str) -> str | None:
        """The tool an entry was measured with (pockets: an entry written before tools were
        recorded is DoGSiteScorer 2, the only tool then), or None without an entry."""
        entry = self.data.get(stage)
        if not entry:
            return None
        return entry.get("tool") or LEGACY_TOOLS.get(stage)

    def cost(self, stage: str, tool: str | None = None) -> float | None:
        """Seconds per unit; with `tool`, only an entry measured with that tool counts
        (another tool's calibration is never used silently)."""
        entry = self.data.get(stage)
        if not entry or (tool is not None and self.tool_of(stage) != tool):
            return None
        return float(entry["seconds_per_unit"])

    def update(self, stage: str, seconds_per_unit: float, tool: str | None = None) -> None:
        """Fold one measurement into the running mean. Non-positive or
        non-finite measurements are ignored. With `tool`, an entry of another tool is
        replaced (not averaged) and the entry records the tool."""
        if not (isinstance(seconds_per_unit, (int, float))
                and math.isfinite(seconds_per_unit) and seconds_per_unit > 0):
            return
        if tool is not None and self.tool_of(stage) not in (None, tool):
            self.data.pop(stage, None)
        old = self.data.get(stage)
        try:
            n = max(1, int(old["n_observations"])) if old else 0
        except (KeyError, TypeError, ValueError):
            n = 1
        mean = (float(old["seconds_per_unit"]) * n + seconds_per_unit) / (n + 1) if old \
            else float(seconds_per_unit)
        self.data[stage] = {
            "seconds_per_unit": mean,
            "n_observations": n + 1,
            "nproc": os.cpu_count(),
            "cpu_model": _cpu_model(),
            "updated": datetime.now().astimezone().isoformat(timespec="seconds"),
        }
        if tool is not None:
            self.data[stage]["tool"] = tool

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_name(f"{self.path.name}.{os.getpid()}.tmp")
        tmp.write_text(json.dumps(self.data, indent=2, sort_keys=True) + "\n",
                       encoding="utf-8")
        tmp.replace(self.path)


# --- estimates -------------------------------------------------------------

UNCALIBRATED = "uncalibrated — run polyphar benchmark"


@dataclass(frozen=True)
class Estimate:
    seconds: float | None
    basis: str


def estimate_stage(stage: str, units: WorkUnits, nprocs: int, calib: Calibration,
                   parallel_efficiency: dict | None = None, *,
                   threads_per_job: int = 1, tool: str | None = None,
                   slots: float | None = None) -> Estimate:
    """Wall-clock estimate: max(total work / slots, longest single job). With `tool`, only
    a calibration measured with that tool is used.

    `slots`, when given, OVERRIDES the PARALLEL_EFFICIENCY table. The table is engine-blind;
    a stage that knows better (dock: one unit at a time on a GPU) reports its own count and
    `cli.gate` passes it here and to Progress, so the ETA and the calibration write-back
    both use the number the stage will actually run at."""
    if slots is None:
        slots = slots_for(stage, nprocs, threads_per_job, parallel_efficiency)
    if slots is None:
        return Estimate(0.0, "trivial")
    if units.count == 0:
        return Estimate(0.0, "no work")
    key = CALIBRATION_ALIAS.get(stage, stage)
    c = calib.cost(key, tool)
    if c is None:
        other = calib.tool_of(key) if tool is not None else None
        if other is not None:
            hint = RECALIBRATE_HINT.get(stage, "polyphar benchmark --stages {stage}")
            return Estimate(None, f"uncalibrated for {tool} (calibration.json has {other}) — "
                                  f"run " + hint.format(tool=tool, stage=stage))
        return Estimate(None, UNCALIBRATED)
    basis = f"calibrated ({units.cost_basis}, {slots:g} slots)"
    if units.sizes:
        return Estimate(max(sum(units.sizes) * c / slots, max(units.sizes) * c), basis)
    return Estimate(units.count * c / slots, basis)


def format_duration(seconds: float) -> str:
    if seconds < 60:
        return f"~{int(round(seconds))} s"
    minutes = int(round(seconds / 60))
    hours, minutes = divmod(minutes, 60)
    if hours == 0:
        return f"~{minutes} min"
    return f"~{hours} h {minutes} min"


def format_eta(est: Estimate, now: datetime) -> str:
    """e.g. "~3 h 20 min (finish ≈ Sat 14:05)"."""
    if est.seconds is None:
        return est.basis
    if est.seconds == 0:
        return f"instant ({est.basis})"
    finish = now + timedelta(seconds=est.seconds)
    fmt = "%a %d %b %H:%M" if (finish.date() - now.date()).days >= 7 else "%a %H:%M"
    return f"{format_duration(est.seconds)} (finish ≈ {finish.strftime(fmt)})"


# --- progress --------------------------------------------------------------

class Progress:
    """Live progress for one stage run. The stage calls advance(); the CLI
    calls done() or fail(). Notifier failures never propagate."""

    def __init__(self, stage: str, total: int, est: Estimate, notifier,
                 milestones=MILESTONES, *, total_size: float | None = None,
                 calib: Calibration | None = None, slots: float = 1.0,
                 log_path=None, clock=time.monotonic, tool: str | None = None):
        self.stage = stage
        self.tool = tool
        self.total = total
        self.total_size = total_size
        self.est = est
        self.notifier = notifier
        self.milestones = sorted(milestones)
        self.calib = calib
        self.slots = slots
        self.log_path = Path(log_path) if log_path else None
        self.clock = clock
        self.t0 = clock()
        self.units_done = 0
        self.work_done = 0.0
        self._next_ms = 0
        self._last_log = None
        self._recent = deque(maxlen=10)
        self._log(f"start: {total} units, estimate {format_eta(est, datetime.now())}")

    # -- internals --
    @property
    def _work_total(self) -> float:
        return float(self.total_size if self.total_size is not None else self.total)

    def _log(self, msg: str) -> None:
        line = f"{datetime.now().strftime('%H:%M:%S')} {self.stage}: {msg}"
        self._recent.append(line)
        print(line, flush=True)
        if self.log_path is not None:
            try:
                self.log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(self.log_path, "a", encoding="utf-8") as fh:
                    fh.write(line + "\n")
            except OSError:
                pass

    def _notify(self, title: str, message: str, **kw) -> None:
        try:
            self.notifier(title, message, **kw)
        except Exception:
            pass

    def _live_eta(self, fraction: float) -> str:
        remaining = (self.clock() - self.t0) * (1 - fraction) / fraction
        return format_eta(Estimate(remaining, "live rate"), datetime.now())

    # -- API --
    def advance(self, n: int = 1, size: float | None = None) -> None:
        self.units_done += n
        if self.total_size is None:
            self.work_done += n
        elif size is not None:
            self.work_done += size
        elif self.total > 0:   # no size given: assume the mean unit size
            self.work_done += n * self.total_size / self.total
        total = self._work_total
        if total <= 0:
            return
        fraction = self.work_done / total
        now = self.clock()
        if self._last_log is None or now - self._last_log >= LOG_INTERVAL_S:
            self._last_log = now
            self._log(f"{self.units_done}/{self.total} units ({fraction:.0%})")
        while (self._next_ms < len(self.milestones)
               and self.work_done >= self.milestones[self._next_ms] * total - 1e-9):
            m = self.milestones[self._next_ms]
            self._next_ms += 1
            eta = self._live_eta(max(fraction, 1e-9))
            label = "revised ETA" if m == self.milestones[0] else "ETA"
            msg = f"{m:.0%} done after {format_duration(now - self.t0)}; {label}: {eta}"
            self._log(msg)
            self._notify(f"{self.stage} {m:.0%}", msg)

    def done(self, counts: dict | None = None) -> None:
        """Post-stage bookkeeping. The stage already succeeded, so nothing here
        may turn the run into a failure."""
        try:
            self._done(counts or {})
        except Exception as exc:
            try:
                self._log(f"warning: post-run bookkeeping failed: "
                          f"{type(exc).__name__}: {exc}")
            except Exception:
                pass

    def _done(self, counts: dict) -> None:
        elapsed = self.clock() - self.t0
        failed = counts.get("failed", counts.get("failures", 0))
        est = ("n/a" if self.est.seconds is None
               else format_duration(self.est.seconds))
        detail = ", ".join(f"{k}={v}" for k, v in counts.items())
        msg = (f"done in {format_duration(elapsed)} (estimated {est}); "
               f"counts: {detail or 'none'}; failures: {failed}")
        self._log(msg)
        flagged = any(isinstance(counts.get(k), (int, float)) and not isinstance(counts.get(k), bool)
                      and counts.get(k) > 0 for k in ("failed", "incomplete_targets"))
        self._notify(f"{self.stage} done", msg, priority=1 if flagged else 0)
        if self.calib is not None and self.work_done > 0 and elapsed > 0:
            # seconds per unit on ONE process: elapsed wall time * slots / work
            self.calib.update(self.stage, elapsed * self.slots / self.work_done, self.tool)
            try:
                self.calib.save()
            except OSError as exc:
                self._log(f"could not save calibration: {exc}")

    def fail(self, exc: BaseException) -> None:
        self._log(f"FAILED: {type(exc).__name__}: {exc}")
        self._notify(f"{self.stage} FAILED",
                     "\n".join(self._recent), priority=1)
