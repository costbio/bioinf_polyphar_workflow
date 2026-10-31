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
    "dock": lambda nprocs, tpj=1: max(1, nprocs // max(1, tpj)),
}


def slots_for(stage: str, nprocs: int, threads_per_job: int = 1,
              table: dict | None = None):
    """Parallel slots for a stage, or None when the stage is trivial."""
    table = PARALLEL_EFFICIENCY if table is None else table
    fn = table.get(stage, lambda n, t=1: n)
    if fn is None:
        return None
    return max(1.0, float(fn(nprocs, threads_per_job)))


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

    def cost(self, stage: str) -> float | None:
        entry = self.data.get(stage)
        return float(entry["seconds_per_unit"]) if entry else None

    def update(self, stage: str, seconds_per_unit: float) -> None:
        """Fold one measurement into the running mean. Non-positive or
        non-finite measurements are ignored."""
        if not (isinstance(seconds_per_unit, (int, float))
                and math.isfinite(seconds_per_unit) and seconds_per_unit > 0):
            return
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
                   threads_per_job: int = 1) -> Estimate:
    """Wall-clock estimate: max(total work / slots, longest single job)."""
    slots = slots_for(stage, nprocs, threads_per_job, parallel_efficiency)
    if slots is None:
        return Estimate(0.0, "trivial")
    if units.count == 0:
        return Estimate(0.0, "no work")
    c = calib.cost(stage)
    if c is None:
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
                 log_path=None, clock=time.monotonic):
        self.stage = stage
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
        self._notify(f"{self.stage} done", msg)
        if self.calib is not None and self.work_done > 0 and elapsed > 0:
            # seconds per unit on ONE process: elapsed wall time * slots / work
            self.calib.update(self.stage, elapsed * self.slots / self.work_done)
            try:
                self.calib.save()
            except OSError as exc:
                self._log(f"could not save calibration: {exc}")

    def fail(self, exc: BaseException) -> None:
        self._log(f"FAILED: {type(exc).__name__}: {exc}")
        self._notify(f"{self.stage} FAILED",
                     "\n".join(self._recent), priority=1)
