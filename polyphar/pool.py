"""Shared spawn-process pool for the CPU stages (fix, similarity, dock).

run_pool(tasks, nprocs, on_result, on_failure) runs [(item, fn, args)] on spawn process
pools. When a worker dies hard (segfault, OOM kill: BrokenProcessPool), only the tasks that
were in flight go back to a fresh full-size pool; a task in flight during a second break
runs alone, so a crash is attributed to its own task and the other tasks keep running in
parallel. Pools that keep breaking without finishing anything, or solo runs that all die,
end the retries (the rest are failed, not retried). Any BaseException in the parent (Ctrl-C,
SIGTERM, SIGHUP) stops the workers and propagates.

Worker functions wrap their body in task(): a failure crosses the process boundary as
TaskError (a plain message, always pickles), and a worker being stopped by SIGTERM finishes
its cleanup and exits instead of going back to idle.

Heavy imports stay inside functions.
"""

from __future__ import annotations

import os
import signal
import time
from collections import deque
from dataclasses import dataclass

REASON_DIED = "worker process died (BrokenProcessPool)"
SOLO_RETRY_PROBES = 3
TERMINATE_GRACE_S = 15.0


class TaskError(RuntimeError):
    """A task's failure carried across the process boundary as a plain message."""


_BUSY = False   # worker process: inside a task (temp files / a child process to clean up)
_STOP = False   # worker process: SIGTERM received while busy


def stopping() -> bool:
    """True in a worker process that is being stopped (SIGTERM while busy)."""
    return _STOP


def worker_init() -> None:
    # SIGTERM (the parent stopping the pool): an idle worker exits at once. A busy one
    # unwinds through subprocess.run (which kills its child) and through its temp file
    # cleanup, then exits in task(). Further SIGTERMs (the pool's terminate_broken sends
    # one) do not cut that cleanup short.
    def handler(signum, frame):
        global _STOP
        if not _BUSY:
            os._exit(128 + signum)
        if not _STOP:
            _STOP = True
            raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, handler)


def task(fn):
    """Run a worker task: mark the process busy, carry failures across the process
    boundary as TaskError (always pickles), exit once stopped by SIGTERM."""
    global _BUSY
    _BUSY = True
    try:
        return fn()
    except TaskError:
        raise
    except Exception as exc:
        raise TaskError(f"{type(exc).__name__}: {exc}") from None
    except BaseException:
        os._exit(143)   # SIGTERM / Ctrl-C: the cleanup inside fn has run
    finally:
        _BUSY = False
        if _STOP:
            os._exit(143)   # the parent is stopping the pool; never go back to idle


def reason(exc) -> str:
    return str(exc) if isinstance(exc, TaskError) else f"{type(exc).__name__}: {exc}"


def terminate(ex) -> None:
    """Cancel pending tasks and stop running ones now (an interrupted run leaves no
    workers and no children behind)."""
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


def _fail(failures, on_failure, item, why) -> None:
    failures.append((item, why))
    if on_failure is not None:
        on_failure(item, why)


def _pool_pass(tasks, queue, nprocs: int, on_result, on_failure, failures):
    """Run tasks from `queue` (indices into tasks) on one fresh spawn pool, keeping at most
    nprocs in flight, until the queue is empty or the pool breaks. Returns (indices that
    were in flight when the pool broke, number of tasks completed, whether it broke)."""
    import multiprocessing
    from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
    from concurrent.futures.process import BrokenProcessPool

    ctx = multiprocessing.get_context("spawn")   # workers start clean (no fork of a threaded parent)
    completed = 0
    with ProcessPoolExecutor(max_workers=max(1, min(nprocs, len(queue))), mp_context=ctx,
                             initializer=worker_init) as ex:
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
                        return list(inflight.values()), completed, True
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
                        _fail(failures, on_failure, tasks[i][0], reason(exc))
                    else:
                        if on_result is not None:
                            on_result(tasks[i][0], result)
                    completed += 1
                if broke:   # everything still in flight went down with it
                    return broke + list(inflight.values()), completed, True
            return [], completed, False
        except BaseException:
            terminate(ex)
            raise


def run_pool(tasks, nprocs: int, on_result=None, on_failure=None) -> PoolOutcome:
    """Run tasks [(item, fn, args)] on spawn process pools; returns the failed items with
    their reasons. on_result(item, result) and on_failure(item, reason) run in this (parent)
    process, once per task. See the module docstring for the retry rules."""
    failures: list = []
    queue = deque(range(len(tasks)))
    breaks: dict[int, int] = {}
    solo: list[int] = []
    empty = total = 0
    while queue:
        broken, completed, broke = _pool_pass(tasks, queue, nprocs, on_result, on_failure,
                                              failures)
        total += completed
        if not broke:
            continue
        empty = 0 if completed else empty + 1
        retry = []
        for i in broken:
            breaks[i] = breaks.get(i, 0) + 1
            (solo if breaks[i] >= 2 else retry).append(i)
        queue.extendleft(reversed(sorted(retry)))
        if empty >= SOLO_RETRY_PROBES:
            if total == 0:   # pools cannot run anything at all
                for i in list(queue) + solo:
                    _fail(failures, on_failure, tasks[i][0],
                          REASON_DIED + "; not retried, every pool died")
                return PoolOutcome(failures, [])
            # pools keep breaking before anything finishes (e.g. at submit): the rest
            # run alone, so whatever breaks them is attributed to its own task
            solo += list(queue)
            queue.clear()
    isolated, tried, died = [], 0, 0
    for i in sorted(set(solo)):
        if tried >= SOLO_RETRY_PROBES and died == tried:
            _fail(failures, on_failure, tasks[i][0], REASON_DIED + "; not retried, every retry died")
            continue
        tried += 1
        isolated.append(tasks[i][0])
        broken, _, _ = _pool_pass(tasks, deque([i]), 1, on_result, on_failure, failures)
        if broken:
            _fail(failures, on_failure, tasks[i][0], REASON_DIED)
            died += 1
    return PoolOutcome(failures, isolated)
