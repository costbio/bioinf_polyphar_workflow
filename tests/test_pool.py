"""Shared spawn pool: crashes are attributed to the crashing task only."""
import concurrent.futures as cf
import multiprocessing
import time
from concurrent.futures.process import BrokenProcessPool

import pytest

import _pool_tasks
from polyphar import pool

N, NPROCS = 28, 4


def _run(crashers):
    got = []
    tasks = [(i, _pool_tasks.work, (i, frozenset(crashers))) for i in range(N)]
    t0 = time.monotonic()
    out = pool.run_pool(tasks, NPROCS, lambda item, r: got.append((item, r)))
    return out, sorted(got), time.monotonic() - t0


@pytest.mark.parametrize("crashers", [{1, 5}, {9, 13}, {0, 1, 2, 3}, {0, 1, 2, 3, 4, 5},
                                      {4, 5, 6, 7, 9}],
                         ids=["early", "after-progress", "first-nprocs", "first-nprocs+2",
                              "after-solo-deaths"])
def test_two_crashers_fail_alone_and_healthy_tasks_stay_parallel(crashers):
    out, got, wall = _run(crashers)
    healthy = [i for i in range(N) if i not in crashers]
    assert sorted(out.failures) == [(i, pool.REASON_DIED) for i in sorted(crashers)]
    assert got == [(i, i) for i in healthy]           # every healthy task succeeded, once
    # only the crashers and the tasks that shared both of a crasher's broken pools run
    # alone; every other healthy task ran in a pool of size nprocs
    assert crashers <= set(out.isolated)
    assert len(out.isolated) <= len(crashers) * NPROCS
    pooled = [i for i in healthy if i not in out.isolated]
    assert len(pooled) >= len(healthy) - len(crashers) * (NPROCS - 1)
    serial = len(healthy) * _pool_tasks.HEALTHY_S
    assert wall < serial / 2, wall                    # close to parallel, not solo serial
    assert not multiprocessing.active_children()


def test_every_task_crashing_stops_after_a_few_probes():
    """Every task crashes: the first pool's tasks die alone, then nprocs tasks no pool broke
    under (fresh probes) die alone too, and the rest are failed without running. When only
    the first few tasks crash, a fresh probe succeeds instead (see above)."""
    out, got, _ = _run(set(range(N)))
    assert got == [] and sorted(i for i, _ in out.failures) == list(range(N))
    assert sorted(out.isolated) == list(range(2 * NPROCS))   # nprocs + nprocs fresh probes
    reasons = dict(out.failures)
    assert all(reasons[i] == pool.REASON_DIED for i in range(2 * NPROCS))
    assert all(reasons[i] == pool.REASON_DIED + "; not retried, every pool died"
               for i in range(2 * NPROCS, N))
    assert not multiprocessing.active_children()


class _SubmitBreaks(cf.ProcessPoolExecutor):
    """Multi-worker pools break at submit (nothing in flight) from the `after`-th submit."""
    after = 0
    n = 0

    def submit(self, *a, **kw):
        if self._max_workers > 1:
            type(self).n += 1
            if type(self).n > type(self).after:
                raise BrokenProcessPool("simulated")
        return super().submit(*a, **kw)


@pytest.mark.parametrize("after", [0, 3])
def test_pools_breaking_at_submit_route_the_rest_to_solo(monkeypatch, after):
    """A pool that breaks with nothing in flight attributes the break to no task: after a
    few such empty pools the remaining tasks run alone instead of being failed."""
    monkeypatch.setattr(_SubmitBreaks, "after", after)
    monkeypatch.setattr(_SubmitBreaks, "n", 0)
    monkeypatch.setattr(cf, "ProcessPoolExecutor", _SubmitBreaks)
    got = []
    tasks = [(i, _pool_tasks.work, (i, frozenset())) for i in range(6)]
    out = pool.run_pool(tasks, 2, lambda item, r: got.append(item))
    assert out.failures == [] and sorted(got) == list(range(6))
    assert _SubmitBreaks.n > after                  # the break really happened
    assert sorted(out.isolated) == list(range(after, 6))   # the units never submitted
    assert not multiprocessing.active_children()
