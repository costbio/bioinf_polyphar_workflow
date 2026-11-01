"""Worker functions for tests/test_pool.py (a light module: spawn children import it)."""
import os
import time

from polyphar.pool import task

HEALTHY_S = 0.5


def work(i, crashers):
    def f():
        if i in crashers:
            time.sleep(0.05)
            os._exit(1)   # a hard worker death (BrokenProcessPool), like a segfault
        time.sleep(HEALTHY_S)
        return i
    return task(f)
