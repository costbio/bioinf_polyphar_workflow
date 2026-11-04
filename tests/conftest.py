import os
import shutil
import subprocess

import pytest


def _rdkit_ok() -> bool:
    try:
        import rdkit  # noqa: F401
    except ImportError:
        return False
    return True


def _openbabel_ok() -> bool:
    try:
        from openbabel import pybel  # noqa: F401
    except Exception:
        return False
    return True


def _docker_ok() -> bool:
    if shutil.which("docker") is None:
        return False
    try:
        return subprocess.run(["docker", "info"], capture_output=True, timeout=15).returncode == 0
    except Exception:
        return False


# marker -> (predicate that is True when the test may run, skip reason)
_REQUIREMENTS = {
    "regression": (lambda: bool(os.environ.get("ADEP_DATA")), "env ADEP_DATA not set"),
    "dogsite": (lambda: bool(os.environ.get("DOGSITE_BIN")), "env DOGSITE_BIN not set"),
    "pocketmatch": (lambda: bool(os.environ.get("PM_BIN_DIR")), "env PM_BIN_DIR not set"),
    "docking": (lambda: shutil.which("smina") is not None, "smina not on PATH"),
    "docker": (_docker_ok, "docker CLI not available"),
    "rdkit": (_rdkit_ok, "rdkit not importable"),
    "openbabel": (_openbabel_ok, "openbabel not importable"),
    "pushover": (
        lambda: bool(os.environ.get("PUSHOVER_APP_TOKEN") and os.environ.get("PUSHOVER_USER_KEY")),
        "PUSHOVER_APP_TOKEN / PUSHOVER_USER_KEY not set",
    ),
}
# 'network' tests are not skipped: they run by default and only need internet.


def pytest_collection_modifyitems(config, items):
    cache = {}
    for item in items:
        for name, (ok, reason) in _REQUIREMENTS.items():
            if name in item.keywords:
                if name not in cache:
                    cache[name] = ok()
                if not cache[name]:
                    item.add_marker(pytest.mark.skip(reason=reason))
