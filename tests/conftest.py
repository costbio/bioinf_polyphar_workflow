import os
import shutil
import subprocess
import warnings
from pathlib import Path

import pytest

# Import prody HERE, once, with warnings silenced.
#
# `pytest -W error` is the command the README and the runbook tell you to run, and it had
# stopped working: prody's import emits a pyparsing deprecation
# ("'setParseAction' deprecated"), so under -W error every test that reached prody failed
# as an error -- 20 of them, in a suite that otherwise passed. The quality gate meant to
# catch OUR warnings was failing on a third party's, which is the same shape of problem as
# a test nobody runs. Importing it once at collection, outside any test's warning context,
# makes the module cached and the deprecation harmless, and leaves -W error free to fail
# on warnings this project actually raises. The ini filterwarnings entry cannot do this:
# a command-line -W takes precedence over it.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        import prody  # noqa: F401
    except Exception:
        pass


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


def _repo_files_ok() -> bool:
    """The image ships polyphar/ and tests/ only (docker/Dockerfile), not docs/ or bin/.
    Tests that assert on a document or on the host wrapper are repo-hygiene checks: they
    belong to a checkout, not to the container, and must skip rather than fail there."""
    root = Path(__file__).resolve().parents[1]
    return (root / "docs").is_dir() and (root / "bin").is_dir()


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
    "repo_files": (_repo_files_ok, "docs/ and bin/ are not in the image, only a checkout"),
    # pgrep is a test-only tool here: polyphar kills a timed-out engine by process group,
    # never by name, so the slim image has no procps and these assertions cannot run there.
    "procps": (lambda: shutil.which("pgrep") is not None, "pgrep not available (no procps)"),
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
