"""Docker image + host wrapper tests. Build the image first:

    docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) \
        -f docker/Dockerfile -t polyphar:dev .
"""
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.docker

REPO = Path(__file__).resolve().parent.parent
IMAGE = "polyphar:dev"
WRAPPER = REPO / "bin" / "polyphar-docker"


@pytest.fixture(scope="session")
def image():
    r = subprocess.run(["docker", "image", "inspect", IMAGE], capture_output=True)
    if r.returncode != 0:
        pytest.skip(f"image {IMAGE} not built (see module docstring for the build command)")
    return IMAGE


def _sh(image, script, **kw):
    return subprocess.run(
        ["docker", "run", "--rm", "--entrypoint", "bash", image, "-c", script],
        capture_output=True, text=True, **kw)


# The legacy upstream pipeline (Step0-cabbage.sh -type3 then Step3-PM_serial),
# run on a scratch copy so the image stays untouched.
UPSTREAM_PIPELINE = r"""
set -e
w=$(mktemp -d); cd "$w"
for f in $(ls /opt/pocketmatch/Sample_pockets); do
  /opt/pocketmatch/Step0-cabbage-core /opt/pocketmatch/Sample_pockets/$f >> outfile.cabbage
done
/opt/pocketmatch/Step0-END-FILE >> outfile.cabbage
/opt/pocketmatch/Step3-PM_serial outfile.cabbage > pm.log 2>&1
cat PocketMatch_score.txt
"""


def _scores(text):
    """pair key (name1, name2) -> (Pmin, Pmax) as the literal strings."""
    out = {}
    for line in text.splitlines():
        f = line.split()
        if len(f) >= 4:
            out[(f[0], f[1])] = (f[2], f[3])
    return out


def test_pocketmatch_reproduces_upstream_sample(image):
    run = _sh(image, UPSTREAM_PIPELINE)
    assert run.returncode == 0, run.stdout[-500:] + run.stderr
    ref_text = _sh(image, "cat /opt/pocketmatch/sample_outputs/PocketMatch_score.txt").stdout
    got, ref = _scores(run.stdout), _scores(ref_text)
    assert len(ref) == 120
    assert got.keys() == ref.keys()
    assert all(got[k] == ref[k] for k in ref), "Pmin/Pmax differ from upstream sample_outputs"


def test_wrapper_writes_into_local_folder(image, tmp_path):
    (tmp_path / "inputs").mkdir()
    (tmp_path / "inputs" / "ids.txt").write_text("P25815\nQ05586\n")
    env = {**os.environ, "POLYPHAR_RUN_DIR": str(tmp_path), "POLYPHAR_VERSION": "dev",
           "POLYPHAR_HOME": str(tmp_path / "home")}
    r = subprocess.run([str(WRAPPER), "ids", "--txt", "inputs/ids.txt", "--no-notify"],
                       cwd=tmp_path, env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    out = tmp_path / "00_ids" / "ids.txt"
    assert out.is_file()
    assert out.stat().st_uid == os.getuid()
    assert out.read_text().split() == ["P25815", "Q05586"]


def test_image_has_no_dogsite_or_credentials(image):
    r = _sh(image, "find / -xdev -iname 'dogsite*' -not -path '/proc/*' 2>/dev/null; "
                   "echo ---; env")
    assert r.returncode == 0, r.stderr
    found, env = r.stdout.split("---\n", 1)
    assert found.strip() == "", f"dogsite files in image: {found}"
    assert "PUSHOVER_" not in env

