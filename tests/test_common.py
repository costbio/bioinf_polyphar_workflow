import json
import subprocess
import sys
from pathlib import Path

import pytest

from polyphar.common import (
    PM_NAME_LIMIT,
    RunLayout,
    WorkUnits,
    append_failure,
    sha256_file,
    write_manifest,
)

SUBCOMMANDS = [
    "ids", "uniprot", "structures", "fix", "pockets", "similarity",
    "select", "dock", "run", "estimate", "benchmark",
]


def test_layout_paths(tmp_path):
    lay = RunLayout(tmp_path)
    assert lay.pocketsim_csv == tmp_path / "05_similarity/Pocketsim.csv"
    assert lay.inputs == tmp_path / "inputs"
    assert lay.ids_txt == tmp_path / "00_ids/ids.txt"
    assert lay.uniprot_csv == tmp_path / "01_uniprot/uniprot.csv"
    assert lay.structures == tmp_path / "02_structures"
    assert lay.protein_info == tmp_path / "02_structures/Protein_info.csv"
    assert lay.pdbfixer == tmp_path / "03_pdbfixer"
    assert lay.pockets == tmp_path / "04_pockets"
    assert lay.drugscore_csv == tmp_path / "04_pockets/DrugScore.csv"
    assert lay.similarity == tmp_path / "05_similarity"
    assert lay.targets == tmp_path / "06_targets"
    assert lay.final_targets_csv == tmp_path / "06_targets/final_targets.csv"
    assert lay.docking == tmp_path / "07_docking"
    assert lay.logs == tmp_path / "logs"
    assert lay.stage_dir("01_uniprot") == tmp_path / "01_uniprot"


def test_pm_name_limit():
    assert PM_NAME_LIMIT == 32


def test_workunits_frozen():
    wu = WorkUnits(count=3, cost_basis="files", sizes=[1.0, 2.0, 3.0])
    assert wu.count == 3 and wu.sizes == [1.0, 2.0, 3.0]
    with pytest.raises(Exception):
        wu.count = 4


def test_sha256_file(tmp_path):
    p = tmp_path / "a.txt"
    p.write_bytes(b"abc")
    assert sha256_file(p) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_manifest_records_inputs_versions_timing(tmp_path):
    inp = tmp_path / "in.txt"
    inp.write_text("hello")
    sd = tmp_path / "stage"
    path = write_manifest(
        sd, "uniprot", {"nprocs": 2}, [inp], {"python": "3.12"},
        {"done": 5}, {"estimated_s": 10.0, "actual_s": 9.5},
    )
    assert path == sd / "manifest.json"
    m = json.loads(path.read_text())
    assert m["stage"] == "uniprot"
    assert m["args"] == {"nprocs": 2}
    assert m["inputs"][0]["sha256"] == sha256_file(inp)
    assert m["inputs"][0]["path"] == str(inp)
    assert m["versions"] == {"python": "3.12"}
    assert m["counts"] == {"done": 5}
    assert m["timing"]["estimated_s"] == 10.0
    assert m["timing"]["actual_s"] == 9.5
    assert m["git_commit"]
    assert "timestamp" in m


def test_append_failure_header_once(tmp_path):
    append_failure(tmp_path, "P1", "uniprot", "not found")
    append_failure(tmp_path, "P2", "uniprot", "timeout\twith tab")
    lines = (tmp_path / "failed.tsv").read_text().splitlines()
    assert lines[0] == "item\tstage\treason"
    assert sum(1 for l in lines if l.startswith("item\t")) == 1
    assert lines[1] == "P1\tuniprot\tnot found"
    assert len(lines) == 3 and lines[2].count("\t") == 2


def test_http_session_retries_and_ua():
    from polyphar.common import http_session

    s = http_session(retries=3, backoff=1.5)
    assert s.headers["User-Agent"] == (
        "polyphar (+https://github.com/costbio/bioinf_polyphar_workflow)"
    )
    retry = s.get_adapter("https://example.org").max_retries
    assert retry.total == 3 and retry.backoff_factor == 1.5
    assert 429 in retry.status_forcelist and 503 in retry.status_forcelist


def _run_cli(*args):
    return subprocess.run(
        [sys.executable, "-m", "polyphar.cli", *args], capture_output=True, text=True
    )


def test_cli_lists_subcommands():
    r = _run_cli("--help")
    assert r.returncode == 0
    for name in SUBCOMMANDS:
        assert name in r.stdout


def test_cli_unimplemented_exits_2(monkeypatch, capsys):
    from polyphar import cli

    monkeypatch.setitem(cli.STAGES, "uniprot", ("polyphar.stages._does_not_exist", "x"))
    assert cli.main(["uniprot", "--no-notify"]) == 2
    assert "not implemented" in capsys.readouterr().err


def test_cli_common_flags(monkeypatch):
    from polyphar.cli import build_parser

    monkeypatch.setenv("NPROCS", "7")
    a = build_parser().parse_args(["uniprot", "--dry-run", "--no-notify"])
    assert a.nprocs == 7 and a.dry_run and a.no_notify


def test_nprocs_default_and_empty_env(monkeypatch):
    from polyphar.cli import build_parser

    monkeypatch.delenv("NPROCS", raising=False)
    assert build_parser().parse_args(["uniprot"]).nprocs == 4
    monkeypatch.setenv("NPROCS", "")
    assert build_parser().parse_args(["uniprot"]).nprocs == 4


@pytest.mark.parametrize("bad", ["0", "-2", "abc"])
def test_nprocs_flag_rejects_invalid(bad):
    from polyphar.cli import build_parser

    with pytest.raises(SystemExit) as e:
        build_parser().parse_args(["uniprot", "--nprocs", bad])
    assert e.value.code == 2


def test_nprocs_bad_env_is_clear_error_and_help_works(monkeypatch, capsys):
    from polyphar.cli import build_parser

    monkeypatch.setenv("NPROCS", "0")
    p = build_parser()  # must not raise
    with pytest.raises(SystemExit) as e:
        p.parse_args(["uniprot"])
    assert e.value.code == 2
    assert "NPROCS" in capsys.readouterr().err
    with pytest.raises(SystemExit) as e:
        p.parse_args(["--help"])
    assert e.value.code == 0
