import json
import os
from pathlib import Path

import pytest

from polyphar import benchmark, cli
from polyphar.common import RunLayout, write_manifest
from polyphar.eta import Calibration
from polyphar.stages.uniprot import parse_uniprot_tsv

FIX = Path(__file__).parent / "fixtures"


def _calib(tmp_path, **entries):
    c = Calibration(tmp_path / "home" / "calibration.json")
    for stage, spc in entries.items():
        c.update(stage, spc)
    return c


# --- benchmark -----------------------------------------------------------------------

@pytest.mark.pocketmatch
def test_benchmark_writes_calibration(tmp_path, monkeypatch):
    """Real fix + similarity workloads. Inside the container pockets is reported as
    host-calibrated and an existing entry is left exactly as it was."""
    monkeypatch.setenv("POLYPHAR_IN_CONTAINER", "1")
    assert benchmark.in_container()
    calib = _calib(tmp_path, pockets=3.2e-4)
    before = json.dumps(calib.data["pockets"], sort_keys=True)
    log = []
    res = benchmark.run_benchmark(["fix", "similarity", "pockets"], calib, log=log.append)
    assert res["fix"]["status"] == "calibrated" and res["similarity"]["status"] == "calibrated"
    saved = Calibration.load(tmp_path / "home")
    assert saved.cost("fix") > 0 and saved.cost("similarity") > 0
    assert saved.data["similarity"]["bench_pairs"] == 120
    assert res["pockets"]["status"] == "host-calibrated"
    assert "host-calibrated" in " ".join(log)
    assert json.dumps(saved.data["pockets"], sort_keys=True) == before


def test_pockets_is_never_run_here_and_never_created(tmp_path):
    calib = _calib(tmp_path)
    res = benchmark.run_benchmark(["pockets"], calib, log=lambda m: None)
    assert res["pockets"]["status"] == "host-calibrated"
    assert Calibration.load(tmp_path / "home").cost("pockets") is None


def test_benchmark_replaces_entries_skips_and_records_dock_mix(tmp_path, monkeypatch):
    calib = _calib(tmp_path, fix=9.0)
    calib.update("fix", 11.0)   # a mean of two observations: the benchmark replaces it
    monkeypatch.setitem(benchmark.WORKLOADS, "fix", lambda tmp: (2e-5, {"basis": "bytes"}))

    def no_pm(tmp):
        raise benchmark.Skip("PocketMatch binaries not configured")

    monkeypatch.setitem(benchmark.WORKLOADS, "similarity", no_pm)

    def boom(tmp):
        raise RuntimeError("smina exploded")

    monkeypatch.setitem(benchmark.WORKLOADS, "dock", boom)
    res = benchmark.run_benchmark(["fix", "similarity", "dock"], calib, log=lambda m: None)
    saved = Calibration.load(tmp_path / "home")
    assert saved.cost("fix") == pytest.approx(2e-5)
    assert saved.data["fix"]["n_observations"] == 1
    assert res["similarity"] == {"status": "skipped",
                                 "reason": "PocketMatch binaries not configured"}
    assert saved.cost("similarity") is None
    assert res["dock"]["status"] == "failed" and "smina exploded" in res["dock"]["reason"]


def test_dock_calibration_carries_the_library_mix_correction(tmp_path, monkeypatch):
    extra = {"basis": "ligands", "library_mix_factor": benchmark.DOCK_MIX_FACTOR,
             "library_mix_source": benchmark.DOCK_MIX_SOURCE}
    monkeypatch.setitem(benchmark.WORKLOADS, "dock", lambda tmp: (30.0 * 71.5 / 45, extra))
    calib = _calib(tmp_path)
    benchmark.run_benchmark(["dock"], calib, log=lambda m: None)
    e = Calibration.load(tmp_path / "home").data["dock"]
    assert e["seconds_per_unit"] == pytest.approx(30.0 * 71.5 / 45)
    assert e["library_mix_factor"] == pytest.approx(1.5889, abs=1e-3)
    assert "71.5" in e["library_mix_source"]


def test_bench_data_is_packaged():
    assert benchmark.FIX_CHAINS[1] == "2mjw_B.pdb"
    for name in (*benchmark.FIX_CHAINS, benchmark.DOCK_RECEPTOR, benchmark.DOCK_POCKET,
                 benchmark.DOCK_LIGANDS, "accessions.txt", "NOTICE"):
        assert (benchmark.BENCH_DIR / name).stat().st_size > 0
    assert len(list((benchmark.BENCH_DIR / "sample_pockets").glob("*.pdb"))) == 16
    sizes = sorted((benchmark.BENCH_DIR / n).stat().st_size for n in benchmark.FIX_CHAINS)
    assert sizes[-1] > 3 * sizes[0]
    worker = (Path(__file__).parents[1] / "bin" / "polyphar-dogsite-worker").read_text()
    assert all(f'"{n}"' in worker for n in benchmark.FIX_CHAINS)


def test_cli_benchmark_pockets_only(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    assert cli.main(["benchmark", "--run-dir", str(tmp_path), "--stages", "pockets"]) == 0
    assert "host-calibrated" in capsys.readouterr().out


# --- estimate ------------------------------------------------------------------------

def _run_with_uniprot_done(tmp_path):
    layout = RunLayout(tmp_path / "run")
    df = parse_uniprot_tsv((FIX / "uniprot_sample.tsv").read_text())
    layout.uniprot_csv.parent.mkdir(parents=True)
    df.to_csv(layout.uniprot_csv, index=False)
    write_manifest(layout.uniprot_csv.parent, "uniprot", {}, [], {}, {"kept": len(df)}, {})
    return layout, len(df)


def test_estimate_table_marks_projection(tmp_path, monkeypatch):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    layout, n = _run_with_uniprot_done(tmp_path)
    calib = _calib(tmp_path, structures=1.5, fix=2e-5, similarity=0.002, dock=70.0)
    rows, observed = benchmark.estimate_rows(layout, 4, calib)
    by = {r["stage"]: r for r in rows}
    assert [r["stage"] for r in rows] == ["structures", "fix", "pockets", "similarity",
                                          "select", "dock"]
    assert by["structures"]["kind"] == "exact" and by["structures"]["seconds"] > 0
    assert by["dock"]["kind"] == "projected" and by["fix"]["kind"] == "projected"
    assert by["select"]["kind"] == "trivial"
    assert "accessions" in observed   # uniprot.csv rows replace the prior
    assert by["pockets"]["seconds"] is None   # uncalibrated is marked, not guessed
    text = benchmark.format_table(rows, observed)
    line = {l.split()[0]: l for l in text.splitlines()[1:] if l.split()}
    assert " exact " in line["structures"] and " projected " in line["dock"]
    assert "uncalibrated" in line["pockets"]
    assert "not counted: pockets" in line["total:"]


def test_estimate_uses_observed_values_over_priors(tmp_path):
    layout, n = _run_with_uniprot_done(tmp_path)
    state, observed = benchmark.run_state(layout, ["ids", "uniprot"])
    assert state["accessions"] == n and state["targets"] == benchmark.PRIORS["targets"]
    # a finished select replaces the target-count prior with final_targets.csv's rows
    layout.targets.mkdir()
    layout.final_targets_csv.write_text("Protein\nA\nB\nC\n")
    state, observed = benchmark.run_state(layout, ["ids", "uniprot", "select"])
    assert state["targets"] == 3 and "targets" in observed
    dock_units = benchmark.projected_units("dock", state)
    assert dock_units.count == 3 * benchmark.PRIORS["ligands"]


def test_cli_estimate_is_fast_and_offline(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    layout, _ = _run_with_uniprot_done(tmp_path)
    assert cli.main(["estimate", "--run-dir", str(layout.run_dir)]) == 0
    out = capsys.readouterr().out
    assert "structures" in out and "projected" in out and "total:" in out


def test_estimate_on_empty_run_and_finished_run(tmp_path, monkeypatch):
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    layout = RunLayout(tmp_path / "empty")
    layout.run_dir.mkdir()
    rows, observed = benchmark.estimate_rows(layout, 4, _calib(tmp_path, uniprot=1.0))
    assert rows[0]["stage"] == "ids" and rows[0]["kind"] == "trivial"
    assert benchmark.format_table([], []) .startswith("all stages are complete")
