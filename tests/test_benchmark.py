import json
import os
from pathlib import Path

import pytest

from polyphar import benchmark, cli
from polyphar.common import RunLayout, WorkUnits, write_manifest
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
    """Real fix + similarity workloads. pockets is reported as
    host-calibrated and an existing entry is left exactly as it was."""
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
    assert [r["stage"] for r in rows] == ["structures", "fix", "pockets", "annotate",
                                          "select", "similarity", "ligands", "dock",
                                          "stats"]
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
    # dock docks the ligands `ligands` KEPT, not the source library's records: the
    # dockability filter drops ~6% of e-Drug3D_2083, so the raw count over-projects dock.
    dock_units = benchmark.projected_units("dock", state)
    assert benchmark.PRIORS["ligands_kept"] < benchmark.PRIORS["ligands"]
    assert dock_units.count == 3 * benchmark.PRIORS["ligands_kept"]
    assert benchmark.projected_units("ligands", state).count == benchmark.PRIORS["ligands"]

    # a finished ligands stage makes the kept count exact
    (layout.docking / "library").mkdir(parents=True)
    (layout.docking / "library" / "manifest.json").write_text(
        '{"stage": "ligands", "counts": {"records": 2083, "kept": 7}}')
    state, observed = benchmark.run_state(
        layout, ["ids", "uniprot", "select", "ligands"])
    assert state["ligands_kept"] == 7 and "ligands_kept" in observed
    assert benchmark.projected_units("dock", state).count == 3 * 7


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


def test_uniprot_bench_request_is_one_stage_batch():
    from polyphar.stages import uniprot

    assert benchmark.UNIPROT_BATCH == uniprot.BATCH


def test_uniprot_bench_sends_full_batches(tmp_path, monkeypatch):
    from polyphar.stages import uniprot

    seen = []

    def fake_fetch(ids, session, batch=100, on_batch=None):
        seen.append((len(ids), batch))
    monkeypatch.setattr(uniprot, "fetch_uniprot", fake_fetch)
    spc, extra = benchmark._bench_uniprot(tmp_path)
    assert seen == [(uniprot.BATCH, uniprot.BATCH)] * benchmark.UNIPROT_REQUESTS
    assert extra["accessions_per_request"] == uniprot.BATCH and spc > 0


def _snapshot(root):
    return {str(p.relative_to(root)): (p.stat().st_mtime_ns, p.stat().st_size)
            for p in sorted(root.rglob("*"))}


def _write_pockets(layout, names):
    import pandas as pd

    rows = []
    for name in names:
        d = layout.pockets / "P12345"
        d.mkdir(parents=True, exist_ok=True)
        (d / name).write_text("ATOM\n")
        rows.append({"UniProt ID": "P12345", "Gene Name": "G", "Protein": "1abc_A",
                     "Pocket": "P_1", "drugScore": 0.9, "Pocket Filename": name})
    pd.DataFrame(rows).to_csv(layout.drugscore_csv)
    # similarity now runs on the panel's eligible pockets, written by `select`. Every
    # pocket here is eligible, so the name guard is still what trips.
    layout.targets.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"UniProt_ID": "P12345", "Gene_Name": "G", "Protein": "1abc_A",
                   "Pocket": "P_1", "Pocket Filename": n,
                   "Pocket_Path": f"04_pockets/P12345/{n}", "AlphaFold": "False",
                   "drugScore": 0.9, "radius_of_gyration": 6.0} for n in names]
                 ).to_csv(layout.targets / "eligible_pockets.csv", index=False)


def test_estimate_is_read_only_even_when_the_name_guard_trips(tmp_path, monkeypatch):
    """similarity.units would record the guard's findings (failed.tsv, unlink manifest)
    unless dry-run; estimate must not touch the run folder."""
    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    layout, _ = _run_with_uniprot_done(tmp_path)
    layout.structures.mkdir(parents=True, exist_ok=True)
    # annotate and select are done too, so `similarity` is the NEXT stage and estimate
    # asks it for exact units — which is what makes it reach the name guard at all.
    for stage in ("structures", "fix", "pockets", "annotate", "select"):
        write_manifest(layout.run_dir / benchmark.STAGE_DIRS[stage], stage, {}, [], {}, {}, {})
    _write_pockets(layout, ["P12345_" + "x" * 40 + ".pdb"])   # > 32 characters
    before = _snapshot(layout.run_dir)
    rows, _ = benchmark.estimate_rows(layout, 4, _calib(tmp_path, similarity=0.002))
    sim_row = next(r for r in rows if r["stage"] == "similarity")
    assert sim_row["kind"] == "projected"
    assert "exact units unavailable" in sim_row["note"]   # the guard raised ...
    assert _snapshot(layout.run_dir) == before            # ... and recorded nothing
    assert benchmark._stage_args("similarity", layout, 4).dry_run is True
    assert not (layout.similarity / "failed.tsv").exists()


def test_run_state_observes_structures_fix_and_pockets(tmp_path):
    layout, _ = _run_with_uniprot_done(tmp_path)   # 4 accessions
    layout.structures.mkdir()
    (layout.structures / "A").mkdir()
    (layout.structures / "A" / "1abc_A.pdb").write_bytes(b"x" * 1000)
    (layout.structures / "A" / "2abc_A.pdb").write_bytes(b"x" * 3000)
    (layout.structures / "A" / "AF-A-F1.pdb").write_bytes(b"x" * 5000)
    layout.protein_info.write_text(
        "UniProt ID,File,AlphaFold\nA,02_structures/A/1abc_A.pdb,False\n"
        "A,02_structures/A/2abc_A.pdb,False\nA,02_structures/A/AF-A-F1.pdb,True\n")
    fixed = layout.pdbfixer / "A"
    fixed.mkdir(parents=True)
    (fixed / "1abc_A.pdb").write_bytes(b"x" * 2500)
    write_manifest(layout.pockets, "pockets", {}, [], {}, {"models": 4}, {})
    (layout.pockets / "DrugScore.csv").write_text(
        "idx,drugScore\n0,0.9\n1,0.7\n2,0.5\n3,0.4\n")   # 3 kept
    layout.targets.mkdir()
    layout.final_targets_csv.write_text("Protein\nA\nB\n")           # 2 targets
    layout.inputs.mkdir()
    (layout.inputs / "e-Drug3D_2083.sdf").write_text("a\n$$$$\nb\n$$$$\nc\n$$$$\n")   # 3
    state, observed = benchmark.run_state(
        layout, ["ids", "uniprot", "structures", "fix", "pockets", "select"])
    assert state["accessions"] == 4 and state["targets"] == 2 and state["ligands"] == 3
    assert state["af_models_per_accession"] == pytest.approx(1 / 4)
    assert state["pdb_chain_bytes"] == pytest.approx(2000)
    assert state["af_model_bytes"] == pytest.approx(5000)
    assert state["fixed_chain_bytes"] == pytest.approx(2500)
    assert state["pdb_chains_per_accession"] == pytest.approx(2 / 4)   # Protein_info PDB rows
    assert state["kept_pockets_per_model"] == pytest.approx(3 / 4)
    for key in ("accessions", "pdb_chains_per_accession", "targets", "ligands",
                "af_models_per_accession", "pdb_chain_bytes", "af_model_bytes",
                "fixed_chain_bytes", "kept_pockets_per_model"):
        assert key in observed and state[key] != benchmark.PRIORS[key]


def test_exact_row_comes_from_the_stage_units(tmp_path, monkeypatch):
    from polyphar.stages import structures

    monkeypatch.setenv("POLYPHAR_HOME", str(tmp_path / "home"))
    layout, _ = _run_with_uniprot_done(tmp_path)
    monkeypatch.setattr(structures, "units",
                        lambda layout, args: WorkUnits(7777, "sentinel"))
    rows, _ = benchmark.estimate_rows(layout, 4, _calib(tmp_path, structures=2.0))
    assert rows[0]["stage"] == "structures" and rows[0]["kind"] == "exact"
    assert rows[0]["units"] == "7777 (sentinel)"
    assert rows[0]["seconds"] == pytest.approx(7777 * 2.0 / 4)


def test_timed_request_counts_only_the_successful_attempt(monkeypatch):
    import time as _t

    calls = []

    def flaky():
        calls.append(1)
        _t.sleep(0.3 if len(calls) == 1 else 0.01)
        if len(calls) == 1:
            raise RuntimeError("500")

    seconds, failed = benchmark._timed_request(flaky, pause=0.2)
    assert failed == 1 and len(calls) == 2
    assert seconds < 0.2   # neither the failed attempt (0.3 s) nor the pause (0.2 s)
    with pytest.raises(RuntimeError):
        benchmark._timed_request(lambda: (_ for _ in ()).throw(RuntimeError("x")),
                                 attempts=2, pause=0)


def test_uniprot_bench_records_retries_and_excludes_failed_attempts(tmp_path, monkeypatch):
    import time as _t

    from polyphar.stages import uniprot

    state = {"n": 0}

    def flaky_fetch(ids, session, batch=100, on_batch=None):
        state["n"] += 1
        _t.sleep(0.01)
        if state["n"] == 2:   # the second request fails once
            _t.sleep(0.4)
            raise RuntimeError("500 error")

    monkeypatch.setattr(uniprot, "fetch_uniprot", flaky_fetch)
    monkeypatch.setattr(benchmark, "RETRY_PAUSE_S", 0)
    spc, extra = benchmark._bench_uniprot(tmp_path)
    assert extra["retries"] == 1 and state["n"] == benchmark.UNIPROT_REQUESTS + 1
    assert spc < 0.1   # 0.4 s of the failed attempt is not in the cost


def test_run_state_counts_selected_pockets_and_uses_the_tools_prior(tmp_path):
    layout, _ = _run_with_uniprot_done(tmp_path)
    state, _ = benchmark.run_state(layout, ["ids", "uniprot"], tool="dogsite3")
    assert state["kept_pockets_per_model"] == benchmark.KEPT_POCKETS_PRIOR["dogsite3"]
    state, _ = benchmark.run_state(layout, ["ids", "uniprot"], tool="dogsite2")
    assert state["kept_pockets_per_model"] == benchmark.PRIORS["kept_pockets_per_model"]
    write_manifest(layout.pockets, "pockets", {}, [], {}, {"models": 4}, {})
    (layout.pockets / "DrugScore.csv").write_text(
        "idx,volume,Tool,Selected\n0,400,dogsite3,True\n1,350,dogsite3,True\n"
        "2,100,dogsite3,False\n")
    state, observed = benchmark.run_state(layout, ["ids", "uniprot", "pockets"],
                                          tool="dogsite3")
    assert state["kept_pockets_per_model"] == 2 / 4 and "kept_pockets_per_model" in observed
