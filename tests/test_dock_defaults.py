"""The docking engine default. No openbabel/rdkit needed, so it runs everywhere
(tests/test_stage_dock.py is importorskip'd on openbabel and is silent without it).
"""

from polyphar.stages import dock


def test_the_default_engine_is_unidock_and_control_inherits_it():
    """Uni-Dock on a GPU is the production route; smina is opt-in with --engine smina.

    `control` must inherit the same default. Its library-percentile readout scores the
    co-crystal ligand against the PRODUCTION library at the same target, so a control run
    on a different engine than `dock` compares two score scales -- Uni-Dock reads
    +0.326 kcal/mol less negative than smina (measured, n=580, rho=0.951) -- and biases
    the percentile without saying so.
    """
    from polyphar.cli import build_parser

    assert dock.DEFAULT_ENGINE == dock.ENGINE_UNIDOCK

    for stage in ("dock", "dock-export", "control"):
        args = build_parser().parse_args([stage, "--run-dir", "/tmp/x"]
                                         + (["--bundle", "/tmp/b"] if stage == "dock-export" else []))
        assert dock.engine_name(args) == dock.ENGINE_UNIDOCK, stage
        taken = build_parser().parse_args(
            [stage, "--run-dir", "/tmp/x", "--engine", "smina"]
            + (["--bundle", "/tmp/b"] if stage == "dock-export" else []))
        assert dock.engine_name(taken) == dock.ENGINE_SMINA, stage


class _U:
    def __init__(self, label, replicate):
        self.label, self.replicate = label, replicate

    def __repr__(self):
        return f"{self.label}/r{self.replicate}"


def _units(*pairs):
    return [_U(l, r) for l, r in pairs]


def test_a_bundle_never_splits_one_targets_replicates_across_sessions():
    """Target-major export, truncated to whole targets.

    The replicates measure the ENGINE's non-determinism. A session gets whatever GPU the
    provider has free, so a target whose replicates land in two sessions can be drawn on
    two GPU models, and any per-model score difference then sits inside the replicate
    spread looking like engine jitter.
    """
    todo = _units(("A", 0), ("A", 1), ("A", 2), ("B", 0), ("B", 1), ("B", 2))

    # asks for 4: taking 4 would split B, so it stops after A
    assert [repr(u) for u in dock._whole_targets(todo, 4)] == ["A/r0", "A/r1", "A/r2"]
    # an exact target boundary is kept whole
    assert len(dock._whole_targets(todo, 3)) == 3
    assert len(dock._whole_targets(todo, 6)) == 6
    # more than there is: everything
    assert len(dock._whole_targets(todo, 99)) == 6


def test_a_max_units_smaller_than_one_target_still_exports_one_whole_target():
    """Returning nothing would make a too-small --max-units silently export an empty
    bundle, which looks like "no work outstanding" rather than a bad argument."""
    todo = _units(("A", 0), ("A", 1), ("A", 2), ("B", 0))
    got = dock._whole_targets(todo, 1)
    assert [repr(u) for u in got] == ["A/r0", "A/r1", "A/r2"]


def test_outstanding_replicates_of_one_target_stay_together_even_when_partly_done():
    """r0 already recorded: the remaining r1/r2 are still one target's worth and go
    together."""
    todo = _units(("A", 1), ("A", 2), ("B", 1), ("B", 2))
    assert [repr(u) for u in dock._whole_targets(todo, 3)] == ["A/r1", "A/r2"]


def test_control_mode_swaps_the_layout_and_the_library_not_a_separate_stage():
    """`--control` is a flag on dock-export/dock-merge, not a control-export triple.

    `control` owns no docking code: it calls dock.run() with a ControlLayout (a RunLayout
    whose `docking` is 08_control/docking) and the co-crystal SDF as `library`. Export and
    merge need exactly that substitution. dock-worker needs nothing: a bundle is
    self-contained.
    """
    import argparse
    from pathlib import Path

    from polyphar.common import RunLayout
    from polyphar.stages import control as control_stage

    run = Path("/tmp/polyphar-control-view-test")
    layout = RunLayout(run)

    plain, lib = dock.control_view(layout, argparse.Namespace(control=False))
    assert plain is layout and lib is None

    # without the library, the error names the command that writes it
    try:
        dock.control_view(layout, argparse.Namespace(control=True))
    except FileNotFoundError as exc:
        assert "polyphar control" in str(exc)
    else:                                                    # pragma: no cover
        raise AssertionError("a missing control library must be refused")

    # with it: ControlLayout, and docking redirected away from 07_docking
    libp = control_stage.library_path(layout)
    libp.parent.mkdir(parents=True, exist_ok=True)
    libp.write_text("")
    try:
        got, got_lib = dock.control_view(layout, argparse.Namespace(control=True))
        assert isinstance(got, control_stage.ControlLayout)
        assert got_lib == libp
        assert got.docking != layout.docking
        assert "08_control" in str(got.docking)
    finally:
        libp.unlink(missing_ok=True)


def test_control_flag_is_on_export_and_merge_but_not_on_the_worker():
    """The worker runs a bundle and needs no run folder, so it has no --control."""
    from polyphar.cli import build_parser

    for stage, extra in (("dock-export", []), ("dock-merge", ["--bundle", "/tmp/b"])):
        a = build_parser().parse_args([stage, "--run-dir", "/tmp/x", "--control"] + extra)
        assert a.control is True, stage

    worker = build_parser().parse_args(["dock-worker", "--bundle", "/tmp/b"])
    assert not hasattr(worker, "control") or worker.control is False


def test_dock_refuses_a_raw_library_when_ligands_has_not_run_in_a_pipeline_folder(tmp_path):
    """A deleted `ligands` directory looks exactly like a stage that never ran.

    On 2026-11-04 the prepared library was removed by hand as a presumed leftover, and
    dock-export silently fell back to the raw source: 2083 unprepared records per unit
    instead of 1960 prepared ones -- no valence repair, no MMFF94 minimisation, and the
    123 records `ligands` excludes back in the library. prepared_state() documents "no
    manifest: not a warning", which is right for a fixture and wrong for a run folder
    whose panel already exists.
    """
    import argparse

    from polyphar.common import RunLayout
    from polyphar.stages import dock as dock_mod

    run = tmp_path / "run"
    (run / "inputs").mkdir(parents=True)
    raw = run / "inputs" / "e-Drug3D_2083.sdf"
    raw.write_text("$$$$\n")
    layout = RunLayout(run)
    args = argparse.Namespace(library=dock_mod.DEFAULT_LIBRARY)

    # no panel yet: nothing to refuse, the raw library is all there is
    assert dock_mod._library_path(layout, args, copy=False) == raw

    # a panel exists -> this is a pipeline run, and ligands belongs before dock
    layout.final_targets_csv.parent.mkdir(parents=True, exist_ok=True)
    layout.final_targets_csv.write_text("Node\n")
    try:
        dock_mod._library_path(layout, args, copy=False)
    except FileNotFoundError as exc:
        assert "polyphar ligands" in str(exc)
        assert "--library" in str(exc)          # the escape hatch is named
    else:                                        # pragma: no cover
        raise AssertionError("dock must refuse to dock an unprepared library here")

    # naming a library explicitly is a deliberate choice and still allowed
    other = run / "inputs" / "mine.sdf"
    other.write_text("$$$$\n")
    assert dock_mod._library_path(layout, argparse.Namespace(library="inputs/mine.sdf"),
                                  copy=False) == other
