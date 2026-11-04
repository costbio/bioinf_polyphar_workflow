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
