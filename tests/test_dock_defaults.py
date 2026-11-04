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
