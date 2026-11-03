"""`polyphar dock-worker --bundle DIR`: run a bundle's docking units where the engine is.

It needs the bundle and the engine, and nothing else: no run folder, no final_targets.csv,
no key arithmetic. The unit keys were decided by `dock-export` and are read from
units.tsv verbatim, so what the worker writes can be matched against the plan later.

Resume is the bundle's own: a unit whose marker holds its key and whose result file is
intact is SKIPPED -- it was drawn, and drawing it again would produce a different sample,
not the same one. A unit whose record is missing or truncated is discarded, logged and
re-drawn under the same replicate index.
"""

from __future__ import annotations

from polyphar.stages import dock

STAGE = "dock-worker"


def add_arguments(parser) -> None:
    parser.add_argument("--bundle", required=True,
                        help="bundle directory written by `polyphar dock-export`")
    parser.add_argument("--threads-per-job", type=dock._positive_int, default=dock.DEFAULT_TPJ,
                        help=f"engine threads per unit (default {dock.DEFAULT_TPJ})")
    parser.add_argument("--unit-timeout", type=dock._positive_int, default=None,
                        help="seconds before one unit's engine call is killed (default: the "
                             "bundle's)")


def _bundle(args):
    # argparse makes --bundle required, so this only fires for a programmatic caller.
    if not getattr(args, "bundle", None):
        raise ValueError("dock-worker needs --bundle <directory> "
                         "(written by `polyphar dock-export`)")
    return args.bundle


def units(layout, args):
    from polyphar.common import WorkUnits

    bundle = _bundle(args)
    key = dock.read_run_key(bundle)
    outstanding = [b for b in dock.read_bundle_units(bundle)
                   if dock.unit_state(dock.BundleLayout(bundle), b.unit) != dock.STATE_DONE]
    return WorkUnits(len(outstanding) * int(key.get("n_ligands") or 0), "ligand")


def run(layout, args, progress) -> dict:
    return dock.work_bundle(_bundle(args), args, progress)
