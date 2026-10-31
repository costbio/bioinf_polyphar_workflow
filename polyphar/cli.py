"""polyphar command line.

Each subcommand is bound to a stage module in STAGES. A stage module exposes
    units(layout, args) -> WorkUnits
    run(layout, args, progress) -> dict
To wire a new stage, add or change one line in STAGES; subcommands whose module
does not exist yet print "not implemented" and exit 2.
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path

from polyphar import __version__
from polyphar.common import RunLayout

# subcommand -> (stage module, help text)
STAGES = {
    "ids": ("polyphar.stages.ids", "build 00_ids/ids.txt from an xlsx or txt source"),
    "uniprot": ("polyphar.stages.uniprot", "fetch UniProt records"),
    "structures": ("polyphar.stages.structures", "download PDB/AlphaFold structures"),
    "fix": ("polyphar.stages.fix", "repair structures with pdbfixer"),
    "pockets": ("polyphar.stages.pockets", "detect pockets (DoGSiteScorer, host only)"),
    "similarity": ("polyphar.stages.similarity", "PocketMatch pocket similarity"),
    "select": ("polyphar.stages.select", "select polypharmacology targets"),
    "dock": ("polyphar.stages.dock", "dock the ligand library with smina"),
    "run": ("polyphar.stages.run", "run all stages in order"),
    "estimate": ("polyphar.stages.estimate", "print ETA for the next and later stages"),
    "benchmark": ("polyphar.stages.benchmark", "per-machine calibration"),
}


def default_run_dir() -> str:
    return "/work" if Path("/work").is_dir() else os.getcwd()


def default_nprocs() -> int:
    return int(os.environ.get("NPROCS", 4))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="polyphar", description="Polypharmacology workflow.")
    parser.add_argument("--version", action="version", version=f"polyphar {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>", required=True)
    for name, (_module, help_text) in STAGES.items():
        p = sub.add_parser(name, help=help_text, description=help_text)
        p.add_argument("--run-dir", default=default_run_dir(),
                       help="run folder (default: /work if it exists, else cwd)")
        p.add_argument("--nprocs", type=int, default=default_nprocs(),
                       help="worker processes (default: $NPROCS or 4)")
        p.add_argument("--dry-run", action="store_true",
                       help="print units and ETA, then exit")
        p.add_argument("--no-notify", action="store_true",
                       help="do not send notifications")
    return parser


def _load_stage(name: str):
    module_name = STAGES[name][0]
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        # Only "this stage does not exist yet"; a missing dependency inside an
        # existing stage must still raise.
        if exc.name and (module_name == exc.name or module_name.startswith(exc.name + ".")):
            return None
        raise


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    stage = _load_stage(args.command)
    if stage is None:
        print(f"polyphar {args.command}: not implemented", file=sys.stderr)
        return 2
    layout = RunLayout(args.run_dir)
    if args.dry_run:
        wu = stage.units(layout, args)
        print(f"{args.command}: {wu.count} units ({wu.cost_basis})")
        return 0
    stage.run(layout, args, print)
    return 0


if __name__ == "__main__":
    sys.exit(main())
