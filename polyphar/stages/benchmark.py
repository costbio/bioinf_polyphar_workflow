"""`polyphar benchmark`: per-machine calibration (see polyphar.benchmark)."""

from __future__ import annotations

from polyphar.common import WorkUnits


def add_arguments(parser) -> None:
    from polyphar.benchmark import BENCH_STAGES

    parser.add_argument("--stages", nargs="+", default=list(BENCH_STAGES), metavar="STAGE",
                        choices=list(BENCH_STAGES),
                        help="stages to calibrate (default: all; pockets is host-calibrated)")


def units(layout, args) -> WorkUnits:
    return WorkUnits(len(args.stages), "stages")


def run(layout, args, progress) -> dict:
    from polyphar.benchmark import run_benchmark
    from polyphar.eta import Calibration, polyphar_home

    calib = Calibration.load(polyphar_home())
    results = run_benchmark(args.stages, calib, log=progress if callable(progress) else print)
    print(f"calibration: {calib.path}")
    return {"calibrated": sum(r["status"] == "calibrated" for r in results.values()),
            "skipped": sum(r["status"] in ("skipped", "host-calibrated")
                           for r in results.values()),
            "failed": sum(r["status"] == "failed" for r in results.values())}
