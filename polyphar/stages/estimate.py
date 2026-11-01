"""`polyphar estimate`: ETA table for the next stage (exact) and later ones (projected)."""

from __future__ import annotations

from polyphar.common import WorkUnits


def units(layout, args) -> WorkUnits:
    return WorkUnits(1, "estimate")


def run(layout, args, progress) -> dict:
    from polyphar.benchmark import estimate_rows, format_table
    from polyphar.eta import Calibration, polyphar_home

    rows, observed = estimate_rows(layout, args.nprocs, Calibration.load(polyphar_home()))
    print(format_table(rows, observed))
    return {"stages": len(rows), "failed": 0}
