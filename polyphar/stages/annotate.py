"""Per-pocket geometric annotation, derived entirely from existing stage outputs.

Why this is a separate stage and not part of `pockets`: `pockets` cannot be re-run
without DoGSiteScorer, whose licence is host-bound and time-limited. Everything here is
computed from files already on disk (`04_pockets/<acc>/*_res_P_*.pdb` and, for AlphaFold
models, `02_structures/<acc>/<model>.pdb`), so it can be re-run freely, and it never
touches DoGSiteScorer's own output.

Outputs in 04_pockets/annotation/: pocket_annotation.csv (one row per selected pocket,
keyed by `Pocket Filename`), failed.tsv, manifest.json (written last = done).

Note on the output path: the design named `04_pockets/pocket_annotation.csv`. It lives in
a subdirectory instead so this stage can own a manifest.json and failed.tsv without
clashing with the ones `pockets` writes beside DrugScore.csv.

Columns
-------
UniProt ID, Gene Name, Protein, Pocket, Pocket Filename, AlphaFold
    identity, copied from DrugScore.csv.
n_atoms, n_residues
    lining atoms and distinct lining residues in the `_res_P_` file.
radius_of_gyration, segments, span
    DESCRIPTIVE ONLY — none of these filters the panel. Both were tried as compactness
    criteria and both failed: `segments <= 3` is provenance-dependent (AlphaFold pockets
    pass at 39.4% against 17.5% for experimental ones, because pdbfixer renumbers
    experimental chains 1..N with missing loops modelled in while AlphaFold pockets
    carry UniProt positions), and `Rg <= 8.5` kept only 8% of pockets holding a
    co-crystal drug-like ligand against 69% of the rest — real binding sites have LARGE
    lining-residue envelopes. The panel filters on volume, surf/vol and drugScore
    instead; see polyphar/panel.py. These columns are kept because they are cheap and
    were the evidence for that conclusion.
box_center_{x,y,z}, box_size_{x,y,z}, box_max_dim, box_volume
    the unpadded bounding box of the lining atoms, i.e. exactly what
    `dock.pocket_box` computes, so the panel can be filtered on the box that will
    actually be searched.
grid_points_0375, grid_legal_0375
    Uni-Dock's grid at 0.375 A spacing, and whether it is within the engine's limits
    (<=128 points per axis and <=531,441 in total). Uni-Dock exits on a violation;
    versions before 1.2.0 silently corrupted results instead.
parent_volume
    for a subpocket `P_n_m`, the volume of its parent `P_n`, so a subpocket carved out
    of an oversized cleft can be rejected. Empty for a top-level pocket.
residues
    the lining residue numbers, range-encoded (`10-15,30,45-50`). Carried so the TM
    pocket-overlap rule can be evaluated in `select` without re-reading 11,599 pocket
    files. For an AlphaFold model these are UniProt positions; for an experimental
    chain they are pdbfixer numbering, which is why they are usable for the AF-only
    overlap check but not for a cross-class comparison.
pocket_plddt_mean, pocket_plddt_min
    mean and minimum pLDDT over lining residues, for AlphaFold models only, read back
    from the model in 02_structures because DoGSiteScorer zeroes the B-factor column in
    its pocket files. Recorded as a sanity floor; it does **not** detect the oversized
    interdomain cavities (top-degree AF pockets average 84.8 against 83.3 overall),
    which is a job for PAE, not pLDDT.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

from polyphar import pool
from polyphar.common import WorkUnits, append_failure, write_manifest

STAGE = "annotate"
OUTPUT = "pocket_annotation.csv"
OUTPUTS = (OUTPUT,)
SPACING = 0.375          # Uni-Dock default grid spacing
MAX_AXIS = 128           # Uni-Dock: grid points per axis
MAX_POINTS = 531441      # Uni-Dock: total grid points (81^3)
SEGMENT_GAP = 5          # a gap of more than this many residues starts a new segment

COLUMNS = [
    "UniProt ID", "Gene Name", "Protein", "Pocket", "Pocket Filename", "AlphaFold",
    "n_atoms", "n_residues", "radius_of_gyration", "segments", "span",
    "box_center_x", "box_center_y", "box_center_z",
    "box_size_x", "box_size_y", "box_size_z", "box_max_dim", "box_volume",
    "grid_points_0375", "grid_legal_0375", "parent_volume", "residues",
    "pocket_plddt_mean", "pocket_plddt_min",
]


def add_arguments(parser) -> None:
    parser.add_argument("--spacing", type=float, default=SPACING, metavar="A",
                        help=f"grid spacing for the Uni-Dock legality check "
                             f"(default: {SPACING})")


# --- geometry ------------------------------------------------------------------------

def read_pocket(path):
    """(atom coordinates, sorted distinct residue numbers) from a pocket PDB."""
    xyz, res = [], set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
            try:
                res.add(int(line[22:26]))
            except ValueError:
                pass
    return xyz, sorted(res)


def radius_of_gyration(xyz) -> float:
    """RMS distance of the points from their centroid."""
    n = len(xyz)
    c = [sum(p[k] for p in xyz) / n for k in range(3)]
    return math.sqrt(sum(sum((p[k] - c[k]) ** 2 for k in range(3)) for p in xyz) / n)


def bounding_box(xyz):
    """(centre, size) of the unpadded bounding box — the same box `dock.pocket_box` uses."""
    lo = [min(p[k] for p in xyz) for k in range(3)]
    hi = [max(p[k] for p in xyz) for k in range(3)]
    return ([(hi[k] + lo[k]) / 2 for k in range(3)], [hi[k] - lo[k] for k in range(3)])


def count_segments(res, gap: int = SEGMENT_GAP) -> int:
    """Contiguous runs of lining residues, a gap of more than `gap` starting a new one."""
    if not res:
        return 0
    segs, prev = 1, res[0]
    for r in res[1:]:
        if r - prev > gap:
            segs += 1
        prev = r
    return segs


def grid(size, spacing: float = SPACING):
    """(points per axis, total, legal) for Uni-Dock at `spacing`."""
    n = [math.floor(s / spacing) + 1 for s in size]
    total = n[0] * n[1] * n[2]
    return n, total, (max(n) <= MAX_AXIS and total <= MAX_POINTS)


def residue_plddt(model_pdb):
    """{residue number: pLDDT} from an AlphaFold model's B-factor column."""
    out = {}
    with open(model_pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                try:
                    out[int(line[22:26])] = float(line[60:66])
                except ValueError:
                    pass
    return out


def encode_residues(res) -> str:
    """[10,11,12,30,45,46] -> '10-12,30,45-46'. Compact, and round-trips exactly."""
    if not res:
        return ""
    out, start, prev = [], res[0], res[0]
    for r in res[1:]:
        if r == prev + 1:
            prev = r
            continue
        out.append(f"{start}-{prev}" if start != prev else f"{start}")
        start = prev = r
    out.append(f"{start}-{prev}" if start != prev else f"{start}")
    return ",".join(out)


def decode_residues(text: str) -> list[int]:
    """Inverse of encode_residues."""
    out = []
    for part in str(text).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part[1:]:
            a, b = part.rsplit("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def receptor_stem(pocket_filename: str) -> str:
    """`8vut_A_res_P_1.pdb` -> `8vut_A`; `AF-P1-F1_res_P_2.pdb` -> `AF-P1-F1`."""
    return pocket_filename.split("_res_P")[0]


# --- one pocket ----------------------------------------------------------------------

def annotate_one(row: dict, pockets_dir: str, structures_dir: str,
                 parent_volume, spacing: float = SPACING) -> dict:
    """Annotate one pocket. Raises on a missing or unreadable pocket file."""
    acc = row["UniProt ID"]
    fname = row["Pocket Filename"]
    path = Path(pockets_dir) / acc / fname
    xyz, res = read_pocket(path)
    if len(xyz) < 4:
        raise ValueError(f"{fname}: {len(xyz)} atoms, too few to annotate")

    centre, size = bounding_box(xyz)
    n, total, legal = grid(size, spacing)
    is_af = str(row.get("Protein", "")).startswith("AF-")

    pl_mean = pl_min = ""
    if is_af:
        model = Path(structures_dir) / acc / f"{receptor_stem(fname)}.pdb"
        if model.is_file():
            per_res = residue_plddt(model)
            vals = [per_res[r] for r in res if r in per_res]
            if vals:
                pl_mean, pl_min = round(sum(vals) / len(vals), 2), round(min(vals), 2)

    return {
        "UniProt ID": acc,
        "Gene Name": row.get("Gene Name", ""),
        "Protein": row.get("Protein", ""),
        "Pocket": row.get("Pocket", ""),
        "Pocket Filename": fname,
        "AlphaFold": is_af,
        "n_atoms": len(xyz),
        "n_residues": len(res),
        "radius_of_gyration": round(radius_of_gyration(xyz), 4),
        "segments": count_segments(res),
        "span": (res[-1] - res[0] + 1) if res else 0,
        "box_center_x": round(centre[0], 4),
        "box_center_y": round(centre[1], 4),
        "box_center_z": round(centre[2], 4),
        "box_size_x": round(size[0], 3),
        "box_size_y": round(size[1], 3),
        "box_size_z": round(size[2], 3),
        "box_max_dim": round(max(size), 3),
        "box_volume": round(size[0] * size[1] * size[2], 1),
        "grid_points_0375": total,
        "grid_legal_0375": legal,
        "parent_volume": "" if parent_volume is None else round(parent_volume, 2),
        "residues": encode_residues(res),
        "pocket_plddt_mean": pl_mean,
        "pocket_plddt_min": pl_min,
    }


def _worker(row, pockets_dir, structures_dir, parent_volume, spacing):
    return annotate_one(row, pockets_dir, structures_dir, parent_volume, spacing)


# --- stage ---------------------------------------------------------------------------

def _selected(layout) -> list[dict]:
    """Selected rows of DrugScore.csv, in file order."""
    path = layout.drugscore_csv
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found: run the pockets stage first")
    with open(path, newline="", encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh)
                if str(r.get("Selected", "")).strip().lower() in ("true", "1", "yes")]


def _parent_volumes(layout) -> dict:
    """{(accession, protein, pocket id): volume} over every pocket, for the subpocket
    rule. The key MUST include the structure: `P_1` exists for every structure of an
    accession, so keying on (accession, pocket) alone silently returns some other
    structure's pocket volume."""
    out = {}
    with open(layout.drugscore_csv, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                out[(r["UniProt ID"], r["Protein"], r["Pocket"])] = float(r["volume"])
            except (KeyError, TypeError, ValueError):
                continue
    return out


def parent_volume_of(acc: str, protein: str, pocket: str, volumes: dict):
    """Volume of `P_n` for a subpocket `P_n_m` of the SAME structure; None for a
    top-level pocket."""
    parts = pocket.split("_")
    if len(parts) <= 2:
        return None
    return volumes.get((acc, protein, "_".join(parts[:-1])))


def _size(layout, row) -> float:
    """Pocket file size, used as the ETA weight (this stage is I/O bound)."""
    p = layout.pockets / row["UniProt ID"] / row["Pocket Filename"]
    try:
        return float(p.stat().st_size)
    except OSError:
        return 0.0


def units(layout, args) -> WorkUnits:
    rows = _selected(layout)
    return WorkUnits(len(rows), "pockets", [_size(layout, r) for r in rows])


def _versions() -> dict:
    return {"spacing": SPACING, "segment_gap": SEGMENT_GAP,
            "max_axis": MAX_AXIS, "max_points": MAX_POINTS}


def run(layout, args, progress) -> dict:
    spacing = float(getattr(args, "spacing", SPACING))
    d = layout.pockets / "annotation"
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # the outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")

    rows = _selected(layout)
    volumes = _parent_volumes(layout)
    results: dict[str, dict] = {}
    failures: list[tuple[str, str]] = []
    weight = {f'{r["UniProt ID"]}/{r["Pocket Filename"]}': _size(layout, r) for r in rows}

    def on_result(item, result):
        results[item] = result
        progress.advance(1, size=weight.get(item, 0.0))

    def on_failure(item, why):
        failures.append((item, why))
        progress.advance(1, size=weight.get(item, 0.0))

    jobs = []
    for r in rows:
        item = f'{r["UniProt ID"]}/{r["Pocket Filename"]}'
        pv = parent_volume_of(r["UniProt ID"], r.get("Protein", ""),
                              r.get("Pocket", ""), volumes)
        jobs.append((item, _worker,
                     (dict(r), str(layout.pockets), str(layout.structures), pv, spacing)))
    try:
        pool.run_pool(jobs, max(1, getattr(args, "nprocs", 1)),
                      on_result=on_result, on_failure=on_failure)
    finally:
        for item, reason in sorted(failures):
            append_failure(d, item, STAGE, reason)

    # written in DrugScore.csv order so the output is reproducible whatever the pool did
    out = d / OUTPUT
    tmp = out.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            got = results.get(f'{r["UniProt ID"]}/{r["Pocket Filename"]}')
            if got:
                w.writerow(got)
    os.replace(tmp, out)

    af = sum(1 for v in results.values() if v["AlphaFold"])
    counts = {"selected": len(rows), "annotated": len(results),
              "alphafold": af, "experimental": len(results) - af,
              "grid_illegal": sum(1 for v in results.values()
                                  if not v["grid_legal_0375"]),
              "failed": len(failures)}
    # The manifest is written last: its presence marks the stage as finished.
    write_manifest(d, STAGE,
                   args={"nprocs": getattr(args, "nprocs", 1), "spacing": spacing},
                   inputs=[layout.drugscore_csv], versions=_versions(),
                   counts=counts, timing={},
                   extra={"outputs": {OUTPUT: (d / OUTPUT).stat().st_size}})
    return counts
