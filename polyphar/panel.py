"""Panel selection: the eligible-pocket filters and one-row-per-protein rule.

This replaces hub-and-neighbours selection. The hub was the highest-degree node of the
pocket-similarity graph, which made the target count a property of how many structures
happened to be processed rather than of biology, and anchored everything on one pocket
chosen by network topology. Pocket similarity is no longer used to select targets at
all; it is applied post hoc to the docking results.

The filters, in the order the funnel reports them:

  drugScore >= 0.5          DoGSiteScorer 2 selection, as in 2023
  300 <= volume <= 2000     the UPPER bound is new. drugScore rewards large pockets and
                            has no ceiling, so oversized interdomain clefts dominated the
                            old network: the top-degree AlphaFold pockets averaged 90
                            lining residues across 7 sequence segments at 2,809 A^3,
                            against 16 / 3 / 363 for a typical one.
  subpocket parent <= 2000  DoGSiteScorer emits P_n and P_n_m. Without this a 1,500 A^3
                            subpocket of a 25,840 A^3 cleft passes the volume ceiling and
                            becomes its protein's representative, defeating it entirely.
  AF pocket pLDDT >= 70     a sanity floor only. pLDDT does NOT detect the oversized
                            cavities (top-degree AF pockets average 84.8 against 83.3
                            overall) because it is local per-residue confidence, not
                            domain packing; PAE would be the right metric.
  TM policy                 see tm_reject(). The artifact is multi-pass proteins, whose
                            single extracted chains yield lipid-facing grooves and pore
                            lumens.
  radius of gyration <= 8.5 the compactness criterion, measured in SPACE. A
                            sequence-segment count was tried first and is provenance-
                            biased: pdbfixer renumbers experimental chains 1..N with
                            missing loops modelled in, while AlphaFold pockets carry
                            UniProt positions, so `segments <= 3` passed AF pockets at
                            39.4% against 17.5% for experimental ones. Rg closes that
                            21.9-point gap to 1.5, and bounds the docking box as a side
                            effect, so no separate box cap is needed.
  Uni-Dock grid legal       <=128 grid points per axis and <=531,441 in total at 0.375 A.
                            The engine exits on a violation; versions before 1.2.0
                            silently corrupted results instead. Removes nothing once Rg
                            is applied, and is kept as a guard.

One row per protein: the eligible pocket with the highest drugScore, ties broken by the
smallest pocket name so the result never depends on row order.
"""

from __future__ import annotations

import math

DRUGSCORE_MIN = 0.5
VOLUME_MIN = 300.0
VOLUME_MAX = 2000.0
PARENT_VOLUME_MAX = 2000.0
PLDDT_MIN = 70.0
RG_MAX = 8.5
TM_OVERLAP_MAX = 0.10
SPACING = 0.375
MAX_AXIS = 128
MAX_POINTS = 531441
TM_POLICIES = ("A", "B", "C")

# The funnel steps, in order. Each pocket is attributed to the FIRST step it fails, so
# the counts form a chain and sum correctly.
STEPS = [
    "drugscore", "volume", "subpocket_parent", "plddt", "tm", "rg", "grid",
]
STEP_REASON = {
    "drugscore": "drugScore below the cutoff",
    "volume": "pocket volume outside the allowed range",
    "subpocket_parent": "subpocket of a parent above the volume ceiling",
    "plddt": "AlphaFold pocket pLDDT below the floor",
    "tm": "rejected by the transmembrane policy",
    "rg": "radius of gyration above the compactness cutoff",
    "grid": "docking box exceeds Uni-Dock's grid limits",
}


class Options:
    """Filter thresholds. Defaults are the decided panel definition."""

    def __init__(self, drugscore_min=DRUGSCORE_MIN, volume_min=VOLUME_MIN,
                 volume_max=VOLUME_MAX, parent_volume_max=PARENT_VOLUME_MAX,
                 plddt_min=PLDDT_MIN, rg_max=RG_MAX, tm_policy="C",
                 tm_overlap_max=TM_OVERLAP_MAX, spacing=SPACING):
        if tm_policy not in TM_POLICIES:
            raise ValueError(f"tm_policy must be one of {TM_POLICIES}, got {tm_policy!r}")
        self.drugscore_min = float(drugscore_min)
        self.volume_min = float(volume_min)
        self.volume_max = float(volume_max)
        self.parent_volume_max = float(parent_volume_max)
        self.plddt_min = float(plddt_min)
        self.rg_max = float(rg_max)
        self.tm_policy = tm_policy
        self.tm_overlap_max = float(tm_overlap_max)
        self.spacing = float(spacing)

    def as_dict(self) -> dict:
        return {k: v for k, v in vars(self).items()}


def _f(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def grid_legal(size, spacing: float = SPACING) -> bool:
    n = [math.floor(s / spacing) + 1 for s in size]
    return max(n) <= MAX_AXIS and n[0] * n[1] * n[2] <= MAX_POINTS


def tm_reject(is_af, tm_segments, residues, opts) -> bool:
    """Should this pocket be rejected by the transmembrane policy?

    A: drop proteins with >= 2 TM segments (multi-pass).
    B: drop any protein with >= 1 TM segment.
    C (default): A, plus drop a pocket whose lining residues overlap a TM segment by
       more than `tm_overlap_max`.

    The overlap check applies to **AlphaFold models only**, where the PDB residue number
    is the UniProt position. Experimental chains carry pdbfixer numbering, so the same
    comparison would be meaningless; for them only the protein-level rule applies, on
    the assumption that a crystallised construct of a single-pass protein is its soluble
    domain. That assumption belongs in Methods.

    Measured on net1, policy C costs 4 pockets over policy A, and policy B costs 21
    proteins including CAH9, whose druggable site is its extracellular catalytic domain.
    """
    n = len(tm_segments)
    if opts.tm_policy == "B":
        return n >= 1
    if n >= 2:
        return True
    if opts.tm_policy == "C" and is_af and n >= 1 and residues:
        covered = {r for a, b in tm_segments for r in range(a, b + 1)}
        return len(set(residues) & covered) / len(residues) > opts.tm_overlap_max
    return False


def classify(row, ann, tm_segments, parent_volume, opts) -> str | None:
    """None if the pocket is eligible, else the first step it fails.

    `row` is a DrugScore.csv row, `ann` its pocket_annotation.csv row.
    """
    if _f(row.get("drugScore"), -1.0) < opts.drugscore_min:
        return "drugscore"
    vol = _f(row.get("volume"))
    if vol is None or not (opts.volume_min <= vol <= opts.volume_max):
        return "volume"
    pv = _f(parent_volume)
    if pv is not None and pv > opts.parent_volume_max:
        return "subpocket_parent"
    is_af = str(ann.get("AlphaFold", "")).strip().lower() == "true"
    if is_af:
        pl = _f(ann.get("pocket_plddt_mean"))
        if pl is None or pl < opts.plddt_min:
            return "plddt"
    residues = ann.get("_residues") or []
    if tm_reject(is_af, tm_segments, residues, opts):
        return "tm"
    rg = _f(ann.get("radius_of_gyration"))
    if rg is None or rg > opts.rg_max:
        return "rg"
    size = [_f(ann.get(f"box_size_{k}"), 0.0) for k in "xyz"]
    if not grid_legal(size, opts.spacing):
        return "grid"
    return None


def parent_key(row):
    """The (accession, protein, parent pocket) a subpocket belongs to, or None for a
    top-level pocket. The structure is part of the key because `P_1` exists for every
    structure of an accession."""
    pocket = str(row.get("Pocket", ""))
    parts = pocket.split("_")
    if len(parts) <= 2:
        return None
    return (row["UniProt ID"], row.get("Protein", ""), "_".join(parts[:-1]))


def select_panel(rows, annotations, tm, parent_volumes, opts):
    """(panel, eligible, funnel, rejections).

    panel      {accession: (DrugScore row, annotation row)} — one pocket per protein,
               highest drugScore, ties to the smallest pocket name.
    eligible   [(DrugScore row, annotation row)] — every pocket that passed the filters,
               not just the one per protein. This is what `similarity` runs on: the
               post hoc test needs pairwise similarity among panel pockets, and the
               wider eligible set also leaves a descriptive network for the paper.
    funnel     [(step, pockets remaining after it)], a chain that sums correctly.
    rejections [(item, step)] for every rejected pocket, for failed.tsv.
    """
    eligible, rejections = [], []
    failed_at = {s: 0 for s in STEPS}
    for row in rows:
        acc = row["UniProt ID"]
        fname = row.get("Pocket Filename", "")
        ann = annotations.get(fname)
        if ann is None:
            rejections.append((f"{acc}/{fname}", "missing annotation"))
            continue
        pk = parent_key(row)
        step = classify(row, ann, tm.get(acc, []),
                        None if pk is None else parent_volumes.get(pk), opts)
        if step is None:
            eligible.append((row, ann))
        else:
            failed_at[step] += 1
            rejections.append((f"{acc}/{fname}", step))

    total = len(rows)
    funnel, remaining = [("selected", total)], total
    for s in STEPS:
        remaining -= failed_at[s]
        funnel.append((s, remaining))

    panel: dict[str, tuple] = {}
    for row, ann in sorted(eligible,
                           key=lambda ra: (-_f(ra[0].get("drugScore"), 0.0),
                                           str(ra[0].get("Pocket", "")))):
        panel.setdefault(row["UniProt ID"], (row, ann))
    return panel, eligible, funnel, rejections
