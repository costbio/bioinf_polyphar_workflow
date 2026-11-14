"""Panel selection: the eligible-pocket filters and one-row-per-protein rule.

This replaces hub-and-neighbours selection. The hub was the highest-degree node of the
pocket-similarity graph, which made the target count a property of how many structures
happened to be processed rather than of biology, and anchored everything on one pocket
chosen by network topology. Pocket similarity is no longer used to select targets at
all; it is applied post hoc to the docking results.

Choosing the filters
--------------------
An earlier version filtered on radius of gyration and on a count of sequence segments.
Both were invented thresholds laid on top of DoGSiteScorer's trained druggability score,
and neither was checked against anything. Calibrating them against pockets holding a
co-crystallised drug-like ligand (687 of 8,005 experimental net1 pockets) showed
`Rg <= 8.5` kept only **8%** of pockets known to bind a drug while keeping 69% of the
rest — anti-enriching roughly eightfold, because the lining-residue envelope of a real
binding site is large, not small. Both were removed, along with a separate box cap.

The filters below are chosen on physical grounds and stated plainly, rather than tuned to
net1, which is one subtype with a known-bad history.

  drugScore          DoGSiteScorer's trained druggability SVM: the field-standard answer
                     to this question. NOTE: DrugScore.csv's `Selected` column already
                     means drugScore >= 0.5, so at the default this step is inherited
                     from `pockets` and constrains nothing further. Raise it above 0.5
                     for it to do work.
  volume 300-1500    the cited druggable range is ~300-1000 A^3, with stricter criteria
                     near 1000; 1500 is a deliberate allowance above it. The ceiling is
                     what rejects the oversized interdomain cavities that dominated the
                     old network — those run ~2,600 A^3 against ~270 for a typical
                     pocket. At this setting 2 of the 100 worst survive.
  subpocket parent   without it, a subpocket carved from a 25,840 A^3 cleft passes the
                     ceiling and becomes its protein's representative.
  AF pocket pLDDT    a sanity floor only. pLDDT does NOT detect the oversized cavities
                     (the worst average 84.8 against 83.3 overall): it is local
                     per-residue confidence, not domain packing. PAE would be the right
                     metric.
  (surf/vol)         REMOVED. It was added as an enclosure criterion, then measured to
                     have a 25.9-point AlphaFold-versus-experimental pass-rate gap
                     (12.2% against 38.1%) — worse than the 21.9 that got
                     `segments <= 3` rejected. The AlphaFold share falling 53% -> 32%
                     under it was that bias, not an improvement. Its enrichment claim
                     was also confounded: surface-to-volume falls with volume by
                     geometry, so "ligand-bound pockets have lower surf/vol" may simply
                     restate "ligand-bound pockets are bigger". Third filter in a row to
                     fail the same class-bias check; removed (user, 2026-11-03).
  TM policy          see tm_reject(). The artifact is multi-pass proteins, whose single
                     extracted chains yield lipid-facing grooves and pore lumens.
  Uni-Dock grid      a guard, not a selector: the engine exits on a violation, and
                     versions before 1.2.0 silently corrupted results instead.

One row per protein: the eligible pocket with the highest drugScore, ties broken by the
smallest pocket name so the result never depends on row order.
"""

from __future__ import annotations

import math

DRUGSCORE_MIN = 0.5
VOLUME_MIN = 300.0
VOLUME_MAX = 1500.0
PARENT_VOLUME_MAX = 1500.0
PLDDT_MIN = 70.0
TM_OVERLAP_MAX = 0.10
SPACING = 0.375
MAX_AXIS = 128
MAX_POINTS = 531441       # what Uni-Dock formally refuses: 81^3
# What it can actually WORK. Uni-Dock accepts a box just under its own limit and then
# crawls: on net3 a unit at ~530k points took 2 h 17 min against the ~10 min a normal one
# takes, so a single such unit consumed a whole Colab session and the rest of the bundle
# was never drawn. `incomplete_targets` could therefore never reach zero and "re-run until
# nothing is outstanding" never terminated. Measured over 212 net3 targets: of the 197 that
# COMPLETED, the largest grid was 508,518 points; of the 15 that never finished, 7 were
# above it. The separation is clean, and the 4 % margin below MAX_POINTS is where it sits.
#
# The performance is the symptom. These boxes are 26-40 A per axis, against 10-15 A for a
# drug site: the pockets have twice the lining residues of a normal one (46 vs 25) in as
# many as 8 disconnected sequence segments, so they are surface grooves, not cavities. The
# volume filter does not catch them because DoGSite's volume is the CAVITY (1,028 A^3
# median, inside the 300-1500 range) while the docking box is the bounding box of the
# lining residues -- a different quantity, and until now an unconstrained one. A
# within-pocket percentile over a 30 A cube does not mean what the statistic assumes.
WORKABLE_POINTS = 508518
TM_POLICIES = ("A", "B", "C")
PER_PROTEIN = 2   # pockets docked per protein, best drugScore first

# The funnel steps, in order. Each pocket is attributed to the FIRST step it fails, so
# the counts form a chain and sum correctly.
STEPS = ["drugscore", "volume", "subpocket_parent", "plddt", "tm", "grid"]
STEP_REASON = {
    "drugscore": "drugScore below the cutoff",
    "volume": "pocket volume outside the allowed range",
    "subpocket_parent": "subpocket of a parent above the volume ceiling",
    "plddt": "AlphaFold pocket pLDDT below the floor",
    "tm": "rejected by the transmembrane policy",
    "grid": "docking box too large to dock (over the workable grid ceiling)",
}


class Options:
    """Filter thresholds. Defaults are the decided panel definition."""

    def __init__(self, drugscore_min=DRUGSCORE_MIN, volume_min=VOLUME_MIN,
                 volume_max=VOLUME_MAX, parent_volume_max=PARENT_VOLUME_MAX,
                 plddt_min=PLDDT_MIN, tm_policy="C",
                 tm_overlap_max=TM_OVERLAP_MAX, spacing=SPACING,
                 per_protein=PER_PROTEIN):
        if tm_policy not in TM_POLICIES:
            raise ValueError(f"tm_policy must be one of {TM_POLICIES}, got {tm_policy!r}")
        self.drugscore_min = float(drugscore_min)
        self.volume_min = float(volume_min)
        self.volume_max = float(volume_max)
        self.parent_volume_max = float(parent_volume_max)
        self.plddt_min = float(plddt_min)
        self.tm_policy = tm_policy
        self.tm_overlap_max = float(tm_overlap_max)
        self.spacing = float(spacing)
        # More than one pocket per protein because picking one is a coin flip: measured
        # against pockets with a co-crystal ligand, the highest-drugScore pick was the
        # known site in only 17 of 29 cases (59%). Docking the best few removes the
        # choice instead of trying to make it correctly. Only pockets that pass every
        # filter are taken, so a protein with one eligible pocket contributes one.
        self.per_protein = max(1, int(per_protein))

    def as_dict(self) -> dict:
        return {k: v for k, v in vars(self).items()}


def _f(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def grid_legal(size, spacing: float = SPACING, limit: int = WORKABLE_POINTS) -> bool:
    """Whether a box of this size can be docked in practice, not merely accepted.

    `limit` defaults to WORKABLE_POINTS, not MAX_POINTS: a box between the two is legal to
    Uni-Dock and unusable to us. Pass MAX_POINTS to ask the engine's own question.
    """
    n = [math.floor(s / spacing) + 1 for s in size]
    return max(n) <= MAX_AXIS and n[0] * n[1] * n[2] <= limit


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
    size = [_f(ann.get(f"box_size_{k}")) for k in "xyz"]
    # a missing dimension must FAIL: defaulting to 0.0 made a malformed annotation into
    # a 1-point-per-axis box, which is "legal" and silently becomes a dockable target
    if any(s is None or s <= 0 for s in size) or not grid_legal(size, opts.spacing):
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

    panel      {accession: [(DrugScore row, annotation row), ...]} — up to
               `opts.per_protein` pockets per protein, highest drugScore first, ties to
               the smallest pocket name. `annotations` is keyed by (accession, filename).
    eligible   [(DrugScore row, annotation row)] — every pocket that passed the filters,
               not just the one per protein. This is what `similarity` runs on: the
               post hoc test needs pairwise similarity among panel pockets, and the
               wider eligible set also leaves a descriptive network for the paper.
    funnel     [(step, pockets remaining after it)], a chain that sums correctly.
    rejections [(item, step)] for every rejected pocket, for failed.tsv.
    """
    eligible, rejections = [], []
    failed_at = {s: 0 for s in STEPS}
    missing_annotation = 0
    for row in rows:
        acc = row["UniProt ID"]
        fname = row.get("Pocket Filename", "")
        # keyed on the PAIR: a chimeric or shared chain can give the same pocket file
        # name under two accessions, and a filename-only key silently overwrites one
        # annotation with the other
        ann = annotations.get((acc, fname))
        if ann is None:
            missing_annotation += 1
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
    funnel = [("selected", total)]
    remaining = total - missing_annotation
    if missing_annotation:
        funnel.append(("annotation", remaining))
    for s in STEPS:
        remaining -= failed_at[s]
        funnel.append((s, remaining))
    assert remaining == len(eligible), "funnel must end at the eligible count"

    panel: dict[str, list] = {}
    for row, ann in sorted(eligible,
                           key=lambda ra: (-_f(ra[0].get("drugScore"), 0.0),
                                           str(ra[0].get("Pocket", "")))):
        picks = panel.setdefault(row["UniProt ID"], [])
        if len(picks) < opts.per_protein:
            picks.append((row, ann))
    return panel, eligible, funnel, rejections
