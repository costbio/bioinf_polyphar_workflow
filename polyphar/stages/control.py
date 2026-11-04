"""Stage 08: the positive control -- does a docking score mean anything at all?

Implements item 5 of `docs/superpowers/specs/2026-11-03-statistics-design.md`
("Positive control: does the score mean anything at all?"). `polyphar stats` already
CONSUMES this stage's result: it reads `inputs/positive_control.json`, range-checks the
two readouts and stamps every primary row `yes` / `no: ... UNINTERPRETABLE` /
`unknown: not run`. This stage is what PRODUCES that file.

Why it exists
    The receptors are ~63% AlphaFold models and the engine carries an uncharacterised
    +0.33 kcal/mol offset from smina. Nothing else in the pipeline anchors a score to real
    binding, and the usual anchor -- the approved HCC drugs' own targets -- is unavailable,
    because the subtype gene lists do not contain them. What IS available is the panel's
    own crystallography: a pocket with a drug-like ligand sitting in it is a site a
    molecule is KNOWN to bind, in this very dataset, at these very coordinates.

The three readouts
    1. Self-recovery rank. Every labelled ligand is docked against the WHOLE panel. Its
       own pocket should rank far above chance (chance = 0.5). Reported as a percentile
       in [0, 1] so it is directly comparable with stats' `SELF_RANK_MAX = 0.05`, plus the
       fraction of ligands landing in the top 1% / 5%.
    2. Redocking RMSD. For the SELF pairs -- the labelled ligand against the pocket it was
       crystallised in -- the symmetry-insensitive RMSD between replicate 0's best pose and
       the deposited coordinates, gated at the field-standard `REDOCK_RMSD_MAX = 2.0 A`.
       This replaces pose-in-pocket as a gate (2026-11-04): `dock.pocket_box` is the
       UNPADDED bounding box of the pocket's own lining atoms and the engine keeps the
       ligand inside it, so "centroid inside the box, away from the faces" nearly restates
       a constraint the engine already enforces. Measured on real net2 pockets, random
       NON-CLASHING centroids confined to that box pass pose-in-pocket about 79% of the
       time against 85% for the true deposited poses -- a gate at 0.70 sits BELOW its own
       null. The same random placements pass RMSD <= 2 A 0.0% of the time. See CAVEAT.
    3. Library percentile at home. The co-crystal ligand's score percentile among the
       PRODUCTION library at its own target, from the same receptor, the same box and the
       same engine (`07_docking/df_results_<home>.csv` beside
       `08_control/docking/df_results_<home>.csv`). Readout 1 is within-ligand across
       targets; `M` is within-target across ligands. Those are rows and columns of the same
       matrix and one does not imply the other, so the axis `M` is built on is validated
       directly: a known binder that does not reach the top decile of the library at the
       pocket it was crystallised in means `M` carries no site-discriminative content.

    Pose-in-pocket is still computed and reported for every (ligand, target) pair, as a
    DESCRIPTIVE diagnostic. It is not a gate.

The label, and what it is NOT
    `docs/superpowers/specs/2026-11-03-panel-selection-spikes/rg_calibrate.py` labelled 687
    experimental pockets this way for the Rg calibration; its criterion is reused verbatim
    here (SKIP, CONTACT, MIN_FRAC, MIN_HEAVY, MAX_HEAVY) rather than reinvented. Two
    deliberate differences: the pockets labelled are the PANEL's targets
    (06_targets/final_targets.csv), because only a pocket in the panel can have a
    self-rank at all; and the extracted chain is found through the target's own
    `Receptor_Path` (03_pdbfixer/<acc>/<file> -> 02_structures/<acc>/<file>) instead of by
    splitting the pocket file name on "_res_P", which is a DoGSiteScorer-2-only spelling
    and silently misses chimeric chains, whose pocket files carry an extra `<acc>_` prefix.

    Coordinates line up because pdbfixer leaves existing heavy atoms where they are:
    `pockets.py` `_author_atoms`/`_lookup` map pocket residues back to the extracted chain
    within 0.001 A, and the claim was re-checked for this stage (2026-11-04) on
    tests/fixtures -- all 788 heavy atoms of 2mjw_A and all 719 of 2mjw_B's first model
    come back through `fix_one` byte-identical in the PDB's 3-decimal representation, so
    the 0.001 A tolerance is not merely satisfied, it is never exercised.

    CAVEAT (it travels with the output, in the JSON and in summary.md, because the number
    will be quoted by people who never read this file): see `CAVEAT`.

How it reuses `dock`
    There is ONE docking path in this repo and this stage does not add a second. It builds
    a ligand SDF, then calls `dock.units` / `dock.run` with
      * `library=` -- the co-crystal SDF, docked as given (`dock._library_path`, and with
        it the `ligands` stage's prepared e-Drug3D, is bypassed), and
      * a `ControlLayout` -- the same run folder with `docking` redirected to
        08_control/docking, so the production results in 07_docking are untouched.
    Receptor preparation, boxes, unit keys, resume, discards, provenance, failed.tsv and
    the per-target summaries are therefore literally dock's, and a (target, replicate)
    here is the same kind of object as a (target, replicate) there.

    `--poses` is forced to `first` (replicate 0's poses are what readout 2 reads) and the
    replicate count is this stage's own `--control-replicates` (default 1): the control
    asks a different question from the production run and does not need its jitter budget.

Outputs in 08_control/
    labels.csv              one row per panel target: its co-crystal ligand, or why not
    docking/                dock's own output tree for the co-crystal library
    self_rank.csv           per ligand: home target, self-rank percentile, n panel targets
    pose_in_pocket.csv      per (ligand, target): centroid, verdict, why, self-pair RMSD
    library_percentile.csv  per ligand: its score percentile among the production library
                            at its own target
    positive_control.json   THE RESULT (also installed at inputs/positive_control.json)
    summary.md              the three gated numbers with the caveat attached
    failed.tsv, manifest.json (written last = done)

Determinism: labelling, the SDF and both analyses run in the parent process over sorted
inputs, and dock's own parallel part writes one file per unit, so `--nprocs` cannot change
a byte of the output. Resume: `units()` returns 0 once the manifest's fingerprint matches.
"""

from __future__ import annotations

import csv
import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from polyphar.common import RunLayout, WorkUnits, append_failure, sha256_file, write_manifest
from polyphar.stages import dock

STAGE = "control"
DIRNAME = "08_control"
DOCKING_DIRNAME = "docking"
# Inside <run>/inputs so `dock._raw_library_path` never tries to copy it anywhere (and so
# the ligand set that produced the control sits beside the positive set that `stats` reads).
LIBRARY = "inputs/control_cocrystal_ligands.sdf"
POSITIVE_CONTROL = "positive_control.json"


def _inputs_positive_control() -> str:
    """What `stats` reads by default -- imported FROM stats, not spelled again: two string
    literals that must stay equal are one literal in the module that owns the contract."""
    from polyphar.stages import stats

    return stats.POSITIVE_CONTROL


INPUTS_POSITIVE_CONTROL = _inputs_positive_control()
DEFAULT_REPLICATES = 1
# The molblock's program line, written in place of Open Babel's MMDDYYHHmm stamp: the SDF
# must be a function of the coordinates alone, or the library hash (and with it every unit
# key) changes at the next minute boundary.
PROGRAM_LINE = " polyphar-control      3D"

# --- the co-crystal label (rg_calibrate.py's criterion, unchanged) ----------------------
CONTACT = 4.5       # A, ligand heavy atom to pocket lining atom
MIN_FRAC = 0.70     # share of the ligand's heavy atoms that must be in contact
MIN_HEAVY = 12      # drug-like size range, heavy atoms
MAX_HEAVY = 60

# water, ions, buffers, cryoprotectants, glycans, lipids and detergents: present in
# crystals but not evidence of a drug-binding site. Verbatim from rg_calibrate.py.
SKIP = set("""HOH DOD UNX UNK
NAG MAN BMA FUC GAL GLC NDG SIA XYS RAM BGC A2G
GOL EDO PEG PGE PG4 1PE 2PE P6G MPD DMS DMF TRS ACT ACY FMT OXL
SO4 PO4 NO3 CO3 EPE MES CIT TAR FLC IMD BME DTT TCE MLI
CL NA K MG CA ZN FE FE2 MN CU CU1 NI CD HG IOD BR CS SR BA CO PT AU AG
CLR POV PLM OLA MYR STE LDA C8E BOG DDQ LMT PC1 PEE PCW HEM HEC""".split())

# --- pose-in-pocket geometry ------------------------------------------------------------
# "Inside the cavity, not a surface patch at the box edge" made explicit. These three
# numbers define what the reported fraction MEANS, so they are recorded in the JSON.
SHELL_RADIUS = 8.0      # A: pocket lining atoms that count as surrounding the centroid
MIN_LINING_ATOMS = 4    # fewer than this within the shell = open space, not a cavity
# ||mean(lining atoms in the shell) - centroid|| / SHELL_RADIUS. 0 = the centroid sits at
# the middle of its lining; a pose lying ON a surface patch has all of its lining on one
# side and scores ~0.4-0.5. The cut is deliberately lenient: this is a sanity check.
MAX_ASYMMETRY = 0.5
EDGE_MARGIN = 2.0       # A the centroid must stay inside every face of the docking box,
#                         capped at a quarter of that axis so a thin box is not impossible

POSE_IN_POCKET_DEFINITION = (
    f"DESCRIPTIVE DIAGNOSTIC, NOT A GATE. Centroid of the best replicate-0 pose: "
    f">= {MIN_LINING_ATOMS} pocket lining atoms within {SHELL_RADIUS:g} A, their mean "
    f"displaced from the centroid by at most {MAX_ASYMMETRY:g} x {SHELL_RADIUS:g} A (so "
    f"the lining surrounds it rather than lying to one side of it), and the centroid at "
    f"least min({EDGE_MARGIN:g} A, a quarter of the axis) inside every face of the "
    f"unpadded docking box. Its null is high -- random NON-CLASHING centroids confined to "
    f"the same box pass about 79% of the time against 85% for the true deposited poses on "
    f"real net2 pockets -- because the box IS the pocket's own unpadded bounding box and "
    f"the engine already keeps the ligand inside it. Computed for every (ligand, target) "
    f"pair and most informative for the NON-COGNATE ones, which have no deposited pose to "
    f"measure an RMSD against; for the self pairs the RMSD above is the readout that "
    f"matters")

# --- redocking RMSD (readout 2, the gate) ------------------------------------------------
REDOCK_RMSD_MAX = 2.0      # A; the field-standard redocking-success criterion
REDOCK_RMSD_DEFINITION = (
    f"symmetry-insensitive RMSD, in place and with no superposition, between replicate 0's "
    f"best-scoring pose and the ligand's DEPOSITED coordinates, over the SELF pairs only "
    f"(the labelled ligand against the pocket it was crystallised in). Heavy atoms are "
    f"matched by optimal assignment WITHIN EACH ELEMENT rather than by file order, so a "
    f"symmetric group has no single correct correspondence to get wrong; that is <= the "
    f"order-preserving RMSD, so a pose is never reported as worse than it is. Success is "
    f"RMSD <= {REDOCK_RMSD_MAX:g} A -- a real redocking-success rate with a published base "
    f"rate (Vina-family cognate redocking is usually reported at 60-80%), whose null is "
    f"0.0%: random non-clashing rigid placements of the same ligand in the same box pass "
    f"0 of 1,871 draws on real net2 pockets")

# --- library percentile at home (readout 3, the gate on the axis M is built on) ----------
TOP_DECILE = 0.10
LIBRARY_PERCENTILE_DEFINITION = (
    f"mid-ranked percentile of the co-crystal ligand's own score among the PRODUCTION "
    f"library's scores at the SAME target (07_docking/df_results_<home>.csv), lower = "
    f"better, so chance is 0.5 and the top decile is <= {TOP_DECILE:g}. Same receptor "
    f"PDBQT, same box, same engine, same replicate policy -- only the ligand set differs. "
    f"This is the axis `M` is built on (within-target, across-ligand); readout 1 is its "
    f"transpose (within-ligand, across-target) and does not imply it")

CAVEAT = (
    "SANITY CHECK, NOT AN OPTIMISATION TARGET. "
    "(1) THE EASIEST BENCHMARK IN THE FIELD, and this is the largest bias here: the home "
    "pocket is the HOLO conformation, induced-fit around this very ligand. DoGSiteScorer "
    "ran on the pdbfixer structure, i.e. the holo backbone with the heterogens stripped, "
    "so self-docking recovers a pose into a receptor whose side chains are already "
    "arranged for it. A passing self-recovery rank and a passing redocking RMSD say the "
    "engine and the pocket are not broken; they say nothing about cross-docking, about "
    "the ~63% of the panel that is AlphaFold, or about absolute affinity. "
    "(2) The co-crystal label says a drug-like ligand was observed in that pocket; its "
    "absence says nothing about druggability, because most such structures were solved "
    "apo. The label is biased toward LARGE pockets (a 12-60 heavy-atom ligand needs room) "
    "and toward WELL-STUDIED proteins (crystallographers solve complexes for sites they "
    "already care about). "
    "(3) \"DRUG-LIKE LIGAND\" OVERSTATES THE CRITERION, which is only a 12-60 heavy-atom "
    "non-solvent HETATM residue in contact with the lining. Across 1,142 net2 chains the "
    "commonest labelled codes are PSU 49 (pseudouridine), 7PE 22 (a PEG oligomer the SKIP "
    "list misses), NAP 20, OMG 19, FAD 18, ATP 12, ADP 10: about 14% modified RNA residues "
    "and about 10% nucleotide cofactors, neither of which is a drug. "
    "(4) WHAT A PASS LICENSES, exactly: within THIS panel and THIS engine, a molecule "
    "scores better at the site it is known to occupy than at ~95% of other panel sites, "
    "recovers that site's deposited pose, and outranks the production library there -- so "
    "a within-target percentile carries site-discriminative content. Nothing about "
    "absolute affinity, and nothing that calibrates the threshold x min_targets surface. "
    "In particular the self-rank is a relative within-ligand comparison and is therefore "
    "INVARIANT to any monotone transform applied uniformly across targets, including the "
    "+0.33 kcal/mol smina offset, which cancels exactly: the control cannot detect the "
    "defect its own rationale names. "
    "Nothing in the pipeline may be tuned to raise these numbers: they exist to say "
    "whether a score and a pocket mean anything at all, and a tuned control answers "
    "nothing.")

LABEL_COLUMNS = ["Target", "UniProt_ID", "Pocket", "Ligand_ID", "Ligand", "Chain",
                 "Residue", "Heavy_Atoms", "Contact_Fraction", "Redocked", "Reason"]
SELF_RANK_COLUMNS = ["Ligand_ID", "Ligand", "Home_Target", "Home_Score", "Best_Score",
                     "Panel_Targets", "Self_Rank", "Self_Rank_Percentile"]
POSE_COLUMNS = ["Target", "Ligand_ID", "Score", "Centroid_X", "Centroid_Y", "Centroid_Z",
                "Lining_Atoms", "Asymmetry", "In_Pocket", "Contact_Fraction",
                "Self_Pair", "RMSD", "Redock_Success", "Reason"]
LIBRARY_COLUMNS = ["Ligand_ID", "Ligand", "Home_Target", "Control_Score",
                   "Library_Ligands", "Better_In_Library", "Percentile", "Top_Decile",
                   "Reason"]


# --- small geometry helpers --------------------------------------------------------------

def _dist2(a, b) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def centroid(xyz) -> tuple:
    n = len(xyz)
    return tuple(sum(p[k] for p in xyz) / n for k in range(3))


def contact_fraction(ligand_xyz, pocket_xyz, cutoff: float = CONTACT) -> float:
    """Share of the ligand's heavy atoms within `cutoff` of any pocket atom."""
    if not ligand_xyz:
        return 0.0
    c2 = cutoff * cutoff
    hit = sum(1 for a in ligand_xyz if any(_dist2(a, b) <= c2 for b in pocket_xyz))
    return hit / len(ligand_xyz)


def pose_in_pocket(pose_centroid, pocket_xyz, box_center, box_size) -> tuple:
    """(verdict, lining atom count, asymmetry, reason) for one pose centroid.

    See POSE_IN_POCKET_DEFINITION. The reason is "" for a pose in the pocket and names the
    failing condition otherwise, so pose_in_pocket.csv says WHY a pose was counted out."""
    near = [p for p in pocket_xyz if _dist2(p, pose_centroid) <= SHELL_RADIUS ** 2]
    if len(near) < MIN_LINING_ATOMS:
        return False, len(near), float("nan"), (
            f"only {len(near)} lining atom(s) within {SHELL_RADIUS:g} A: open space, "
            f"not a cavity")
    asym = math.sqrt(_dist2(centroid(near), pose_centroid)) / SHELL_RADIUS
    if asym > MAX_ASYMMETRY:
        return False, len(near), asym, (
            f"lining asymmetry {asym:.2f} > {MAX_ASYMMETRY:g}: the pocket atoms lie to one "
            f"side of the centroid, i.e. a surface patch")
    for k, axis in enumerate("xyz"):
        inset = min(EDGE_MARGIN, 0.25 * box_size[k])
        lo = box_center[k] - box_size[k] / 2 + inset
        hi = box_center[k] + box_size[k] / 2 - inset
        if not lo <= pose_centroid[k] <= hi:
            return False, len(near), asym, (
                f"centroid {axis}={pose_centroid[k]:.2f} is within {inset:.2f} A of the "
                f"box edge [{lo - inset:.2f}, {hi + inset:.2f}]")
    return True, len(near), asym, ""


# --- symmetry-insensitive RMSD (readout 2) -----------------------------------------------

def _hungarian(cost) -> list:
    """Minimum-cost perfect assignment on a square cost matrix (Jonker-Volgenant
    shortest-augmenting-path form, O(n^3)). Returns the column chosen for each row.

    Written out rather than taken from scipy: scipy is not in environment.yml, and a
    50-line exact algorithm is preferable to a new dependency in the image."""
    n = len(cost)
    if n == 0:
        return []
    inf = float("inf")
    u = [0.0] * (n + 1)
    v = [0.0] * (n + 1)
    p = [0] * (n + 1)          # p[j] = the row currently assigned to column j
    way = [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [inf] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], inf, 0
            for j in range(1, n + 1):
                if used[j]:
                    continue
                cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                if cur < minv[j]:
                    minv[j], way[j] = cur, j0
                if minv[j] < delta:
                    delta, j1 = minv[j], j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
    out = [0] * n
    for j in range(1, n + 1):
        out[p[j] - 1] = j - 1
    return out


def _by_element(atoms) -> dict:
    """{element: [coordinates]} for (x, y, z, element) atoms, element upper-cased."""
    out: dict = {}
    for a in atoms:
        out.setdefault(str(a[3]).strip().upper(), []).append(a)
    return out


def symmetry_insensitive_rmsd(pose, deposited) -> float:
    """RMSD in A between a docked pose and the deposited ligand. See REDOCK_RMSD_DEFINITION.

    In place: NO superposition is applied, because the question is whether the engine put
    the ligand back where the crystallographer found it, in the receptor's own frame.
    Raises ValueError when the two are not the same heavy-atom composition, which means the
    pose is not of this ligand and no RMSD is defined."""
    a, b = _by_element(pose), _by_element(deposited)
    if {k: len(v) for k, v in a.items()} != {k: len(v) for k, v in b.items()}:
        raise ValueError(
            f"composition differs: pose {sorted((k, len(v)) for k, v in a.items())} vs "
            f"deposited {sorted((k, len(v)) for k, v in b.items())}")
    total, n = 0.0, 0
    for el, pa in sorted(a.items()):
        pb = b[el]
        cost = [[_dist2(x, y) for y in pb] for x in pa]
        for i, j in enumerate(_hungarian(cost)):
            total += cost[i][j]
            n += 1
    if not n:
        raise ValueError("no heavy atoms")
    return math.sqrt(total / n)


def self_rank_percentile(home, scores: dict) -> float:
    """Percentile of `home` among `scores` ({target: score}, lower = better), mid-ranked.

    The denominator is len(scores) -- THE WHOLE PANEL the ligand was docked against, never
    a subset: a ligand missing from any panel target is dropped by `self_ranks` rather than
    ranked against the targets that happen to have a number. The best of N comes out at
    0.5/N, chance at ~0.5, the worst at ~1."""
    if home not in scores:
        raise KeyError(f"{home!r} is not among the {len(scores)} scored targets")
    n = len(scores)
    s = scores[home]
    # NaN compares false against everything, so a non-finite HOME score would come out
    # `better = 0, ties = 1` -> percentile 0.0, BETTER than the best achievable 0.5/N, and
    # an all-NaN matrix would read as a perfect control. `_summary_scores` already drops
    # non-finite rows through `stats.read_df_results`; this is the assertion that says so.
    bad = sorted(t for t, v in scores.items() if not math.isfinite(v))
    if bad:
        raise ValueError(f"non-finite score(s) at {', '.join(bad[:3])}"
                         f"{'...' if len(bad) > 3 else ''}: a self-rank over them would "
                         f"put the home target above the best achievable rank")
    better = sum(1 for v in scores.values() if v < s)
    ties = sum(1 for v in scores.values() if v == s)
    midrank = better + (ties + 1) / 2
    return (midrank - 0.5) / n


def median(values) -> float:
    v = sorted(values)
    n = len(v)
    if not n:
        raise ValueError("median of no values")
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


# --- PDB reading --------------------------------------------------------------------------

def _xyz(line) -> tuple | None:
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except ValueError:
        return None


def _is_heavy(line) -> bool:
    el = (line[76:78].strip() or line[12:14].strip()).upper()
    return not el.startswith(("H", "D"))


def het_groups(path) -> dict:
    """{(resname, chain, residue): [heavy-atom coordinates]} of the drug-like-candidate
    HETATM residues of a PDB file's FIRST model. SKIP residues are left out here, size is
    checked by the caller (the counts belong in labels.csv)."""
    out: dict = {}
    try:
        fh = open(path, encoding="utf-8", errors="replace")
    except OSError:
        return {}
    with fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("HETATM") or not _is_heavy(line):
                continue
            resn = line[17:20].strip()
            if resn in SKIP:
                continue
            c = _xyz(line)
            if c is not None:
                out.setdefault((resn, line[21], line[22:27].strip()), []).append(c)
    return out


def _element(line) -> str:
    el = (line[76:78].strip() or line[12:14].strip()).upper()
    return "".join(ch for ch in el if ch.isalpha()) or "X"


def het_atoms(path, key) -> list:
    """[(x, y, z, element)] for the heavy atoms of one HETATM residue of a PDB file's first
    model: the DEPOSITED coordinates readout 2 measures the redocked pose against."""
    out = []
    for line in het_lines(path, key):
        if not _is_heavy(line):
            continue
        c = _xyz(line)
        if c is not None:
            out.append((c[0], c[1], c[2], _element(line)))
    return out


def het_lines(path, key) -> list:
    """The verbatim HETATM lines of one residue of a PDB file's first model, hydrogens
    included (Open Babel reads what the crystallographer deposited)."""
    out = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if (line.startswith("HETATM")
                    and (line[17:20].strip(), line[21], line[22:27].strip()) == key):
                out.append(line.rstrip("\n"))
    return out


def pocket_atoms(path, hetatm: bool = False) -> list:
    """Heavy-atom coordinates of a DoGSiteScorer pocket residue file.

    ATOM only by default -- the protein lining, which is what rg_calibrate.py contacted the
    ligand against, so the label here is the same label. `hetatm=True` adds HETATM, which
    is the atom set `dock.pocket_box` builds the docking box from, so the geometry of
    readout 2 is measured against the same atoms that defined the box."""
    out = []
    tags = ("ATOM", "HETATM") if hetatm else ("ATOM",)
    try:
        fh = open(path, encoding="utf-8", errors="replace")
    except OSError:
        return []
    with fh:
        for line in fh:
            if line.startswith(tags) and _is_heavy(line):
                c = _xyz(line)
                if c is not None:
                    out.append(c)
    return out


# --- labelling ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Label:
    target: str            # dock's target label, f"{Protein}_{Pocket}"
    accession: str
    pocket_rel: str
    ligand_id: str         # the SDF <ID>; unique per (pocket, ligand residue)
    resname: str
    chain: str
    residue: str
    heavy_atoms: int
    contact: float
    chain_pdb: Path
    key: tuple


def chain_pdb_of(layout, target) -> Path | None:
    """The EXTRACTED CHAIN behind a target's receptor, which still has its HETATM records
    (pdbfixer strips heterogens, so 03_pdbfixer is no use here), or None for an AlphaFold
    model -- a prediction has no co-crystal ligand to find.

    Derived from Receptor_Path, not from the pocket file name: `pocket_filename` spells
    DoGSite3 output differently from DoGSiteScorer 2, and a chimeric chain's pocket file
    carries an extra `<accession>_` prefix that the name-splitting approach mistakes for
    part of the chain's file name."""
    if target.alphafold:
        return None
    parts = Path(target.receptor_rel).parts
    if len(parts) < 2 or parts[0] != layout.pdbfixer.name:
        return None
    return layout.structures.joinpath(*parts[1:])


def _ligand_id(accession: str, target: str, resname: str, chain: str, residue: str) -> str:
    stem = f"{accession}_{target}_{resname}_{chain or 'x'}{residue}"
    return "".join(ch if (ch.isalnum() or ch in "-_") else "-" for ch in stem)


def label_target(layout, target) -> tuple:
    """(Label or None, reason it has none). The reason is "" for a labelled target."""
    chain = chain_pdb_of(layout, target)
    if chain is None:
        return None, "AlphaFold model: no co-crystal ligand exists"
    if not chain.is_file():
        return None, f"extracted chain {dock._rel(chain, layout.run_dir)} not found"
    pocket = layout.run_dir / target.pocket_rel
    lining = pocket_atoms(pocket)
    if not lining:
        return None, f"pocket {target.pocket_rel} has no heavy ATOM records"
    groups = het_groups(chain)
    if not groups:
        return None, "no non-solvent HETATM residue in the extracted chain"
    sized = {k: v for k, v in groups.items() if MIN_HEAVY <= len(v) <= MAX_HEAVY}
    if not sized:
        return None, (f"no HETATM residue of {MIN_HEAVY}-{MAX_HEAVY} heavy atoms "
                      f"(largest: {max(len(v) for v in groups.values())})")
    best, best_frac = None, 0.0
    for key in sorted(sized):               # sorted: the pick cannot depend on file order
        frac = contact_fraction(sized[key], lining)
        if frac > best_frac:
            best, best_frac = key, frac
    if best is None or best_frac < MIN_FRAC:
        return None, (f"no ligand with >= {MIN_FRAC:.0%} of its heavy atoms within "
                      f"{CONTACT:g} A of the pocket (best {best_frac:.0%})")
    resname, ch, residue = best
    acc = _accession_of(target.pocket_rel)
    return Label(target.label, acc, target.pocket_rel,
                 _ligand_id(acc, target.label, resname, ch, residue),
                 resname, ch, residue, len(sized[best]), best_frac, chain, best), ""


def _accession_of(pocket_rel: str) -> str:
    """The accession folder of `04_pockets/<accession>/<pocket file>`."""
    parts = Path(pocket_rel).parts
    return parts[1] if len(parts) > 1 else ""


def _physical(label) -> tuple:
    """What makes two Labels the SAME physical ligand: one residue of one extracted chain."""
    return (str(label.chain_pdb), label.key)


def label_targets(layout, targets) -> tuple:
    """([Label], [(target label, reason)]) over every panel target, in target order.

    DE-DUPLICATED PER PHYSICAL LIGAND. DoGSiteScorer reports nested subpockets, and a
    ligand filling a parent pocket usually also fills its subpocket (verified: ADP reaches
    contact 1.00 against both `_res_P_1` and `_res_P_1_1`). Keyed on
    (accession, target, resname, chain, residue), `_ligand_id` would give that one molecule
    two ligand IDs, two SDF records, two self-ranks and two redocking RMSDs -- n inflated
    with no information added, and the nested pair correlated almost perfectly. One row per
    physical ligand, the pocket it is most completely enclosed by (highest contact
    fraction, ties to the lexicographically first target) winning; the others are reported
    as unlabelled with the pocket that took them."""
    best: dict = {}
    reasons = []
    for t in sorted(targets, key=lambda t: t.label):
        label, why = label_target(layout, t)
        if label is None:
            reasons.append((t.label, why))
            continue
        key = _physical(label)
        kept = best.get(key)
        if kept is None or label.contact > kept.contact:
            if kept is not None:
                reasons.append((kept.target, _DUPLICATE.format(
                    res=kept.resname, other=label.target, frac=label.contact)))
            best[key] = label
        else:
            reasons.append((t.label, _DUPLICATE.format(
                res=label.resname, other=kept.target, frac=kept.contact)))
    labels = sorted(best.values(), key=lambda x: x.target)
    return labels, sorted(reasons)


_DUPLICATE = ("the same physical ligand {res} is already labelled on {other} (contact "
              "{frac:.0%}), a nested pocket of the same chain: counted once, there")


# --- the ligand SDF -----------------------------------------------------------------------

def ligand_record(pdb_lines, ligand_id: str, name: str) -> str:
    """One SDF record for a co-crystal ligand, from its deposited HETATM lines.

    Open Babel perceives the bonds from the coordinates: a crystallographic HETATM block
    carries no connectivity. Connectivity from geometry is reliable for a 12-60 heavy-atom
    organic ligand; bond ORDER and protonation are the parts that can come out wrong, and
    that is recorded in the JSON `source`. Hydrogens are not added here -- the engine adds
    ligand hydrogens itself (dock.LIGAND_ADDH), exactly as for the production library."""
    from openbabel import pybel

    pybel.ob.obErrorLog.SetOutputLevel(0)
    text = "\n".join([*pdb_lines, "END", ""])
    mol = pybel.readstring("pdb", text)
    if not mol.atoms:
        raise ValueError(f"{ligand_id}: Open Babel read no atoms")
    mol.title = ligand_id
    mol.data["ID"] = ligand_id
    mol.data["name"] = name
    out = mol.write("sdf")
    # Open Babel stamps the molblock's program line with MMDDYYHHmm. Left in, the library's
    # sha256 changes on every minute boundary, which changes dock's unit keys and re-docks
    # the whole control. The line is replaced with a fixed one.
    lines = out.split("\n")
    if len(lines) > 1:
        lines[1] = PROGRAM_LINE
    out = "\n".join(lines)
    if dock.END not in out:
        out = out.rstrip("\n") + f"\n{dock.END}\n"
    return out


def write_control_library(labels, path: Path) -> list:
    """Write the co-crystal SDF (atomically, only when its bytes change, so a rerun does
    not invalidate dock's library hash). Returns [(ligand id, reason)] for the ligands Open
    Babel could not convert."""
    records, failed = [], []
    for label in labels:
        lines = het_lines(label.chain_pdb, label.key)
        try:
            records.append(ligand_record(lines, label.ligand_id, label.resname))
        except Exception as exc:                     # noqa: BLE001 - reported, not raised
            failed.append((label.ligand_id, f"SDF conversion failed: "
                                            f"{type(exc).__name__}: {exc}"))
    text = "".join(records)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not (path.is_file() and path.read_text(encoding="utf-8") == text):
        _atomic_text(path, text)
    return failed


def _atomic_text(path: Path, text: str) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


# --- the dock delegation ------------------------------------------------------------------

class ControlLayout(RunLayout):
    """The run folder with `docking` pointed at 08_control/docking.

    Everything `dock` writes -- receptor PDBQTs, unit markers, replicates, df_results_*,
    poses, discarded.tsv, its own manifest -- lands there, so the production results in
    07_docking are neither read nor touched by the control."""

    @property
    def docking(self) -> Path:
        return self.run_dir / DIRNAME / DOCKING_DIRNAME


def stage_dir(layout) -> Path:
    return layout.run_dir / DIRNAME


def library_path(layout) -> Path:
    return layout.run_dir / LIBRARY


def dock_args(args):
    """dock's namespace for the control run: this stage's replicate count, `--poses first`
    (replicate 0's poses are readout 2), and the engine options as given."""
    return SimpleNamespace(
        nprocs=int(getattr(args, "nprocs", 1)),
        dry_run=False,
        no_notify=bool(getattr(args, "no_notify", False)),
        library=getattr(args, "library", dock.DEFAULT_LIBRARY),
        seed=int(getattr(args, "seed", dock.DEFAULT_SEED)),
        replicates=int(getattr(args, "control_replicates", None) or DEFAULT_REPLICATES),
        threads_per_job=int(getattr(args, "threads_per_job", dock.DEFAULT_TPJ)),
        unit_timeout=int(getattr(args, "unit_timeout", None) or dock.DEFAULT_UNIT_TIMEOUT),
        engine=dock.engine_name(args),
        search_mode=getattr(args, "search_mode", None),
        batch_size=getattr(args, "batch_size", None),
        poses=dock.POSES_FIRST,
    )


def slots(args, threads_per_job: int = 1):
    """The control runs dock's units, so it parallelises exactly as dock does."""
    return dock.slots(dock_args(args), threads_per_job)


def calibration_tool(args) -> str:
    """dock's calibration entry: the unit is one ligand through the same engine."""
    return dock.calibration_tool(dock_args(args))


# --- the readouts --------------------------------------------------------------------------

def _summary_scores(layout, targets) -> tuple:
    """({target: {ligand id: score_mean}}, [(target, reason)]) from the control docking.

    Read through `stats.read_df_results`, which is the ONE df_results reader in this repo:
    it rejects non-finite and unparseable `score_mean` values and returns them. A NaN here
    would otherwise reach `self_rank_percentile`, where it compares false against
    everything and comes out BETTER than the best achievable rank (see that function)."""
    from polyphar.stages import stats

    scores, reasons = {}, []
    ctrl = ControlLayout(layout.run_dir)
    for t in sorted(targets, key=lambda t: t.label):
        path = ctrl.docking / f"df_results_{t.label}.csv"
        if not path.is_file():
            reasons.append((t.label, "no control docking summary (the target is incomplete; "
                                     f"see {DIRNAME}/{DOCKING_DIRNAME}/failed.tsv)"))
            continue
        got, dropped = stats.read_df_results(path)
        for lid, why in dropped:
            reasons.append((f"{t.label}/{lid}", f"dropped from the control readouts: {why}"))
        if not got:
            reasons.append((t.label, "every row of its control docking summary is "
                                     "unusable: the target is not scored"))
            continue
        scores[t.label] = got
    return scores, reasons


def self_ranks(labels, scores: dict) -> tuple:
    """([row], [(ligand id, reason)]) for self_rank.csv.

    A ligand is ranked against EVERY panel target that produced a summary, and only when it
    was scored on all of them: a shrinking denominator would quietly make a self-rank look
    better than it is."""
    rows, reasons = [], []
    panel = sorted(scores)
    for label in sorted(labels, key=lambda x: x.ligand_id):
        lid = label.ligand_id
        missing = [t for t in panel if lid not in scores[t]]
        if missing:
            reasons.append((lid, f"not scored on {len(missing)} of {len(panel)} panel "
                                 f"targets ({', '.join(missing[:3])}...): no self-rank"))
            continue
        if label.target not in scores:
            reasons.append((lid, f"its own target {label.target} produced no summary"))
            continue
        per_target = {t: scores[t][lid] for t in panel}
        pct = self_rank_percentile(label.target, per_target)
        rank = round(pct * len(panel) + 0.5)
        rows.append({"Ligand_ID": lid, "Ligand": label.resname,
                     "Home_Target": label.target,
                     "Home_Score": per_target[label.target],
                     "Best_Score": min(per_target.values()),
                     "Panel_Targets": len(panel), "Self_Rank": rank,
                     "Self_Rank_Percentile": pct})
    return rows, reasons


def best_poses(poses_sdf, score_field: str | None = None) -> dict:
    """{ligand id: (score, [(x, y, z, element)], reason)} for the BEST pose of each ligand
    in one engine output, chosen on the score rather than on record order. `reason` is ""
    for a readable molblock and says why not otherwise -- an unreadable pose is still the
    best pose and still counts."""
    best: dict = {}
    for rec in dock._records(dock._read_text(poses_sdf)):
        lid = rec.fields.get("ID")
        if not lid:
            title = rec.lines[0].strip()
            lid = title[len(dock.TITLE_PREFIX):] if title.startswith(dock.TITLE_PREFIX) \
                else title
        score = dock._pose_score(rec, poses_sdf, f"a pose of {lid}", score_field)
        if lid in best and best[lid][0] <= score:
            continue
        xyz, why = _molblock_heavy_atoms(rec.lines)
        best[lid] = (score, xyz, why)
    return best


def _molblock_heavy_atoms(lines) -> tuple:
    """([(x, y, z, element)], reason) for an SDF record's molblock.

    Only the V2000 connection table is parsed, which is what both engines write. A record
    this cannot read returns an empty list AND a reason, and the caller counts it AGAINST
    both readouts rather than dropping it from the denominator -- an unreadable pose is a
    pose that did not land in the pocket as far as anything downstream can tell."""
    if len(lines) < 4:
        return [], "molblock shorter than a counts line"
    counts = lines[3]
    if "V3000" in counts.upper():
        return [], ("V3000 connection table: the counts line carries no atom count and "
                    "the atom block is not parsed here")
    try:
        natoms = int(counts[0:3])
    except (ValueError, IndexError):
        return [], f"unparseable counts line {counts.strip()!r}"
    if natoms <= 0:
        return [], f"counts line declares {natoms} atoms"
    out, bad = [], 0
    for line in lines[4:4 + natoms]:
        try:
            x, y, z = float(line[0:10]), float(line[10:20]), float(line[20:30])
        except ValueError:
            bad += 1
            continue
        el = line[31:34].strip().upper()
        if el in ("H", "D"):
            continue
        out.append((x, y, z, el))
    if not out:
        return [], (f"no heavy atom among the {natoms} declared"
                    + (f" ({bad} unparseable)" if bad else ""))
    return out, ""


def pose_rows(layout, targets, engine: str, labels=()) -> tuple:
    """([row] for pose_in_pocket.csv, [(item, reason)]). One row per (target, ligand): the
    best replicate-0 pose, its centroid, whether that centroid is in the cavity
    (descriptive) and, for the SELF pairs, the redocking RMSD to the deposited ligand
    (the gate).

    An unreadable pose is NOT dropped: it gets a row with In_Pocket False and, if it is a
    self pair, Redock_Success False. Excluding it from the denominator is how a readout is
    made to look good by losing exactly the cases that would have counted against it."""
    ctrl = ControlLayout(layout.run_dir)
    field = dock.SCORE_FIELD_OF.get(engine)
    home_of = {(lab.target, lab.ligand_id): lab for lab in labels}
    deposited: dict = {}
    rows, reasons = [], []
    for t in sorted(targets, key=lambda t: t.label):
        poses = ctrl.docking / t.label / "poses.sdf"
        if not poses.is_file():
            continue          # the target is incomplete; _summary_scores already said so
        pocket = layout.run_dir / t.pocket_rel
        lining = pocket_atoms(pocket, hetatm=True)
        try:
            box_center, box_size = dock.pocket_box(pocket)
            best = best_poses(poses, field)
        except Exception as exc:                     # noqa: BLE001 - reported, not raised
            reasons.append((t.label, f"cannot read replicate-0 poses: "
                                     f"{type(exc).__name__}: {exc}"))
            continue
        for lid in sorted(best):
            score, xyz, bad = best[lid]
            lab = home_of.get((t.label, lid))
            row = {"Target": t.label, "Ligand_ID": lid, "Score": score,
                   "Lining_Atoms": "", "Asymmetry": "", "In_Pocket": False,
                   "Contact_Fraction": "", "Self_Pair": lab is not None,
                   "RMSD": "", "Redock_Success": False if lab is not None else "",
                   "Reason": bad}
            if not xyz:
                reasons.append((f"{t.label}/{lid}", f"pose not readable ({bad}): counted "
                                                    f"against both readouts"))
                rows.append(row)
                continue
            c = centroid(xyz)
            ok, n_near, asym, why = pose_in_pocket(c, lining, box_center, box_size)
            row.update({"Centroid_X": round(c[0], 3), "Centroid_Y": round(c[1], 3),
                        "Centroid_Z": round(c[2], 3), "Lining_Atoms": n_near,
                        "Asymmetry": "" if asym != asym else round(asym, 4),
                        "In_Pocket": bool(ok),
                        "Contact_Fraction": round(contact_fraction(xyz, lining), 4),
                        "Reason": why})
            if lab is not None:
                key = (str(lab.chain_pdb), lab.key)
                if key not in deposited:
                    deposited[key] = het_atoms(lab.chain_pdb, lab.key)
                try:
                    rmsd = symmetry_insensitive_rmsd(xyz, deposited[key])
                except Exception as exc:             # noqa: BLE001 - reported, not raised
                    row["Reason"] = (row["Reason"] + "; " if row["Reason"] else "") + \
                        f"no RMSD ({type(exc).__name__}: {exc}): counted as a failure"
                    reasons.append((f"{t.label}/{lid}", row["Reason"]))
                else:
                    row["RMSD"] = round(rmsd, 4)
                    row["Redock_Success"] = bool(rmsd <= REDOCK_RMSD_MAX)
            rows.append(row)
    return rows, reasons


def library_rows(layout, labels, ctrl_scores: dict) -> tuple:
    """([row] for library_percentile.csv, [(ligand id, reason)]) -- readout 3.

    The co-crystal ligand's score percentile among the PRODUCTION library at its OWN
    target: `07_docking/df_results_<home>.csv` against the control's own score for the
    same (receptor, box, engine). Lower = better, mid-ranked, so chance is 0.5."""
    from polyphar.stages import stats

    rows, reasons = [], []
    for lab in sorted(labels, key=lambda x: x.ligand_id):
        row = {"Ligand_ID": lab.ligand_id, "Ligand": lab.resname,
               "Home_Target": lab.target, "Top_Decile": False}
        home = ctrl_scores.get(lab.target, {}).get(lab.ligand_id)
        if home is None or not math.isfinite(home):
            why = (f"no usable control score at its own target {lab.target}: no library "
                   f"percentile")
            reasons.append((lab.ligand_id, why))
            rows.append({**row, "Reason": why})
            continue
        path = stats.df_results_path(layout.docking, lab.target)
        if not path.is_file():
            why = (f"the production docking has no {dock._rel(path, layout.run_dir)}: "
                   f"readout 3 compares against the production library at the same "
                   f"target and cannot be formed without it")
            reasons.append((lab.ligand_id, why))
            rows.append({**row, "Control_Score": home, "Reason": why})
            continue
        lib, _dropped = stats.read_df_results(path)
        lib.pop(lab.ligand_id, None)      # the control ligand is not part of the library
        if not lib:
            why = f"{dock._rel(path, layout.run_dir)} scored no library ligand"
            reasons.append((lab.ligand_id, why))
            rows.append({**row, "Control_Score": home, "Reason": why})
            continue
        n = len(lib) + 1
        better = sum(1 for v in lib.values() if v < home)
        ties = sum(1 for v in lib.values() if v == home) + 1
        pct = (better + (ties + 1) / 2 - 0.5) / n
        rows.append({**row, "Control_Score": home, "Library_Ligands": len(lib),
                     "Better_In_Library": better, "Percentile": pct,
                     "Top_Decile": bool(pct <= TOP_DECILE), "Reason": ""})
    return rows, reasons


# --- outputs --------------------------------------------------------------------------------

def _write_csv(path: Path, columns, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=columns)
            w.writeheader()
            for row in rows:
                w.writerow({c: row.get(c, "") for c in columns})
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def positive_control_json(rank_rows, pose_rows_, labels, n_panel, source: str,
                          lib_rows=(), n_panel_total=None) -> dict:
    """The file `stats` reads. The three gated readouts first, then everything needed to
    judge them -- including the caveat, which is part of the RESULT, not of the docs.

    `n_panel` is the number of panel targets that produced a SCORE; `n_panel_total` is how
    many there are. Both are here, not only in the manifest, because `stats` gates on the
    ratio: a target that failed is excluded from every ligand's denominator, and target
    attrition correlated with score turns a genuine failure into a pass (simulated: a
    ligand truly 30th of 482 reads 0.0589 when 100 targets fail at random and correctly
    fails the gate, but 0.0013 -- a comfortable pass -- when the 100 BEST-scoring ones
    fail, which is the plausible direction, since the most generous pockets are the
    largest, the largest have the biggest boxes, and the biggest boxes time out first)."""
    pcts = [r["Self_Rank_Percentile"] for r in rank_rows]
    n_pose = len(pose_rows_)
    in_pocket = sum(1 for r in pose_rows_ if r["In_Pocket"])
    contacts = [r["Contact_Fraction"] for r in pose_rows_
                if isinstance(r.get("Contact_Fraction"), (int, float))]
    self_pairs = [r for r in pose_rows_ if r.get("Self_Pair")]
    rmsds = [r["RMSD"] for r in self_pairs if isinstance(r.get("RMSD"), (int, float))]
    redock_ok = sum(1 for r in self_pairs if r.get("Redock_Success"))
    scored = [r for r in lib_rows if isinstance(r.get("Percentile"), (int, float))]
    lib_pcts = [r["Percentile"] for r in scored]
    top_decile = sum(1 for r in scored if r.get("Top_Decile"))
    total = int(n_panel_total if n_panel_total is not None else n_panel)
    return {
        "median_self_rank_percentile": median(pcts),
        "redock_rmsd_success_fraction": (redock_ok / len(self_pairs)) if self_pairs else 0.0,
        "library_top_decile_fraction": (top_decile / len(scored)) if scored else 0.0,
        # descriptive, NOT gated: see POSE_IN_POCKET_DEFINITION
        "pose_in_pocket_fraction": in_pocket / n_pose if n_pose else 0.0,
        "source": source,
        "n_redocked": len(rank_rows),
        "n_redock_rmsd": len(self_pairs),
        "n_library_percentile": len(scored),
        "panel_targets_total": total,
        "panel_targets_scored": int(n_panel),
        "panel_scored_fraction": (n_panel / total) if total else 0.0,
        "caveat": CAVEAT,
        "self_rank": {
            "n_ligands": len(rank_rows),
            "panel_targets": n_panel,
            "panel_targets_total": total,
            "median_percentile": median(pcts),
            "fraction_top_1pct": sum(1 for p in pcts if p <= 0.01) / len(pcts),
            "fraction_top_5pct": sum(1 for p in pcts if p <= 0.05) / len(pcts),
            "chance_percentile": 0.5,
            # what the median would be if EVERY unscored target had beaten the home pocket.
            # Not gated -- it is the worst case, not the estimate -- but it is the number
            # that says how much room score-correlated target attrition had.
            "median_percentile_adversarial": (
                (median(pcts) * n_panel + (total - n_panel)) / total) if total else 1.0,
            "definition": ("rank of the ligand's OWN pocket among every panel target it "
                           "was docked against, mid-ranked, as a percentile in [0, 1]; "
                           "0.5 is chance"),
        },
        "redock_rmsd": {
            "n_self_pairs": len(self_pairs),
            "n_within": redock_ok,
            "fraction": (redock_ok / len(self_pairs)) if self_pairs else 0.0,
            "threshold_angstrom": REDOCK_RMSD_MAX,
            "median_rmsd": median(rmsds) if rmsds else None,
            "n_unmeasurable": len(self_pairs) - len(rmsds),
            "definition": REDOCK_RMSD_DEFINITION,
        },
        "library_percentile": {
            "n_ligands": len(scored),
            "n_top_decile": top_decile,
            "fraction_top_decile": (top_decile / len(scored)) if scored else 0.0,
            "median_percentile": median(lib_pcts) if lib_pcts else None,
            "chance_percentile": 0.5,
            "top_decile": TOP_DECILE,
            "definition": LIBRARY_PERCENTILE_DEFINITION,
        },
        "pose_in_pocket": {
            "gated": False,
            "n_poses": n_pose,
            "n_in_pocket": in_pocket,
            "fraction": in_pocket / n_pose if n_pose else 0.0,
            "mean_contact_fraction": (sum(contacts) / len(contacts)) if contacts else 0.0,
            "self_pair_fraction": (sum(1 for r in self_pairs if r["In_Pocket"])
                                   / len(self_pairs)) if self_pairs else 0.0,
            "n_self_pairs": len(self_pairs),
            "definition": POSE_IN_POCKET_DEFINITION,
        },
        "label": {
            "n_labelled_pockets": len(labels),
            "criterion": (f"a HETATM residue of {MIN_HEAVY}-{MAX_HEAVY} heavy atoms, not "
                          f"solvent/ion/buffer/glycan/lipid, with >= {MIN_FRAC:.0%} of its "
                          f"heavy atoms within {CONTACT:g} A of the pocket's lining "
                          f"(rg_calibrate.py's criterion, unchanged), one row per PHYSICAL "
                          f"ligand -- a ligand filling a pocket and its nested subpocket is "
                          f"counted once"),
        },
    }


def summary_markdown(result: dict, verdicts: dict) -> str:
    sr, pp = result["self_rank"], result["pose_in_pocket"]
    rr, lp = result["redock_rmsd"], result["library_percentile"]
    return "\n".join([
        "# Positive control: does the score mean anything at all?",
        "",
        f"- **median self-rank percentile {result['median_self_rank_percentile']:.4f}** "
        f"(pass: <= {verdicts['self_rank_max']:g} over at least "
        f"{verdicts['min_self_rank_n']} ligands; chance 0.5) over "
        f"{result['n_redocked']} redocked co-crystal ligands against "
        f"{sr['panel_targets']} of the panel's {sr['panel_targets_total']} targets",
        f"- top 1%: {sr['fraction_top_1pct']:.1%} of ligands; "
        f"top 5%: {sr['fraction_top_5pct']:.1%}; worst case if every unscored target had "
        f"beaten the home pocket: {sr['median_percentile_adversarial']:.4f}",
        f"- **redocking RMSD <= {rr['threshold_angstrom']:g} A: "
        f"{result['redock_rmsd_success_fraction']:.4f}** "
        f"(pass: >= {verdicts['redock_success_min']:g} over at least "
        f"{verdicts['min_redock_n']} self pairs) -- {rr['n_within']} of "
        f"{rr['n_self_pairs']} co-crystal ligands put back where they were found"
        + (f"; median RMSD {rr['median_rmsd']:.2f} A" if rr["median_rmsd"] is not None
           else ""),
        f"- **library top decile at the home pocket: "
        f"{result['library_top_decile_fraction']:.4f}** "
        f"(pass: >= {verdicts['library_top_decile_min']:g} over at least "
        f"{verdicts['min_library_n']} ligands; chance {lp['top_decile']:g}) -- "
        f"{lp['n_top_decile']} of {lp['n_ligands']} co-crystal ligands beat 90% of the "
        f"production library at the pocket they were crystallised in",
        f"- panel coverage {result['panel_scored_fraction']:.4f} "
        f"(pass: >= {verdicts['min_panel_scored_fraction']:g})",
        f"- descriptive, NOT a gate: pose-in-pocket "
        f"{result['pose_in_pocket_fraction']:.4f} over {pp['n_poses']} replicate-0 best "
        f"poses ({pp['n_in_pocket']} inside the cavity); mean contact fraction "
        f"{pp['mean_contact_fraction']:.3f}",
        "",
        f"Verdict as `polyphar stats` will read it: **{verdicts['verdict']}**",
        "",
        "## What these numbers are",
        "",
        f"Self-rank: {sr['definition']}.",
        "",
        f"Redocking RMSD: {rr['definition']}.",
        "",
        f"Library percentile: {lp['definition']}.",
        "",
        f"Pose-in-pocket: {pp['definition']}.",
        "",
        f"Label: {result['label']['criterion']}.",
        "",
        "## Caveat",
        "",
        CAVEAT,
        "",
        f"Source: {result['source']}",
        "",
    ])


# --- fingerprint / resume -------------------------------------------------------------------

def _dock_fingerprint(layout) -> str | None:
    try:
        data = json.loads((ControlLayout(layout.run_dir).docking
                           / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data.get("fingerprint")


def fingerprint(layout, args) -> str | None:
    """What the finished outputs depend on: the panel, the co-crystal library, the docking
    that was actually done (dock's own fingerprint), and the constants that define the two
    readouts. None when the control docking has not finished, so nothing can match."""
    import hashlib

    lib = library_path(layout)
    dock_fp = _dock_fingerprint(layout)
    if dock_fp is None or not lib.is_file() or not layout.final_targets_csv.is_file():
        return None
    payload = json.dumps({
        "targets": sha256_file(layout.final_targets_csv),
        "library": sha256_file(lib),
        "dock": dock_fp,
        "label": [CONTACT, MIN_FRAC, MIN_HEAVY, MAX_HEAVY, sorted(SKIP)],
        "pose": [SHELL_RADIUS, MIN_LINING_ATOMS, MAX_ASYMMETRY, EDGE_MARGIN],
        "readouts": [REDOCK_RMSD_MAX, TOP_DECILE],
    }, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def _finished(layout, args) -> bool:
    try:
        data = json.loads((stage_dir(layout) / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    fp = fingerprint(layout, args)
    return bool(fp) and data.get("fingerprint") == fp


OUT_ROOTS = ("inputs", DIRNAME)


def out_path(value: str) -> str:
    """Validate `--positive-control-out`. It is a path this stage OVERWRITES and, when the
    control cannot be produced, UNLINKS, so it is not free-form: unvalidated, `polyphar
    control --positive-control-out 07_docking/manifest.json` deletes the docking manifest,
    `../x.json` escapes the run folder, and a first failing run removes a hand-supplied
    `inputs/positive_control.json` the user meant to keep somewhere else.

    Allowed: a relative `.json` path, no `..`, inside `inputs/` or `{DIRNAME}/`, not a
    manifest."""
    import argparse

    p = Path(str(value))
    bad = None
    if p.is_absolute():
        bad = "it must be relative to the run folder"
    elif ".." in p.parts:
        bad = "it must not contain '..'"
    elif p.suffix != ".json":
        bad = "it must name a .json file"
    elif p.name == "manifest.json":
        bad = "a stage manifest is not an output of this stage"
    elif not p.parts or p.parts[0] not in OUT_ROOTS:
        bad = f"it must be inside {' or '.join(r + '/' for r in OUT_ROOTS)}"
    if bad:
        raise argparse.ArgumentTypeError(
            f"--positive-control-out {value!r}: {bad} (this stage overwrites that path, "
            f"and removes it when no control can be produced)")
    return p.as_posix()


def add_arguments(parser) -> None:
    dock.add_engine_options(parser)
    parser.add_argument("--control-replicates", type=dock._positive_int,
                        default=DEFAULT_REPLICATES,
                        help=f"independent draws per target for the CONTROL docking "
                             f"(default {DEFAULT_REPLICATES}); the control asks whether a "
                             f"score means anything, not how much it jitters")
    parser.add_argument("--positive-control-out", type=out_path,
                        default=INPUTS_POSITIVE_CONTROL,
                        help=f"where to install the result for `polyphar stats`, relative "
                             f"to the run folder (default {INPUTS_POSITIVE_CONTROL}); must "
                             f"be a .json under inputs/ or {DIRNAME}/, because the stage "
                             f"overwrites it and removes it when no control can be "
                             f"produced. The canonical copy is always "
                             f"{DIRNAME}/{POSITIVE_CONTROL}")


def units(layout, args):
    """Ligands still to draw for the control docking (0 once the stage has finished on
    exactly these inputs). The basis is dock's -- one ligand through the engine -- so the
    ETA uses dock's measured rate (eta.CALIBRATION_ALIAS)."""
    if _finished(layout, args):
        return WorkUnits(0, "ligands")
    lib = library_path(layout)
    if lib.is_file():
        return dock.units(ControlLayout(layout.run_dir), dock_args(args), library=lib)
    # No library yet (the first run, or one whose library was removed): label in memory to
    # get the exact ligand count. It reads two files per target and writes nothing, so
    # `--dry-run` and `polyphar estimate` stay read-only. An EXISTING library is trusted as
    # it stands: `run` rewrites it before docking, so a count taken against a stale one is
    # only ever the ETA of a pass that is about to relabel anyway.
    try:
        targets = dock.read_targets(layout)
    except (FileNotFoundError, ValueError):
        return WorkUnits(0, "ligands")
    labels, _ = label_targets(layout, targets)
    reps = int(getattr(args, "control_replicates", None) or DEFAULT_REPLICATES)
    return WorkUnits(len(labels) * len(targets) * reps, "ligands")


def run(layout, args, progress) -> dict:
    d = stage_dir(layout)
    # re-validated here and not only in argparse: `run` and the tests build the namespace
    # directly, and the first thing _finish does with this path on a failing control is
    # unlink it.
    out_path(str(getattr(args, "positive_control_out", None) or INPUTS_POSITIVE_CONTROL))
    if _finished(layout, args):
        data = json.loads((d / "manifest.json").read_text(encoding="utf-8"))
        return data.get("counts") or {}
    d.mkdir(parents=True, exist_ok=True)
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    (d / "manifest.json").unlink(missing_ok=True)   # the outputs below are about to change
    rows: list[tuple[str, str]] = []

    targets = dock.read_targets(layout)
    labels, unlabelled = label_targets(layout, targets)

    lib = library_path(layout)
    counts = {"panel_targets": len(targets), "labelled": len(labels),
              "unlabelled": len(unlabelled), "failed": 0}
    try:
        if not labels:
            rows.append((DIRNAME, NO_LABELS))
            counts["failed"] += 1
            _write_labels(d, labels, unlabelled, {})
            return _finish(layout, args, d, counts, rows, None)
        bad = write_control_library(labels, lib)
        rows.extend(bad)
        counts["failed"] += len(bad)
        dropped = dict(bad)
        # labels.csv is written AFTER the conversion, so `Redocked` is what was actually
        # docked rather than what was labelled.
        _write_labels(d, labels, unlabelled, dropped)
        labels = [lab for lab in labels if lab.ligand_id not in dropped]
        counts["ligands"] = len(labels)
        if not labels:
            rows.append((DIRNAME, "every labelled ligand failed SDF conversion"))
            counts["failed"] += 1
            return _finish(layout, args, d, counts, rows, None)

        dargs = dock_args(args)
        # The library is written by now, so a host with no engine can still prepare the
        # control and hand the docking to a GPU elsewhere. Stopping here with a message
        # beats raising: the co-crystal SDF and labels.csv are the inputs
        # `dock-export --control` needs, and they are already on disk.
        try:
            dock.require_engine_host(dargs)
        except Exception as exc:
            progress.note(
                f"control: prepared {len(labels)} co-crystal ligand(s) and wrote {lib.name}, "
                f"but this host cannot run the engine: {exc} "
                f"Export the control's docking instead: `polyphar dock-export --control`, "
                f"work the bundle where the GPU is, `polyphar dock-merge --control`, then "
                f"run `polyphar control` again to compute the readouts.")
            counts["prepared_only"] = True
            return _finish(layout, args, d, counts, rows, None)
        dock_counts = dock.run(ControlLayout(layout.run_dir), dargs, progress, library=lib)
        counts.update({f"dock_{k}": v for k, v in dock_counts.items()})
        # dock's failures are the CONTROL's failures: folded in only as `dock_*` they are
        # invisible to cli.flagged(), so a control whose receptors would not prepare exits
        # 0 and `run`'s summary never flags it, while `dock` exits 3 for the same thing.
        counts["failed"] += int(dock_counts.get("failed", 0) or 0)
        counts["incomplete_targets"] = int(dock_counts.get("incomplete_targets", 0) or 0)

        scores, missing_targets = _summary_scores(layout, targets)
        rows.extend(missing_targets)
        counts["incomplete_targets"] = max(counts["incomplete_targets"],
                                           len(targets) - len(scores))
        rank_rows, rank_reasons = self_ranks(labels, scores)
        rows.extend(rank_reasons)
        _write_csv(d / "self_rank.csv", SELF_RANK_COLUMNS, rank_rows)
        poses, pose_reasons = pose_rows(layout, targets, dargs.engine, labels)
        rows.extend(pose_reasons)
        _write_csv(d / "pose_in_pocket.csv", POSE_COLUMNS, poses)
        lib_rows, lib_reasons = library_rows(layout, labels, scores)
        rows.extend(lib_reasons)
        _write_csv(d / "library_percentile.csv", LIBRARY_COLUMNS, lib_rows)
        counts.update({"self_ranked": len(rank_rows), "poses": len(poses),
                       "panel_targets_scored": len(scores),
                       "self_pairs": sum(1 for r in poses if r.get("Self_Pair")),
                       "library_percentiles": sum(
                           1 for r in lib_rows
                           if isinstance(r.get("Percentile"), (int, float)))})
        if not rank_rows or not poses:
            rows.append((DIRNAME, NO_READOUT.format(ranks=len(rank_rows), poses=len(poses))))
            counts["failed"] += 1
            return _finish(layout, args, d, counts, rows, None, wrote_readouts=True)

        result = positive_control_json(
            rank_rows, poses, labels, len(scores),
            _source(layout, args, dargs, len(rank_rows), len(scores), len(targets)),
            lib_rows=lib_rows, n_panel_total=len(targets))
        counts["median_self_rank_percentile"] = result["median_self_rank_percentile"]
        counts["redock_rmsd_success_fraction"] = result["redock_rmsd_success_fraction"]
        counts["library_top_decile_fraction"] = result["library_top_decile_fraction"]
        counts["pose_in_pocket_fraction"] = result["pose_in_pocket_fraction"]
        return _finish(layout, args, d, counts, rows, result)
    finally:
        for item, reason in rows:
            append_failure(d, item, STAGE, reason)


NO_LABELS = ("no panel target carries a co-crystal ligand, so no positive control can be "
             "produced; inputs/positive_control.json is NOT written and `polyphar stats` "
             "will report the primary result as uninterpretable, which is the correct "
             "outcome -- a control that was not run is not a control that passed")
NO_READOUT = ("the control docking produced {ranks} self-rank(s) and {poses} pose(s): "
              "not enough for either readout, so no positive_control.json is written")


def _write_labels(d: Path, labels, unlabelled, dropped: dict) -> None:
    """labels.csv: one row per panel target, labelled first, then the ones with no label
    and the reason -- which is itself the useful half, because "no co-crystal ligand" is
    not "not druggable" (see CAVEAT)."""
    _write_csv(d / "labels.csv", LABEL_COLUMNS,
               [{"Target": lab.target, "UniProt_ID": lab.accession,
                 "Pocket": lab.pocket_rel, "Ligand_ID": lab.ligand_id,
                 "Ligand": lab.resname, "Chain": lab.chain, "Residue": lab.residue,
                 "Heavy_Atoms": lab.heavy_atoms, "Contact_Fraction": round(lab.contact, 4),
                 "Redocked": lab.ligand_id not in dropped,
                 "Reason": dropped.get(lab.ligand_id, "")} for lab in labels]
               + [{"Target": t, "Redocked": False, "Reason": why}
                  for t, why in unlabelled])


def _source(layout, args, dargs, n_ligands: int, n_scored: int, n_total: int) -> str:
    """How the numbers were made. `n_scored` of `n_total`, never "all {n_scored}": with 100
    of 482 targets failed the old wording read "against all 382 panel targets", and that
    sentence gets quoted."""
    engine = dock.engine_version(ControlLayout(layout.run_dir), dargs)
    whole = "all " if n_scored == n_total else ""
    return (f"polyphar {STAGE}: {n_ligands} co-crystal ligand(s) extracted from the panel's "
            f"own experimental chains (Open Babel bond perception from the deposited "
            f"coordinates; bond order and protonation may differ from the depositor's "
            f"chemistry) and redocked against {whole}{n_scored} of the panel's {n_total} "
            f"targets with engine {dargs.engine} (version {engine}), "
            f"{dargs.replicates} replicate(s), seed {dargs.seed}; redocking RMSD and "
            f"pose-in-pocket from replicate 0. Full provenance: "
            f"{DIRNAME}/manifest.json and {DIRNAME}/{DOCKING_DIRNAME}/manifest.json.")


def _finish(layout, args, d: Path, counts: dict, rows, result,
            wrote_readouts: bool = False) -> dict:
    """Write the result, the summary, the installed copy and the manifest (last).

    `wrote_readouts` says whether self_rank.csv / pose_in_pocket.csv belong to THIS run;
    when they do not and there is no result, they are an earlier run's and are removed."""
    from polyphar.stages import stats

    out_rel = str(getattr(args, "positive_control_out", None) or INPUTS_POSITIVE_CONTROL)
    installed = layout.run_dir / out_rel
    if result is None:
        # A control that cannot be produced must not leave a stale PASS behind -- not for
        # `stats`, which reads the installed JSON, and not for a reader, who reads
        # summary.md. An earlier run's readouts go with them.
        stale = [POSITIVE_CONTROL, "summary.md"]
        if not wrote_readouts:      # an earlier run's readouts, not this one's
            stale += ["self_rank.csv", "pose_in_pocket.csv", "library_percentile.csv"]
        for name in stale:
            (d / name).unlink(missing_ok=True)
        installed.unlink(missing_ok=True)
        counts["installed"] = 0
    else:
        text = json.dumps(result, indent=2, sort_keys=False) + "\n"
        _atomic_text(d / POSITIVE_CONTROL, text)
        installed.parent.mkdir(parents=True, exist_ok=True)
        _atomic_text(installed, text)
        verdict, _detail = stats.positive_control_verdict(installed, display=out_rel)
        counts["installed"] = 1
        counts["verdict"] = verdict.split(":")[0]
        _atomic_text(d / "summary.md",
                     summary_markdown(result, {"verdict": verdict, **stats.THRESHOLDS}))
        if verdict != "yes":
            # NOT counts["failed"]: a control that ran and did not pass is a RESULT, not a
            # broken stage. It is recorded here and in the JSON, and `stats` is what turns
            # it into UNINTERPRETABLE. Only a control that could not be produced flags.
            rows.append((out_rel, f"the positive control does NOT pass -- {verdict}"))
    write_manifest(d, STAGE,
                   args={"nprocs": args.nprocs,
                         "control_replicates": int(getattr(args, "control_replicates", None)
                                                   or DEFAULT_REPLICATES),
                         "engine": dock.engine_name(args),
                         "positive_control_out": out_rel},
                   inputs=[p for p in (layout.final_targets_csv, library_path(layout))
                           if p.is_file()],
                   versions={"openbabel": dock._openbabel_version()},
                   counts=counts, timing={},
                   extra={"fingerprint": fingerprint(layout, args),
                          "caveat": CAVEAT,
                          "thresholds": dict(stats.THRESHOLDS),
                          "outputs": sorted(p.name for p in d.glob("*")
                                            if p.is_file() and p.name != "manifest.json")})
    return counts
