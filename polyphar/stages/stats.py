"""Stage 08: the cumulative weak-binding test.

Implements `docs/superpowers/specs/2026-11-03-statistics-design.md` (items 3-5 of its
implementation order) and is pre-registered by
`docs/superpowers/specs/2026-11-04-preregistration.md`. Read the spec before changing
anything here: every step below exists to make ONE guarantee hold -- that **ligand size
and pocket generosity cannot produce the result** -- and simplifying a step away usually
removes exactly that guarantee.

The claim under test
--------------------
    Approved HCC drugs bind more of the HCC panel, at a useful level, than
    comparable molecules do.

Step 1 -- within-target percentile
    For every target, every ligand scored on it is ranked by `score_mean` and converted
    to a percentile in [0, 1]. More negative = better binding = percentile near 1. Ties
    take the mid-rank, and percentile = 1 - (rank - 0.5)/N, so the percentiles of one
    target have mean exactly 0.5 whatever its score distribution is.

    This REMOVES pocket generosity rather than adjusting for it: a generous cavity scores
    everything well, but its percentiles are uniform like every other target's. No
    descriptors, no fitted model, no calibration run.

Step 2 -- the statistic

        M_i = mean over PROTEINS of ( max over that protein's targets of percentile_i )

    Per protein, best pocket. NOT a mean over all 482 targets: that would weight a protein
    contributing two pockets twice. The panel is 271 proteins / 482 targets, and the claim
    is "how much of the proteome does this molecule touch".

Step 3 -- size matching
    Strata = quintile bins of heavy-atom count x tertile bins of rotatable-bond count.
    The bin edges are computed from the LIBRARY once and written to `strata.json` BEFORE
    any docking score is read (see `run`), and are compared against the pre-registered
    values in `PREREGISTERED_EDGES`.

    The rotatable-bond count is the `ligands` stage's pinned, non-standard definition
    (single, non-ring by bridge detection, non-terminal, amide C-N NOT excluded). It is
    read from `07_docking/library/ligands.csv` and never recomputed with RDKit, whose
    conventional Lipinski count differs for 1,316 of the 2,083 e-Drug3D records.

    **The control set is already "other approved drugs of the same size".** The library is
    e-Drug3D -- approved and withdrawn small-molecule drugs, not a generic screening set --
    so size-matched non-positives are themselves approved drugs. This is load-bearing and
    easy to miss, so it is written into every output: `CONTROL_SET` appears as a column of
    results.csv, in summary.md and in the manifest.

Step 4 -- the test
    Statistic = (mean M over the positive set) - (mean M over the matched controls),
    computed stratum by stratum and combined with the positives' own stratum weights
    (`stratified_difference`): the controls' mean is therefore taken at the positives'
    size composition, which is what "matched" means here.
    Null: permute the drug/non-drug label WITHIN strata, 10,000 times, seeded.
    One-sided (the claim is directional: positives > controls).
    p = (1 + #{permutations >= observed}) / (1 + B) -- the add-one convention, so p is
    never reported as 0 (see `permutation_test`).

Step 5 -- replicate spread
    `M` is recomputed from each replicate's own raw poses (07_docking/<label>/replicates/
    r<i>.sdf) and the arm statistic reported per replicate alongside the permutation p.
    No analytic SE is derived for `M`, and in particular nothing here computes
    mean|delta| / sqrt(N) -- for a normal E|X - Y| = 1.128 sigma, so that figure is biased
    by ~13% and is not the sd of anything (see dock.summarize).

Step 4b -- the pre-registered prediction's OWN test
    HCC-on-HCC > HCC-on-AD is tested as a **paired within-stratum permutation on
    M_primary - M_control** (`paired_difference`, the `contrast` arm). Each primary arm is
    tested against its own matched controls, so reading the prediction off "one arm is
    significant and the other is not" is the difference-of-significance error. The same
    ligand carries an M on both panels, so the difference is paired within ligand.

Step 5b -- what stratification does NOT remove
    Strata remove the BETWEEN-stratum size confound. They do not make the positives and
    their controls the same size INSIDE a stratum, and a within-stratum heavy-atom
    imbalance of ~0.6 sd can produce p < 0.05 from size alone (see `balance_rows`). The
    test is not changed for that -- that would be tuning the analysis to the objection.
    Two things are added: `balance.csv`, the per-stratum size balance of every arm, and a
    pre-registered SENSITIVITY arm on heavy-atom deciles (`HA_QUANTILES_FINE`).

Step 6 -- the arms
    primary      HCC drugs on this panel; HCC drugs on the control panel (--control-panel)
    contrast     the pre-registered prediction itself (step 4b), uncorrected
    subgroup     targeted / cytotoxic split of the HCC set, pre-planned, Holm within family
    sensitivity  the primary HCC arm on finer heavy-atom bins (step 5b), uncorrected
    exploratory  AD drugs on either panel (n about 5; n is on every row)
    Holm corrects the two primary tests; the subgroup family is Holm-corrected within
    itself; contrast, sensitivity and exploratory arms carry no correction. Note that Holm
    makes the NEGATIVE CONTROL harder to reject, which is backwards for an arm one wants
    to fail: both its raw and adjusted p are reported and the raw one is named as the
    honest figure. Per-drug results are exploratory with Benjamini-Hochberg FDR over the
    resolvable subset and are labelled as such in per_drug.csv.

Step 7 -- the descriptive surface
    `surface.csv` / `surface.md` give the score_threshold x min_targets grid with NO
    p-values attached. It is how the field states these results and it should be shown; it
    is not the test. `write_no_pvalues` is an ALLOWLIST -- only the fields in
    SURFACE_COLUMNS may be written, in the column list or smuggled through a row dict --
    and `surface_markdown` is guarded by the same check, so neither file can carry one.
    PVALUE_LIKE remains as a second line applied to the allowlist itself.

The positive-control gate
    The spec pre-registers that if the positive controls fail, the primary result is
    reported as UNINTERPRETABLE rather than reinterpreted. That is enforced here, not left
    to prose: `--positive-control` (default inputs/positive_control.json) is read and every
    primary/subgroup row gets an `interpretable` verdict -- "yes", "no: ..." or
    "unknown: positive control not run" / "unknown: ... does not support a verdict". The
    redocking itself (self-rank, redocking RMSD, library percentile at home) is a separate
    piece of work -- `polyphar control` -- and this stage only consumes and enforces its
    verdict. A readout over fewer observations than the pre-registered minimum, or a panel
    too thinly scored to trust the self-rank denominator, is `unknown`: not a pass, and
    not a failure either.

PocketMatch note
    Nothing in this stage ranks on, reads or reports PocketMatch similarity. Any ranking
    added here later must use **Pmax, not Pmin**: `Pmin == 1` is near-meaningless
    (containment dominates it) and ~49,600 artificial 1.0 rows sit at the top of any
    distribution built on it. See the spec's `parse_scores` section.

Inputs
    06_targets/final_targets.csv            label (Protein_Pocket) -> UniProt_ID
    07_docking/df_results_<label>.csv       Ligand, ID, Pocket, n, score_mean, ...
    07_docking/<label>/replicates/r<i>.sdf  raw per-replicate poses (step 5)
    07_docking/library/ligands.csv          heavy_atoms / rotatable_bonds (pinned)
    inputs/positive_sets.csv                THE POSITIVE SET -- see read_positive_sets
    inputs/positive_control.json            optional; the gate above

Outputs in 08_stats/
    strata.json      bin edges + bin sizes, written before any outcome is read
    m_statistic_<panel>.csv   per ligand: M, M per replicate, size, stratum, membership
    results.csv      one row per arm (the test)
    balance.csv      per arm and stratum: the positives' size against their controls'
    per_drug.csv     exploratory per-drug elevation, BH FDR
    surface.csv      the descriptive threshold x min_targets grid (no p-values)
    surface.md       the same grid as a table
    summary.md       what a reader needs to not misread the above
    failed.tsv, manifest.json (written last = done)

Determinism: one seeded numpy Generator does every permutation in this process, the pool
is used only to read files into dicts, and every table is written in a sorted order, so
`--nprocs` cannot change a byte of the output.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path

from polyphar import pool
from polyphar.common import WorkUnits, append_failure, sha256_file, write_manifest

STAGE = "stats"
DIRNAME = "08_stats"

# --- pre-registered parameters (changing one invalidates the pre-registration) ---------
N_PERMUTATIONS = 10000
SEED = 42
HA_QUANTILES = (0.2, 0.4, 0.6, 0.8)        # quintiles of heavy-atom count
RB_QUANTILES = (1 / 3, 2 / 3)              # tertiles of rotatable-bond count
# The pre-registered SENSITIVITY stratification (step 3b): heavy-atom deciles crossed with
# the same rotatable-bond tertiles. Not the primary analysis. It exists because the primary
# quintile bins are wide -- HA bin 4 is open-ended and spans 32-63 heavy atoms -- and a
# within-stratum size imbalance of about 0.6 sd can produce p < 0.05 from size alone.
HA_QUANTILES_FINE = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# A within-stratum standardised size difference at or above this is reported as a size
# imbalance that could by itself produce the result (see balance_rows / P2 in summary.md).
BALANCE_FLAG_SD = 0.6
# Computed from e-Drug3D_2083's 1,960 dockable records on 2026-11-04, before any docking
# score existed, and recorded in docs/superpowers/specs/2026-11-04-preregistration.md.
# The stage recomputes them from the run's own library and reports a mismatch loudly
# (counts["strata_edges_preregistered"]), it does not silently adopt new edges.
PREREGISTERED_EDGES = {"heavy_atoms": [17.0, 22.0, 26.0, 32.0],
                       "rotatable_bonds": [4.0, 7.0]}
PREREGISTERED_LIBRARY = "e-Drug3D_2083 (1,960 records surviving the dockability filter)"
# The pre-registration document itself. Its sha256 is pinned here so that editing it after
# a number is known fails the test suite until the constant and an `## Amendments` entry
# are updated together (tests/test_stage_stats.py::
# test_the_preregistration_is_the_one_this_code_was_written_against).
PREREGISTRATION = "docs/superpowers/specs/2026-11-04-preregistration.md"
PREREGISTRATION_SHA256 = "4dff92a8c81d7ff8372c3197c30a3d1615c28251c96e58f012e0ba31a6f41747"
# A ligand must be scored on at least this share of the panel's proteins AND of the panel's
# targets to carry an M. 0.98, not 0.9: under informative missingness -- a ligand failing
# where it fits worst, which is the realistic direction -- M is biased by up to
# 0.5 * (1 - coverage), i.e. +0.051 at 0.90 and +0.011 at 0.98, and the effect this test
# calls significant is of that order.
MIN_COVERAGE = 0.98
# Positive-control pass thresholds, pre-registered (§ 6). See the module docstring.
SELF_RANK_MAX = 0.05        # median self-recovery rank within the top 5% of the panel
# Readout 2 is the REDOCKING RMSD, not pose-in-pocket. `dock.pocket_box` is the unpadded
# bounding box of the pocket's own lining atoms and the engine keeps the ligand inside it,
# so "centroid inside the box, away from the faces" nearly restates a constraint the engine
# already enforces: measured on real net2 pockets, random NON-CLASHING centroids pass
# pose-in-pocket ~79% of the time against ~85% for the true deposited poses, i.e. the old
# `POSE_IN_POCKET_MIN = 0.70` sat BELOW its own null and an engine placing ligands at
# random among non-clashing positions in its own box passed it. RMSD <= 2 A is the
# field-standard redocking-success criterion, its null on the same random placements is
# 0.0%, and it is the number a reviewer asks for.
REDOCK_RMSD_MAX = 2.0         # A; what counts as recovering the deposited pose
REDOCK_SUCCESS_MIN = 0.50     # >= 50% of self pairs within it
# Readout 3: the axis M is actually built on. M is within-target across ligands; the
# self-rank is within-ligand across targets -- rows and columns of the same matrix, and one
# does not imply the other. A known binder that does not reach the top decile of the
# production library at the pocket it was crystallised in means a within-target percentile
# carries no site-discriminative content. Chance is 0.10.
LIBRARY_TOP_DECILE_MIN = 0.50
# Minimum n per readout. Without one, ONE labelled ligand ranking first on a 20-target
# panel reads `median 0.025, pose 1.0` and the verdict is `yes`. The realistic n is about
# 24 (482 targets, ~37% experimental, ~13.5% of panel-style experimental pockets labelled),
# not the 687 the pre-registration first quoted, so the minimum has to be low enough to be
# reachable and high enough that one ligand cannot carry it: at n = 10 the binomial
# standard error on a fraction is 0.16 and the median is not set by a single value.
MIN_SELF_RANK_N = 10
MIN_REDOCK_N = 10
MIN_LIBRARY_N = 10
# A target that failed is excluded from EVERY ligand's self-rank denominator, so target
# attrition correlated with score turns a genuine failure into a pass. The same 0.98 as the
# ligand coverage rule (§ 2), and for the same reason: at 0.98 the worst case -- every
# unscored target having beaten the home pocket -- moves the median by at most
# (1 - 0.98) = 0.0207 on a 482-target panel, against 0.21 at 0.79.
MIN_PANEL_SCORED_FRACTION = 0.98

# What the gate prints, and what every manifest records. One dict so the thresholds cannot
# drift apart between the verdict, the summary and the provenance.
THRESHOLDS = {
    "self_rank_max": SELF_RANK_MAX,
    "redock_rmsd_max_angstrom": REDOCK_RMSD_MAX,
    "redock_success_min": REDOCK_SUCCESS_MIN,
    "library_top_decile_min": LIBRARY_TOP_DECILE_MIN,
    "min_self_rank_n": MIN_SELF_RANK_N,
    "min_redock_n": MIN_REDOCK_N,
    "min_library_n": MIN_LIBRARY_N,
    "min_panel_scored_fraction": MIN_PANEL_SCORED_FRACTION,
}

SET_HCC = "hcc"
SET_AD = "ad"
SUBGROUPS = ("targeted", "cytotoxic")
FAMILY_PRIMARY = "primary"
FAMILY_SUBGROUP = "subgroup"
FAMILY_EXPLORATORY = "exploratory"
# The pre-registered prediction's OWN test (step 4b). It is not in the Holm family of two:
# it is one test of one pre-registered claim, reported with its raw p.
FAMILY_CONTRAST = "contrast"
# The pre-registered sensitivity analysis (step 3b), finer heavy-atom bins. Not primary.
FAMILY_SENSITIVITY = "sensitivity"
DIRECTION = "positives > controls (one-sided)"
DIRECTION_CONTRAST = ("M on the primary panel - M on the control panel, positives > "
                      "controls (one-sided, paired within ligand)")

CONTROL_SET = ("size-matched e-Drug3D molecules that are not in the positive set; "
               "e-Drug3D is a library of APPROVED and WITHDRAWN small-molecule drugs, so "
               "the controls are other approved drugs of the same size, not a generic "
               "screening set")

# The descriptive surface's grid (score_mean in kcal/mol, more negative = better).
SURFACE_THRESHOLDS = (-11.0, -10.5, -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0)
SURFACE_MIN_TARGETS = (1, 2, 3, 5, 10, 20, 50)

POSITIVE_SETS = "inputs/positive_sets.csv"
POSITIVE_CONTROL = "inputs/positive_control.json"
POSITIVE_SET_COLUMNS = ("set", "ID", "name", "subgroup", "source")

STRATA = "strata.json"
RESULTS = "results.csv"
PER_DRUG = "per_drug.csv"
SURFACE = "surface.csv"
SURFACE_MD = "surface.md"
SUMMARY_MD = "summary.md"
BALANCE = "balance.csv"

RESULT_COLUMNS = [
    "arm", "drug_set", "panel", "family", "n_positive", "n_control", "n_strata",
    "statistic", "p_value", "p_adjusted", "adjustment", "direction",
    "stat_r0", "stat_r1", "stat_r2", "replicate_spread",
    "replicate_sign_jitter_check",
    "balance_sd", "balance_sd_max_stratum", "balance_flag",
    "interpretable", "control_set", "note",
]
BALANCE_COLUMNS = ["arm", "panel", "drug_set", "family", "stratum", "weight",
                   "n_positive", "n_control",
                   "ha_mean_positive", "ha_mean_control", "ha_median_positive",
                   "ha_median_control", "ha_sd_control", "ha_standardised_difference",
                   "targets_mean_positive", "targets_mean_control", "flag"]
PER_DRUG_COLUMNS = ["panel", "drug_set", "ID", "name", "heavy_atoms", "rotatable_bonds",
                    "stratum", "M", "n_control_in_stratum", "p_value", "q_value_bh",
                    "status"]
SURFACE_COLUMNS = ["panel", "drug_set", "score_threshold", "min_targets",
                   "n_positive", "n_positive_meeting", "frac_positive_meeting",
                   "n_control", "n_control_meeting", "frac_control_meeting"]
# one `M_r<i>` column per replicate found on disk is appended to these (see run)
M_COLUMNS = ["ID", "name", "heavy_atoms", "rotatable_bonds", "stratum", "ha_bin", "rb_bin",
             "stratum_fine",
             "n_proteins", "n_targets", "coverage", "target_coverage", "usable", "M",
             "sets", "subgroups"]

POSITIVE_SET_MISSING = (
    "the positive set {path} is not there. `polyphar stats` cannot invent it: it is the "
    "clinical-guideline (EASL / NCCN / ESMO) systemic therapy list, intersected with the "
    "dockable e-Drug3D records, and building it is guideline reading, not code. Write a "
    "CSV with the columns " + ",".join(POSITIVE_SET_COLUMNS) + " (one row per drug; `set` "
    "is '" + SET_HCC + "' or '" + SET_AD + "'; give `ID` -- the e-Drug3D library ID -- or "
    "`name`, which is matched case-insensitively against the library; `subgroup` is "
    "'targeted' or 'cytotoxic' for the HCC set; `source` is the guideline it came from), "
    "or point --positive-set at one. tests/fixtures/positive_sets.csv is a worked example.")

# Anything that looks like a p-value, a q-value or a significance verdict. This is the
# SECOND line of defence only: `write_no_pvalues` is an ALLOWLIST (a column not in the
# expected set is refused outright), because a denylist on free text always has holes --
# `p-value`, `adj.P.Val`, `Pr(>|t|)`, `pFDR` and a Cyrillic homoglyph all defeated the
# earlier `[_\s]`-bounded version of this pattern.
PVALUE_LIKE = re.compile(
    r"(?i)(^|[^a-z0-9])("
    r"p|pr|pval\w*|pvalue|p_?value|padj|p_?adj\w*|q|qval\w*|q_?value|"
    r"fdr\w*|\w*fdr|bh|holm\w*|bonferroni|benjamini|hochberg|alpha|"
    r"signif\w*|significance|sig|star\w*|asterisk\w*|"
    r"prob\w*|confidence|ci\d*|ci_\w+|nominal|"
    r"\w*log\d*p|p\d+|two_?sided|one_?sided|exact_?p|perm_?p|p_?perm|evidence|verdict"
    r")([^a-z0-9]|$)")

# Characters that make a homoglyph attack on a column name possible. Column names are
# normalised through this before the denylist sees them; the allowlist does not care,
# but the denylist must not be fooled by a zero-width space or a Cyrillic `р`.
_HOMOGLYPHS = str.maketrans({
    "​": "", "‌": "", "‍": "", "﻿": "", "­": "",
    "р": "p", "ρ": "p", "с": "c", "а": "a", "е": "e",
    "о": "o", "х": "x", "у": "y", "і": "i", "ѕ": "s",
    "А": "A", "В": "B", "Е": "E", "К": "K", "М": "M",
    "Н": "H", "О": "O", "Р": "P", "С": "C", "Т": "T",
    "Х": "X", "Ѕ": "S", "І": "I",
})


def normalise_column(name) -> str:
    """A column name with compatibility forms, zero-width characters and the common
    Cyrillic/Greek homoglyphs folded away, so the denylist cannot be spelled around."""
    import unicodedata

    return unicodedata.normalize("NFKC", str(name)).translate(_HOMOGLYPHS)


def pvalue_like(name) -> bool:
    return bool(PVALUE_LIKE.search(normalise_column(name)))


def add_arguments(parser) -> None:
    parser.add_argument("--positive-set", default=POSITIVE_SETS,
                        help=f"positive-set CSV, relative to the run folder "
                             f"(default {POSITIVE_SETS})")
    parser.add_argument("--control-panel", default=None, metavar="RUNDIR",
                        help="another finished run folder (a different disease panel) to "
                             "run the SAME drug sets against: the pre-registered negative "
                             "control. Without it the second primary test is not run and "
                             "every output says so")
    parser.add_argument("--positive-control", default=POSITIVE_CONTROL,
                        help=f"positive-control verdict JSON (default {POSITIVE_CONTROL}); "
                             f"absent means the primary results are reported as "
                             f"uninterpretable-unknown, never as interpretable")
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS, metavar="B",
                        help=f"label permutations within strata (default {N_PERMUTATIONS})")
    # NOT `--seed`: `dock` owns that one (the docking base seed) and `polyphar run`
    # merges the stages' options, so sharing it would make one flag move two things.
    parser.add_argument("--permutation-seed", type=int, default=SEED, metavar="N",
                        help=f"seed for the label permutations (default {SEED})")
    parser.add_argument("--min-coverage", type=float, default=MIN_COVERAGE, metavar="F",
                        help=f"share of the panel's proteins AND of its targets a ligand "
                             f"must be scored on to carry an M (default {MIN_COVERAGE})")


def preregistration_path() -> Path:
    """The pre-registration document, beside this repo's `polyphar/` package."""
    return Path(__file__).resolve().parents[2] / PREREGISTRATION


def preregistration_sha256() -> str:
    """The document's sha256 as it is on disk now, or "" if it is not there (the installed
    image carries the package, not the docs tree).

    Recorded in every manifest next to PREREGISTRATION_SHA256, so a run's outputs say which
    text of the pre-registration produced them. An edit without an `## Amendments` entry
    and a matching constant fails the test suite."""
    p = preregistration_path()
    return sha256_file(p) if p.is_file() else ""


# --- step 1: within-target percentiles -------------------------------------------------

def midranks(values) -> list:
    """Ascending ranks 1..N with ties averaged (no scipy)."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        r = (i + j) / 2.0 + 1.0          # mean of the 1-based ranks i+1 .. j+1
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    return ranks


def percentiles(values) -> list:
    """Within-target percentiles in [0, 1]: 1 = best = MOST NEGATIVE score.

    Mid-ranks for ties and percentile = 1 - (rank - 0.5)/N, so the mean is exactly 0.5 on
    every target whatever its score distribution looks like -- which is the whole point:
    pocket generosity is removed by construction, not adjusted for."""
    n = len(values)
    if n == 0:
        return []
    if n == 1:
        return [0.5]
    return [1.0 - (r - 0.5) / n for r in midranks(values)]


# --- step 2: M -------------------------------------------------------------------------

def m_statistic(per_target: dict, protein_of: dict) -> dict:
    """{ligand ID: (M, n_proteins, n_targets)} from {label: {ID: percentile}}.

    M = mean over PROTEINS of the max percentile over that protein's targets. A protein on
    which the ligand was never scored does not contribute; `n_proteins` says how many did.
    Weighting by target instead would count a protein contributing two pockets twice."""
    best: dict = {}      # ID -> {protein: max percentile}
    hits: dict = {}      # ID -> number of targets
    for label, scores in per_target.items():
        protein = protein_of.get(label)
        if protein is None:
            continue
        for lid, pct in scores.items():
            d = best.setdefault(lid, {})
            if pct > d.get(protein, -1.0):
                d[protein] = pct
            hits[lid] = hits.get(lid, 0) + 1
    out = {}
    for lid, by_protein in best.items():
        vals = list(by_protein.values())
        out[lid] = (sum(vals) / len(vals), len(vals), hits[lid])
    return out


def percentile_matrix(scores_by_label: dict) -> dict:
    """{label: {ID: percentile}} from {label: {ID: score}} (step 1 over every target)."""
    out = {}
    for label, scores in scores_by_label.items():
        ids = sorted(scores)
        out[label] = dict(zip(ids, percentiles([scores[i] for i in ids])))
    return out


# --- step 3: strata --------------------------------------------------------------------

def quantile_edges(values, quantiles) -> list:
    """Bin edges at `quantiles` (numpy's linear interpolation), rounded to 6 decimals.

    Heavy-atom and rotatable-bond counts are integers, so the bins come out unequal; that
    is expected and the realised bin sizes are recorded in strata.json beside the edges."""
    import numpy as np

    if not len(values):
        return [0.0 for _ in quantiles]
    return [round(float(x), 6) for x in np.quantile(np.asarray(values, dtype=float),
                                                    list(quantiles))]


def bin_of(value, edges) -> int:
    """Bin index: the number of edges <= value, i.e. the half-open interval
    [edges[k-1], edges[k])."""
    return sum(1 for e in edges if e <= value)


def stratum_label(ha_bin: int, rb_bin: int) -> str:
    return f"HA{ha_bin}-RB{rb_bin}"


# --- step 4: the test ------------------------------------------------------------------

def stratified_difference(m, positive, stratum):
    """(statistic, per-stratum detail) for the stratified mean difference.

    stat = sum_s w_s * (mean M of positives in s - mean M of controls in s), with
    w_s = n_positive_s / n_positive. The controls are therefore averaged at the POSITIVES'
    size composition -- which is what "matched controls" means -- and the weights are fixed
    under a within-stratum label permutation, so the permutation null distributes exactly
    this statistic. Strata holding no positive, or no control, contribute nothing and are
    excluded from both the observed statistic and the null."""
    by: dict = {}
    for value, pos, s in zip(m, positive, stratum):
        p, c = by.setdefault(s, ([], []))
        (p if pos else c).append(float(value))
    usable = {s: (p, c) for s, (p, c) in by.items() if p and c}
    n_pos = sum(len(p) for p, _c in usable.values())
    if not n_pos:
        return float("nan"), {}
    stat, detail = 0.0, {}
    for s in sorted(usable):
        p, c = usable[s]
        w = len(p) / n_pos
        d = sum(p) / len(p) - sum(c) / len(c)
        stat += w * d
        detail[s] = {"n_positive": len(p), "n_control": len(c), "weight": w,
                     "mean_positive": sum(p) / len(p), "mean_control": sum(c) / len(c),
                     "difference": d}
    return stat, detail


def permutation_test(m, positive, stratum, n_perm: int = N_PERMUTATIONS,
                     seed: int = SEED) -> dict:
    """One-sided within-stratum label permutation test.

    Returns {statistic, p_value, n_permutations, n_positive, n_control, n_strata, strata}.
    p = (1 + #{permuted >= observed}) / (1 + B): the add-one convention, which keeps the
    test valid at the reported level and never returns p = 0 (the design document says
    "proportion of permutations at least as extreme"; with 10,000 draws the difference is
    1e-4 and only ever makes p larger)."""
    import numpy as np

    stat, detail = stratified_difference(m, positive, stratum)
    n_pos = sum(d["n_positive"] for d in detail.values())
    n_ctl = sum(d["n_control"] for d in detail.values())
    if not detail or not math.isfinite(stat):
        return {"statistic": float("nan"), "p_value": float("nan"),
                "n_permutations": 0, "n_positive": n_pos, "n_control": n_ctl,
                "n_strata": 0, "strata": {}}
    rng = np.random.default_rng(seed)
    draws = np.zeros(n_perm, dtype=float)
    values_by_stratum = {}
    for value, pos, s in zip(m, positive, stratum):
        if s in detail:
            values_by_stratum.setdefault(s, []).append(float(value))
    for s in sorted(detail):          # sorted: the rng stream cannot depend on dict order
        d = detail[s]
        vals = np.asarray(values_by_stratum[s], dtype=float)
        k, total = d["n_positive"], float(vals.sum())
        shuffled = rng.permuted(np.tile(vals, (n_perm, 1)), axis=1)
        pos_mean = shuffled[:, :k].mean(axis=1)
        ctl_mean = (total - pos_mean * k) / d["n_control"]
        draws += d["weight"] * (pos_mean - ctl_mean)
    p = (1 + int(np.count_nonzero(draws >= stat))) / (n_perm + 1)
    return {"statistic": stat, "p_value": p, "n_permutations": n_perm,
            "n_positive": n_pos, "n_control": n_ctl, "n_strata": len(detail),
            "strata": detail}


def paired_difference(m_primary, m_control) -> dict:
    """{ID: M on the primary panel - M on the control panel}, over ligands carrying an M
    on BOTH panels.

    This is what makes the pre-registered prediction -- **HCC-on-HCC > HCC-on-AD** -- a
    TEST rather than a comparison of two significance verdicts. Testing each arm against
    its own controls and concluding from "one is significant and the other is not" is the
    difference-of-significance error; the contrast has to be the thing permuted. The same
    ligand carries an M on both panels, so the difference is paired within ligand, every
    ligand-level nuisance (size, flexibility, the molecule's general stickiness) cancels in
    it, and the SAME within-stratum label permutation distributes it exactly as it does
    each arm's own statistic."""
    return {lid: m_primary[lid] - m_control[lid]
            for lid in m_primary if lid in m_control}


def _mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return float("nan")
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def _sd(xs):
    xs = [float(x) for x in xs]
    n = len(xs)
    if n < 2:
        return float("nan")
    mu = sum(xs) / n
    return math.sqrt(sum((x - mu) ** 2 for x in xs) / (n - 1))


def balance_rows(detail, ids, positive, stratum, heavy_atoms, n_targets) -> tuple:
    """(per-stratum balance rows, weighted standardised difference, max |standardised|).

    **Why this exists.** Stratification removes the BETWEEN-stratum size confound; it does
    not remove the WITHIN-stratum one. With the real 1,960-record library and scores that
    are a pure monotone function of heavy-atom count and carry no binding signal at all,
    positives drawn uniformly inside the big strata come out correctly null (median
    p 0.57-0.59) -- but positives taken as the 25 largest molecules give statistic +0.10,
    p = 0.0001, and a within-stratum heavy-atom imbalance of +0.62 sd already reaches
    P(p <= 0.05) = 0.61, +0.79 sd reaches 0.99 and >= +1.15 sd reaches 1.00. The residual
    size-driven component of M is about 0.05 sd in EVERY stratum (the slope dM/dHA, 0.0073
    to 0.0410 per atom, and the within-stratum sd of heavy atoms, 7.1 down to 1.0, trade
    off so the product is flat across all 15 cells), which is 2-3x the effect this test
    detects at 67 % power. HA bin 4 is open-ended and spans 32-63 heavy atoms.

    The test itself is NOT changed by this: that would be tuning the analysis to the
    objection. What is added is the diagnostic that says whether the objection applies to
    this particular positive set, which nothing previously reported. For the plausible real
    HCC set (sorafenib 32, lenvatinib 30, regorafenib 33, cabozantinib 37, apatinib 40,
    doxorubicin 39, sunitinib 29, erlotinib 29) the within-stratum position is near the
    middle and in HA4-RB2 BELOW it, and under the pure-size generator that set gives
    statistic -0.0041, p = 0.59. The real analysis is probably safe -- but by accident of
    which drugs these are, not by construction."""
    by: dict = {}
    for lid, pos, s in zip(ids, positive, stratum):
        if s not in detail:
            continue
        p, c = by.setdefault(s, ([], []))
        (p if pos else c).append(lid)
    rows, weighted, biggest = [], 0.0, 0.0
    for s in sorted(by):
        pids, cids = by[s]
        hp = [float(heavy_atoms.get(i, float("nan"))) for i in pids]
        hc = [float(heavy_atoms.get(i, float("nan"))) for i in cids]
        sd = _sd(hc)
        std = ((_mean(hp) - _mean(hc)) / sd) if sd and math.isfinite(sd) else float("nan")
        w = detail[s]["weight"]
        if math.isfinite(std):
            weighted += w * std
            biggest = max(biggest, abs(std))
        rows.append({
            "stratum": s, "weight": round(w, 6),
            "n_positive": len(pids), "n_control": len(cids),
            "ha_mean_positive": _num(_mean(hp)), "ha_mean_control": _num(_mean(hc)),
            "ha_median_positive": _num(_median(hp)),
            "ha_median_control": _num(_median(hc)),
            "ha_sd_control": _num(sd), "ha_standardised_difference": _num(std),
            "targets_mean_positive": _num(_mean(n_targets.get(i, float("nan"))
                                                for i in pids)),
            "targets_mean_control": _num(_mean(n_targets.get(i, float("nan"))
                                               for i in cids)),
            "flag": ("size imbalance >= %.2g sd" % BALANCE_FLAG_SD)
                    if math.isfinite(std) and abs(std) >= BALANCE_FLAG_SD else ""})
    return rows, weighted, biggest


def balance_flag(weighted, biggest) -> str:
    """The one-line verdict that goes on the arm row and into summary.md."""
    if not math.isfinite(weighted):
        return "balance not computable (no heavy-atom counts)"
    if abs(weighted) >= BALANCE_FLAG_SD:
        return (f"SIZE IMBALANCE: the positives sit {weighted:+.2f} sd from their matched "
                f"controls in heavy atoms WITHIN strata; an imbalance of this size can "
                f"produce p < 0.05 from molecular size alone (see summary.md)")
    if biggest >= 1.0:
        return (f"one stratum is imbalanced by {biggest:.2f} sd in heavy atoms (the "
                f"weighted imbalance is only {weighted:+.2f} sd); check balance.csv")
    return ""


def holm(pvalues) -> list:
    """Holm-Bonferroni adjusted p-values, order preserved, monotone, clipped at 1."""
    vals = [float(p) for p in pvalues]
    m = len(vals)
    order = sorted(range(m), key=lambda i: vals[i])
    out = [0.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * vals[i])
        out[i] = min(1.0, running)
    return out


def benjamini_hochberg(pvalues) -> list:
    """BH FDR q-values, order preserved, monotone, clipped at 1."""
    vals = [float(p) for p in pvalues]
    m = len(vals)
    order = sorted(range(m), key=lambda i: vals[i])
    out = [0.0] * m
    running = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        running = min(running, m * vals[i] / (rank + 1))
        out[i] = min(1.0, running)
    return out


def single_drug_p(value: float, controls) -> float:
    """Exact one-sided p for ONE drug against the controls of its own stratum:
    (1 + #{control M >= drug M}) / (1 + n_control). Exploratory only."""
    controls = list(controls)
    if not controls:
        return float("nan")
    return (1 + sum(1 for c in controls if c >= value)) / (1 + len(controls))


# --- step 7: the descriptive surface ---------------------------------------------------

def check_no_pvalues(where, columns, rows, allowed=None) -> None:
    """Refuse anything that is not an expected descriptive-surface field.

    This is an ALLOWLIST, not a denylist: `allowed` (SURFACE_COLUMNS by default) is the
    complete set of field names the descriptive surface may carry, in the column list AND
    in any row dict, and ANY other name is refused -- whether or not it looks like a
    p-value. A denylist on free text cannot be made tight: `p-value` (the commonest
    spelling), `adj.P.Val` (R/limma), `Pr(>|t|)`, `pFDR`, `log10p` and a Cyrillic `р_value`
    all walked straight through the bounded-token pattern this replaced.

    PVALUE_LIKE is kept as a second line, applied to the allowlist itself, so that an
    expected column can never be renamed into a p-value without this failing first."""
    allow = set(SURFACE_COLUMNS if allowed is None else allowed)
    names = list(columns) + sorted({k for r in rows for k in r})
    unexpected = sorted({str(c) for c in names if str(c) not in allow})
    forbidden = sorted({str(c) for c in allow if pvalue_like(c)})
    if unexpected or forbidden:
        bad = unexpected + forbidden
        raise ValueError(
            f"{where}: the descriptive surface carries no p-values, by design "
            f"(spec 2026-11-03-statistics-design.md step 7, pre-registration section 8). "
            f"It may carry ONLY the fields {', '.join(sorted(allow))}; refusing "
            f"{', '.join(bad)}")


def write_no_pvalues(path, columns, rows, allowed=None) -> None:
    """Write a CSV that CANNOT carry a p-value (see `check_no_pvalues`)."""
    check_no_pvalues(path, columns, rows, allowed)
    _write_csv(path, columns, rows)


def surface_rows(panel: str, drug_set: str, counts_positive: dict, counts_control: dict,
                 thresholds=SURFACE_THRESHOLDS, min_targets=SURFACE_MIN_TARGETS) -> list:
    """The threshold x min_targets grid, descriptively.

    `counts_positive` / `counts_control` are {ligand ID: [score on each target]}."""
    rows = []
    n_pos, n_ctl = len(counts_positive), len(counts_control)
    for t in thresholds:
        hits_pos = {lid: sum(1 for s in v if s <= t) for lid, v in counts_positive.items()}
        hits_ctl = {lid: sum(1 for s in v if s <= t) for lid, v in counts_control.items()}
        for k in min_targets:
            mp = sum(1 for h in hits_pos.values() if h >= k)
            mc = sum(1 for h in hits_ctl.values() if h >= k)
            rows.append({
                "panel": panel, "drug_set": drug_set,
                "score_threshold": f"{t:g}", "min_targets": k,
                "n_positive": n_pos, "n_positive_meeting": mp,
                "frac_positive_meeting": round(mp / n_pos, 6) if n_pos else "",
                "n_control": n_ctl, "n_control_meeting": mc,
                "frac_control_meeting": round(mc / n_ctl, 6) if n_ctl else ""})
    return rows


def surface_markdown(rows) -> str:
    """The same grid as a markdown table, with the "no p-values" statement on it.

    Guarded by the SAME allowlist as surface.csv: the markdown is what gets pasted into a
    draft, so it must not be the unguarded way in."""
    check_no_pvalues(SURFACE_MD, [], rows)
    out = ["# Descriptive surface: score threshold x min_targets",
           "",
           "DESCRIPTIVE ONLY. **No p-values are attached to this surface and none may be.**",
           "It is how the field states these results and it should be shown; the test is "
           "the permutation test in results.csv (see "
           "docs/superpowers/specs/2026-11-03-statistics-design.md, step 7).",
           "",
           f"Controls: {CONTROL_SET}.",
           ""]
    by_set = {}
    for r in rows:
        by_set.setdefault((r["panel"], r["drug_set"]), []).append(r)
    for (panel, drug_set), rs in sorted(by_set.items()):
        out.append(f"## {drug_set} on {panel}")
        out.append("")
        ks = sorted({int(r["min_targets"]) for r in rs})
        out.append("| score <= | " + " | ".join(f"n>={k}" for k in ks) + " |")
        out.append("|---" * (len(ks) + 1) + "|")
        for t in sorted({float(r["score_threshold"]) for r in rs}):
            cells = []
            for k in ks:
                hit = next((r for r in rs if float(r["score_threshold"]) == t
                            and int(r["min_targets"]) == k), None)
                if hit is None:
                    cells.append("")
                    continue
                cell = f'{hit["n_positive_meeting"]}/{hit["n_positive"]}'
                if hit["frac_control_meeting"] != "":
                    cell += f' vs {float(hit["frac_control_meeting"]):.3f}'
                cells.append(cell)
            out.append(f"| {t:g} | " + " | ".join(cells) + " |")
        out.append("")
        out.append("Cells are `positives meeting / positives` vs the matched-control "
                   "fraction meeting the same rule.")
        out.append("")
    return "\n".join(out) + "\n"


# --- inputs ----------------------------------------------------------------------------

def _read_csv(path) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _write_csv(path, columns, rows) -> None:
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(columns))
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in columns})
    os.replace(tmp, path)


def _atomic_text(path, text: str) -> None:
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def read_library(path) -> list:
    """The dockable library records from the `ligands` stage's ligands.csv:
    [{ID, name, heavy_atoms, rotatable_bonds}], kept records only, in file order.

    The rotatable-bond count is the pinned non-standard one (`ligands.profile`); it is
    never recomputed here, and RDKit's Lipinski count must not be substituted -- it
    differs for 1,316 of the 2,083 e-Drug3D records."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found: run `polyphar ligands` first (the strata are built from "
            f"its pinned heavy-atom / rotatable-bond counts)")
    out = []
    for r in _read_csv(path):
        if str(r.get("kept", "")).strip().lower() not in ("true", "1", "yes"):
            continue
        try:
            out.append({"ID": str(r["ID"]).strip(), "name": str(r.get("name", "")).strip(),
                        "heavy_atoms": int(float(r["heavy_atoms"])),
                        "rotatable_bonds": int(float(r["rotatable_bonds"]))})
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{path}: unusable row for ID "
                             f"{r.get('ID')!r}: {exc}") from None
    if not out:
        raise ValueError(f"{path}: no kept ligands")
    return out


def read_positive_sets(path, library) -> tuple:
    """({set name: {ID: {name, subgroup, source}}}, [(item, reason)] for unresolved rows).

    Schema (header required): set,ID,name,subgroup,source
      set       'hcc' or 'ad' (any other value is carried as an extra exploratory set)
      ID        the e-Drug3D library ID; optional if `name` is given
      name      matched case-insensitively against the library's `name` when ID is empty
      subgroup  'targeted' or 'cytotoxic' (HCC set; the halves carry opposite size biases
                and pooling them hides that)
      source    free text: the guideline the drug came from

    A row naming a drug that is NOT in the dockable library is NOT silently dropped: it
    comes back in the unresolved list, is written to failed.tsv and is named in summary.md,
    because "which guideline drug was dropped for undockability" is a number the Methods
    must state."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(POSITIVE_SET_MISSING.format(path=path))
    by_id = {str(r["ID"]): r for r in library}
    by_name: dict = {}
    for r in library:
        by_name.setdefault(r["name"].strip().lower(), []).append(r)
    rows = _read_csv(path)
    if not rows:
        with open(path, newline="", encoding="utf-8") as fh:
            header = csv.DictReader(fh).fieldnames or []
        if "set" in header:
            raise ValueError(f"{path}: the header is there but the file has no rows. The "
                             f"positive set is the clinical-guideline list; an empty one "
                             f"is not an analysis with no positives, it is a missing "
                             f"input")
        raise ValueError(f"{path}: needs a header with at least the columns "
                         f"{', '.join(POSITIVE_SET_COLUMNS)}; found "
                         f"{', '.join(header) if header else 'an empty file'}")
    if "set" not in rows[0]:
        raise ValueError(f"{path}: needs a header with at least the columns "
                         f"{', '.join(POSITIVE_SET_COLUMNS)}; found "
                         f"{', '.join(rows[0])}")
    sets: dict = {}
    unresolved = []
    for i, r in enumerate(rows, start=2):
        name_in = str(r.get("name", "")).strip()
        where = f"{path.name}:{i} {name_in or r.get('ID') or '?'}"
        group = str(r.get("set", "")).strip().lower()
        if not group:
            unresolved.append((where, "no `set` value"))
            continue
        lid = str(r.get("ID", "")).strip()
        if lid and lid in by_id:
            rec = by_id[lid]
        elif lid:
            unresolved.append((where, f"ID {lid} is not a dockable library record "
                                      f"(excluded by the ligands stage, or not in the "
                                      f"library at all)"))
            continue
        else:
            hits = by_name.get(name_in.lower(), [])
            if len(hits) != 1:
                unresolved.append((where, f"name {name_in!r} matches {len(hits)} dockable "
                                          f"library records; give the ID"))
                continue
            rec = hits[0]
        if name_in and rec["name"].strip().lower() != name_in.lower():
            unresolved.append((where, f"ID {rec['ID']} is {rec['name']!r} in the library, "
                                      f"not {name_in!r}; taken as {rec['name']!r}"))
        sets.setdefault(group, {})[rec["ID"]] = {
            "name": rec["name"], "subgroup": str(r.get("subgroup", "")).strip().lower(),
            "source": str(r.get("source", "")).strip()}
    return sets, unresolved


def read_targets(path) -> dict:
    """{label: UniProt accession} from a final_targets.csv."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; run `polyphar select` first")
    out = {}
    for r in _read_csv(path):
        out[f'{r["Protein"]}_{r["Pocket"]}'] = r["UniProt_ID"]
    if not out:
        raise ValueError(f"{path}: no targets")
    return out


def df_results_path(docking: Path, label: str) -> Path:
    return Path(docking) / f"df_results_{label}.csv"


def read_df_results(path) -> tuple:
    """({ligand ID: score_mean}, [(ID, why)]) from one df_results_<label>.csv.

    Unusable rows are COUNTED AND RETURNED, never silently skipped: a file that exists but
    whose every row is unparseable used to yield an empty dict, which looked exactly like
    "this target has no effect" all the way through to a blank statistic and a blank p in
    results.csv.

    `nan` and `inf` are rejected rather than accepted. `percentiles` sorts with NaN
    present, and NaN compares false against everything, so one NaN does not merely rank
    itself wrongly -- it reorders its neighbours ([-10, nan, -5, -8] comes out
    [0.875, 0.625, 0.125, 0.375], with -5 and -8 swapped) and corrupts that target's whole
    ranking. dock_unidock filters non-finite poses; the smina path has no equivalent."""
    out, dropped = {}, []
    for r in _read_csv(path):
        lid = str(r.get("ID", "")).strip()
        raw = r.get("score_mean", None)
        try:
            v = float(raw)
        except (TypeError, ValueError):
            dropped.append((lid or "?", f"score_mean {raw!r} is not a number"))
            continue
        if not math.isfinite(v):
            dropped.append((lid or "?", f"score_mean {raw!r} is not finite"))
            continue
        if not lid:
            dropped.append(("?", "no ID"))
            continue
        out[lid] = v
    return out, dropped


def replicate_paths(docking: Path, label: str) -> list:
    """[(replicate index, sdf)] for one target, replicate order."""
    d = Path(docking) / label / "replicates"
    out = []
    for p in sorted(d.glob("r*.sdf")):
        m = re.fullmatch(r"r(\d+)", p.stem)
        if m and p.with_suffix(".done").is_file():
            out.append((int(m.group(1)), p))
    return sorted(out)


def read_replicate_best(path) -> dict:
    """{ligand ID: best (lowest) pose score} from one replicate's raw SDF.

    Streamed line by line: a replicate file holds every pose of every ligand and the whole
    panel's raw output is tens of gigabytes. The score field is the engine's
    (`minimizedAffinity` for smina, `docking_score` for unidock_tools); a record carrying
    both is refused rather than silently resolved, exactly as dock._pose_score does."""
    from polyphar.stages.dock import SCORE_FIELDS, TITLE_PREFIX

    best: dict = {}
    lid = None
    scores: dict = {}
    title = None
    want = None
    index = 0                       # line index inside the current record
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.rstrip("\r\n")
            index += 1
            if index == 1:          # the title line, exactly: the smina fallback ID
                title = s.strip()
                continue
            if want is not None:
                if want == "ID":
                    lid = s.strip()
                else:
                    try:
                        scores[want] = float(s.strip())
                    except ValueError:
                        pass
                want = None
                continue
            if s.startswith(">"):
                m = re.match(r"^>.*<([^>]+)>", s)
                if m and (m.group(1) == "ID" or m.group(1) in SCORE_FIELDS):
                    want = m.group(1)
                continue
            if s == "$$$$":
                key = lid or (title[len(TITLE_PREFIX):]
                              if title and title.startswith(TITLE_PREFIX) else None)
                present = [f for f in SCORE_FIELDS if f in scores]
                if len(present) > 1:
                    raise ValueError(f"{path}: a pose of {key} carries more than one score "
                                     f"field ({', '.join(present)})")
                if key and present:
                    v = scores[present[0]]
                    if key not in best or v < best[key]:
                        best[key] = v
                lid, scores, title, index = None, {}, None, 0
                continue
    return best


# --- pool workers ----------------------------------------------------------------------

def _worker_df(path):
    return pool.task(lambda: read_df_results(path))


def _worker_replicate(path):
    return pool.task(lambda: read_replicate_best(path))


# --- panels -----------------------------------------------------------------------------

class Panel:
    """One disease panel's docking results, loaded from a finished run folder."""

    def __init__(self, name, run_dir, protein_of, pooled, replicates):
        self.name = name
        self.run_dir = Path(run_dir)
        self.protein_of = protein_of                 # label -> accession
        self.pooled = pooled                         # label -> {ID: score_mean}
        self.replicates = replicates                 # replicate -> {label: {ID: best}}

    @property
    def n_targets(self) -> int:
        return len(self.pooled)

    @property
    def n_proteins(self) -> int:
        return len({self.protein_of[l] for l in self.pooled if l in self.protein_of})

    @property
    def accessions(self) -> set:
        return {self.protein_of[l] for l in self.pooled if l in self.protein_of}


def panel_jobs(run_dir) -> tuple:
    """(protein_of, [(item, fn, args)], {item: (kind, run key, label, replicate, path)})
    for one run folder: one job per df_results file and one per finished replicate SDF.

    The run folder's path is part of every item and of every meta key, so two panels
    sharing a target label (the same protein and pocket id in both diseases) can never
    have one panel's scores read as the other's."""
    run_dir = Path(run_dir)
    key = str(run_dir)
    protein_of = read_targets(run_dir / "06_targets" / "final_targets.csv")
    docking = run_dir / "07_docking"
    jobs, meta = [], {}
    for label in sorted(protein_of):
        p = df_results_path(docking, label)
        if p.is_file():
            item = f"{key}|pooled|{label}"
            jobs.append((item, _worker_df, (str(p),)))
            meta[item] = ("pooled", key, label, None, p)
        for i, sdf in replicate_paths(docking, label):
            item = f"{key}|r{i}|{label}"
            jobs.append((item, _worker_replicate, (str(sdf),)))
            meta[item] = ("replicate", key, label, i, sdf)
    return protein_of, jobs, meta


# --- the positive-control gate ------------------------------------------------------------

def positive_control_verdict(path, display=None) -> tuple:
    """(verdict, detail) for the pre-registered positive-control gate.

    The spec pre-registers that a FAILED positive control makes the primary result
    UNINTERPRETABLE rather than reinterpreted, and that an absent one is not a pass. Both
    are enforced here. JSON schema:
        {"median_self_rank_percentile": 0.02,     # <= 0.05 to pass, over >= 10 ligands
         "redock_rmsd_success_fraction": 0.74,    # >= 0.50 to pass, over >= 10 self pairs
         "library_top_decile_fraction": 0.80,     # >= 0.50 to pass, over >= 10 ligands
         "n_redocked": 24, "n_redock_rmsd": 24, "n_library_percentile": 24,
         "panel_targets_scored": 478, "panel_targets_total": 482,   # ratio >= 0.98
         "pose_in_pocket_fraction": 0.83,         # DESCRIPTIVE, range-checked, NOT gated
         "source": "..."}                         # how the numbers were produced

    Every readout is a fraction and must be a real number in [0, 1]. A boolean, a
    non-finite value, a negative rank or a fraction of 99 is NOT a passing control and NOT
    a near-miss: it comes back as `unknown: ... is out of range`, because an invalid
    control must land in `unknown` and never in `yes`.

    A readout computed over too few observations is `unknown`, not `yes` and not `no`: one
    labelled ligand ranking first on a 20-target panel is not evidence either way. So is a
    control whose panel coverage is below MIN_PANEL_SCORED_FRACTION, because an excluded
    target is excluded from every ligand's denominator and attrition correlated with score
    moves the median the helpful way.

    A file with no `source` still passes on its numbers, but the omission is written to
    failed.tsv (see `run`): the gate cannot check provenance, only say that none was given.

    `display` is how the file is named in the verdict (the option's own value, so the
    result rows stay readable); it defaults to `path`."""
    path = Path(path)
    shown = display or path
    if not path.is_file():
        return ("unknown: positive control not run "
                f"({shown} absent); the primary result is NOT interpretable until it is",
                {})
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return f"unknown: {path} is unreadable ({exc})", {}
    detail = dict(data) if isinstance(data, dict) else {}
    if not isinstance(data, dict):
        return f"unknown: {path} is not a JSON object", {}
    values = {}
    for key in ("median_self_rank_percentile", "redock_rmsd_success_fraction",
                "library_top_decile_fraction", "pose_in_pocket_fraction"):
        raw = data.get(key, None)
        # `True` is an int in Python and float(True) == 1.0; a boolean is not a measured
        # fraction and must not be read as a passing one.
        if raw is None or isinstance(raw, bool) or isinstance(raw, (list, dict)):
            return (f"unknown: {path} lacks a numeric {key}", detail)
        try:
            v = float(raw)
        except (TypeError, ValueError):
            return (f"unknown: {path} lacks a numeric {key}", detail)
        if not math.isfinite(v):
            return (f"unknown: {path} is out of range ({key} = {raw!r} is not finite); "
                    f"an invalid positive control is NOT a pass", detail)
        if not 0.0 <= v <= 1.0:
            return (f"unknown: {path} is out of range ({key} = {v:.6g}, and every readout "
                    f"is a fraction in [0, 1]); an invalid positive control is NOT a pass",
                    detail)
        values[key] = v
    # Too little to judge on: `unknown`, which the spec already treats as not-interpretable
    # rather than as a pass or a failure.
    thin = []
    for key, minimum, what in (
            ("n_redocked", MIN_SELF_RANK_N, "self-rank"),
            ("n_redock_rmsd", MIN_REDOCK_N, "redocking RMSD"),
            ("n_library_percentile", MIN_LIBRARY_N, "library percentile")):
        n = data.get(key, None)
        if n is None or isinstance(n, bool) or not isinstance(n, (int, float)):
            thin.append(f"{key} is missing, so the {what} readout cannot be sized")
        elif n < minimum:
            thin.append(f"{what} is over {n:g} observation(s), below the pre-registered "
                        f"minimum of {minimum}")
    total = data.get("panel_targets_total", None)
    scored = data.get("panel_targets_scored", None)
    if not all(isinstance(v, (int, float)) and not isinstance(v, bool)
               for v in (total, scored)) or not total:
        thin.append("panel_targets_scored / panel_targets_total are missing, so the "
                    "self-rank denominator cannot be checked for target attrition")
    elif scored / total < MIN_PANEL_SCORED_FRACTION:
        thin.append(f"only {scored:g} of {total:g} panel targets were scored "
                    f"({scored / total:.4g}), below the pre-registered minimum "
                    f"{MIN_PANEL_SCORED_FRACTION:g}; an excluded target is excluded from "
                    f"every ligand's denominator")
    if thin:
        return (f"unknown: {shown} does not support a verdict ({'; '.join(thin)}); the "
                f"primary result is NOT interpretable until it does"), detail

    rank = values["median_self_rank_percentile"]
    redock = values["redock_rmsd_success_fraction"]
    lib = values["library_top_decile_fraction"]
    fails = []
    if not rank <= SELF_RANK_MAX:
        fails.append(f"median self-rank {rank:.6g} is not within the top "
                     f"{SELF_RANK_MAX:g}")
    if not redock >= REDOCK_SUCCESS_MIN:
        fails.append(f"redocking success (RMSD <= {REDOCK_RMSD_MAX:g} A) {redock:.6g} is "
                     f"below {REDOCK_SUCCESS_MIN:g}")
    if not lib >= LIBRARY_TOP_DECILE_MIN:
        fails.append(f"only {lib:.6g} of the co-crystal ligands reach the top decile of "
                     f"the production library at their own target, below "
                     f"{LIBRARY_TOP_DECILE_MIN:g}")
    if fails:
        return ("no: " + "; ".join(fails) + " -- the primary result is reported as "
                "UNINTERPRETABLE, not reinterpreted (pre-registered)"), detail
    return "yes", detail


# --- arms ---------------------------------------------------------------------------------

def _arm_rows(panel, members, drug_set, family, m_pooled, m_by_replicate, strata, args,
              interpretable, note, sizes=None, targets=None, stratum_key="of",
              arm_name=None, direction=DIRECTION, balance_out=None) -> dict:
    """One arm: the permutation test on the pooled M plus the same statistic per replicate.

    `stratum_key` selects the stratification: "of" is the pre-registered primary one
    (heavy-atom quintiles x rotatable-bond tertiles); "fine_of" is the pre-registered
    sensitivity analysis (deciles x tertiles)."""
    ids = [lid for lid in strata["order"] if lid in m_pooled]
    m = [m_pooled[lid] for lid in ids]
    positive = [lid in members for lid in ids]
    stratum = [strata[stratum_key][lid] for lid in ids]
    res = permutation_test(m, positive, stratum,
                           getattr(args, "permutations", N_PERMUTATIONS),
                           getattr(args, "permutation_seed", SEED))
    arm = arm_name or f"{drug_set} on {panel.name}"
    notes = [note] if note else []

    # positives that entered a stratum holding no control contribute nothing and used to
    # vanish without a word; `n_positive` is the number that ENTERED the test, which is
    # why it can disagree with len(members).
    entered = res["n_positive"]
    in_panel = sum(1 for lid in ids if lid in members)
    if in_panel > entered:
        notes.append(f"{in_panel - entered} of the {in_panel} positive(s) carrying an M "
                     f"sat in a stratum with no matched control and contributed nothing "
                     f"to this arm")
    if 0 < entered < 3:
        notes.append(f"n_positive = {entered}: with so few positives the permutation p is "
                     f"at or near its floor whatever the data say; read it as descriptive")

    per_rep = []
    for i in sorted(m_by_replicate):
        mr = m_by_replicate[i]
        rids = [lid for lid in ids if lid in mr]
        stat, _ = stratified_difference([mr[lid] for lid in rids],
                                        [lid in members for lid in rids],
                                        [strata[stratum_key][lid] for lid in rids])
        per_rep.append((i, stat))
    finite = [s for _i, s in per_rep if math.isfinite(s)]
    spread = (max(finite) - min(finite)) if len(finite) > 1 else ""
    # NOT corroboration. corr(M_r0, M_r1) = 0.9996 on real data, so this is True on every
    # arm this stage has ever produced, including p = 0.45 and p = 0.80. It can only catch
    # a catastrophic jitter regime in which replicates disagree about the sign at all.
    agrees = ("" if not finite else
              all((s > 0) == (res["statistic"] > 0) for s in finite))
    if len(per_rep) > 3:
        notes.append(f"{len(per_rep)} replicates were read; only r0-r2 have a column "
                     f"(stat_r3+ = "
                     + ", ".join(f"r{i}={_num(s)}" for i, s in per_rep[3:]) + ")")

    brows, weighted, biggest = ([], float("nan"), float("nan"))
    if sizes is not None:
        brows, weighted, biggest = balance_rows(res["strata"], ids, positive, stratum,
                                                sizes, targets or {})
        if balance_out is not None:
            for b in brows:
                balance_out.append({"arm": arm, "panel": panel.name,
                                    "drug_set": drug_set, "family": family, **b})
    flag = balance_flag(weighted, biggest) if sizes is not None else ""
    if flag:
        notes.append(flag)

    row = {"arm": arm, "drug_set": drug_set,
           "panel": panel.name, "family": family,
           "n_positive": entered, "n_control": res["n_control"],
           "n_strata": res["n_strata"],
           "statistic": _num(res["statistic"]), "p_value": _num(res["p_value"]),
           "p_adjusted": "", "adjustment": "",
           "direction": direction,
           "replicate_spread": _num(spread), "replicate_sign_jitter_check": agrees,
           "balance_sd": _num(weighted), "balance_sd_max_stratum": _num(biggest),
           "balance_flag": flag,
           "interpretable": interpretable, "control_set": CONTROL_SET,
           "note": "; ".join(notes)}
    for i, stat in per_rep[:3]:
        row[f"stat_r{i}"] = _num(stat)
    return row


def _num(x):
    if x == "" or x is None:
        return ""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return ""
    return "" if not math.isfinite(v) else round(v, 6)


# --- stage ---------------------------------------------------------------------------------

def stage_dir(layout) -> Path:
    return layout.run_dir / DIRNAME


def _path(layout, value) -> Path:
    p = Path(value)
    return p if p.is_absolute() else layout.run_dir / p


def _inputs(layout, args) -> dict:
    return {"library": layout.docking / "library" / "ligands.csv",
            "targets": layout.final_targets_csv,
            "positive_set": _path(layout, getattr(args, "positive_set",
                                                  POSITIVE_SETS)),
            "positive_control": _path(layout, getattr(args, "positive_control",
                                                      POSITIVE_CONTROL))}


def _params(args) -> dict:
    return {"permutations": int(getattr(args, "permutations", N_PERMUTATIONS)),
            "permutation_seed": int(getattr(args, "permutation_seed", SEED)),
            "min_coverage": float(getattr(args, "min_coverage", MIN_COVERAGE)),
            "ha_quantiles": list(HA_QUANTILES), "rb_quantiles": list(RB_QUANTILES),
            "ha_quantiles_fine": list(HA_QUANTILES_FINE),
            "balance_flag_sd": BALANCE_FLAG_SD,
            "statistic": "mean over proteins of max over that protein's targets of the "
                         "within-target percentile",
            "percentile": "mid-rank; 1 = most negative score_mean",
            "direction": DIRECTION,
            "control_panel": str(getattr(args, "control_panel", None) or ""),
            **THRESHOLDS}


def _signature(layout, args, jobs_meta) -> str:
    """What the result depends on: the parameters, the small inputs' sha256, every
    df_results file's sha256 and every replicate file's (path, size).

    The replicate SDFs are identified by size rather than hashed: the panel's raw poses are
    tens of gigabytes, and their CONTENT is already covered by the df_results hashes, which
    `dock` derives from them."""
    parts = {"params": _params(args)}
    for name, path in _inputs(layout, args).items():
        parts[name] = sha256_file(path) if Path(path).is_file() else ""
    files = []
    for _item, (kind, _key, _label, _rep, path) in sorted(jobs_meta.items()):
        if kind == "pooled":
            files.append([str(path), sha256_file(path)])
        else:
            try:
                files.append([str(path), Path(path).stat().st_size])
            except OSError:
                files.append([str(path), -1])
    parts["files"] = files
    return hashlib.sha256(json.dumps(parts, sort_keys=True,
                                     default=str).encode("utf-8")).hexdigest()


def _finished(layout, args, jobs_meta) -> dict | None:
    try:
        manifest = json.loads((stage_dir(layout) / "manifest.json").read_text())
    except (OSError, ValueError):
        return None
    if not isinstance(manifest, dict) or manifest.get("signature") != _signature(
            layout, args, jobs_meta):
        return None
    counts = manifest.get("counts")
    return dict(counts) if isinstance(counts, dict) else {}


def _all_jobs(layout, args):
    """(primary panel jobs/meta, control panel jobs/meta or None)."""
    protein_of, jobs, meta = panel_jobs(layout.run_dir)
    control = None
    cdir = getattr(args, "control_panel", None)
    if cdir:
        cpath = Path(cdir)
        # the negative control is "the same drugs on ANOTHER disease's panel". Pointing it
        # at this run folder renames the second panel `<name>_control` and Holm-corrects
        # two identical p-values, which fabricates a negative control out of nothing.
        if cpath.resolve() == Path(layout.run_dir).resolve():
            raise ValueError(
                f"--control-panel {cdir} is this run folder. The pre-registered negative "
                f"control is the SAME drug set on a DIFFERENT disease's panel; running a "
                f"panel against itself produces two identical p-values and a negative "
                f"control that cannot fail")
        cprotein, cjobs, cmeta = panel_jobs(cpath)
        control = (cpath, cprotein, cjobs, cmeta)
    return (protein_of, jobs, meta), control


def units(layout, args) -> WorkUnits:
    """One unit per file read (df_results and finished replicate SDFs); 0 once the stage's
    signature is unchanged."""
    inputs = _inputs(layout, args)
    if not inputs["positive_set"].is_file():
        raise FileNotFoundError(POSITIVE_SET_MISSING.format(path=inputs["positive_set"]))
    (_p, jobs, meta), control = _all_jobs(layout, args)
    all_meta = dict(meta)
    if control is not None:
        jobs = jobs + control[2]
        all_meta.update(control[3])
    if _finished(layout, args, all_meta) is not None:
        return WorkUnits(0, "target reads")
    sizes = []
    for item, _fn, _a in jobs:
        try:
            sizes.append(float(all_meta[item][4].stat().st_size))
        except OSError:
            sizes.append(0.0)
    return WorkUnits(len(jobs), "target reads", sizes)


def _load(jobs, meta, args, progress, failures, dropped=None) -> tuple:
    """Run the read jobs. ({(run key, label): {ID: score}},
    {(run key, replicate): {label: {ID: best}}}). Order-independent, so `--nprocs` cannot
    change the result.

    `dropped` collects (item, why) for every unusable df_results row, so an unparseable
    score is reported rather than absorbed into an empty table."""
    pooled, replicates = {}, {}
    sizes = {}
    for item, (_k, _key, _l, _r, path) in meta.items():
        try:
            sizes[item] = float(path.stat().st_size)
        except OSError:
            sizes[item] = 0.0

    def on_result(item, result):
        kind, key, label, rep, path = meta[item]
        if kind == "pooled":
            scores, bad = result
            pooled[(key, label)] = scores
            if bad and dropped is not None:
                shown = "; ".join(f"{lid}: {why}" for lid, why in bad[:5])
                dropped.append((str(path), f"{len(bad)} unusable score_mean row(s) "
                                           f"dropped ({shown}"
                                           f"{', ...' if len(bad) > 5 else ''})",
                                len(bad)))
        else:
            replicates.setdefault((key, rep), {})[label] = result
        if progress is not None:
            progress.advance(1, size=sizes.get(item, 0.0))

    def on_failure(item, why):
        failures.append((item, why))
        if progress is not None:
            progress.advance(1, size=sizes.get(item, 0.0))

    pool.run_pool(jobs, max(1, getattr(args, "nprocs", 1)),
                  on_result=on_result, on_failure=on_failure)
    return pooled, replicates


def build_strata(library, args) -> dict:
    """The size strata: bin edges from the library, bins per ligand, realised sizes.

    Called BEFORE any docking score is read -- see `run`, where strata.json is written
    first. That ordering is the point: the edges cannot be tuned to the outcome if they
    exist before the outcome is loaded."""
    ha = [r["heavy_atoms"] for r in library]
    rb = [r["rotatable_bonds"] for r in library]
    edges = {"heavy_atoms": quantile_edges(ha, HA_QUANTILES),
             "rotatable_bonds": quantile_edges(rb, RB_QUANTILES)}
    fine_edges = {"heavy_atoms": quantile_edges(ha, HA_QUANTILES_FINE),
                  "rotatable_bonds": edges["rotatable_bonds"]}
    of, bins, sizes = {}, {}, {}
    fine_of, fine_sizes = {}, {}
    for r in library:
        hb = bin_of(r["heavy_atoms"], edges["heavy_atoms"])
        rbb = bin_of(r["rotatable_bonds"], edges["rotatable_bonds"])
        s = stratum_label(hb, rbb)
        of[r["ID"]] = s
        bins[r["ID"]] = (hb, rbb)
        sizes[s] = sizes.get(s, 0) + 1
        fs = stratum_label(bin_of(r["heavy_atoms"], fine_edges["heavy_atoms"]), rbb)
        fine_of[r["ID"]] = fs
        fine_sizes[fs] = fine_sizes.get(fs, 0) + 1
    matches = (edges["heavy_atoms"] == PREREGISTERED_EDGES["heavy_atoms"]
               and edges["rotatable_bonds"] == PREREGISTERED_EDGES["rotatable_bonds"])
    return {"edges": edges, "of": of, "bins": bins,
            "sizes": {k: sizes[k] for k in sorted(sizes)},
            "order": [r["ID"] for r in library],
            # the pre-registered SENSITIVITY stratification, not the primary one
            "fine_edges": fine_edges, "fine_of": fine_of,
            "fine_sizes": {k: fine_sizes[k] for k in sorted(fine_sizes)},
            "fine_quantiles": {"heavy_atoms": list(HA_QUANTILES_FINE),
                               "rotatable_bonds": list(RB_QUANTILES)},
            "fine_purpose": "pre-registered sensitivity analysis: heavy-atom deciles x "
                            "rotatable-bond tertiles, because the primary quintile bins "
                            "are wide (HA bin 4 is open-ended, 32-63 heavy atoms) and "
                            "stratification removes the between-stratum size confound, "
                            "not the within-stratum one",
            "preregistered": PREREGISTERED_EDGES,
            "preregistered_library": PREREGISTERED_LIBRARY,
            "matches_preregistration": matches,
            "quantiles": {"heavy_atoms": list(HA_QUANTILES),
                          "rotatable_bonds": list(RB_QUANTILES)},
            "bin_rule": "bin = number of edges <= value, i.e. the half-open interval "
                        "[edge[k-1], edge[k])",
            "rotatable_bond_definition":
                "ligands stage (ligsize.py): single, non-ring (bridge), non-terminal "
                "heavy-atom bonds, amides included; NOT RDKit/Lipinski",
            "control_set": CONTROL_SET,
            "n_library": len(library)}


def _m_table(panel, library, strata, min_coverage, m_raw, m_reps, sets) -> list:
    """One row per library ligand: its M, its M per replicate, its size and its stratum.

    `m_raw` is {ID: (M, n_proteins, n_targets)} over every scored ligand; `m_reps` is
    {replicate: {ID: M}} after the coverage filter."""
    rows = []
    for r in library:
        lid = r["ID"]
        hb, rbb = strata["bins"][lid]
        got = m_raw.get(lid)
        member = sorted(g for g, members in sets.items() if lid in members)
        sub = sorted({members[lid]["subgroup"] for g, members in sets.items()
                      if lid in members and members[lid]["subgroup"]})
        row = {
            "ID": lid, "name": r["name"], "heavy_atoms": r["heavy_atoms"],
            "rotatable_bonds": r["rotatable_bonds"], "stratum": strata["of"][lid],
            "ha_bin": hb, "rb_bin": rbb,
            "stratum_fine": strata["fine_of"][lid],
            "n_proteins": got[1] if got else 0, "n_targets": got[2] if got else 0,
            "coverage": round(got[1] / panel.n_proteins, 6)
            if got and panel.n_proteins else 0,
            "target_coverage": round(got[2] / panel.n_targets, 6)
            if got and panel.n_targets else 0,
            "usable": bool(got) and got[1] >= min_coverage * panel.n_proteins
            and got[2] >= min_coverage * panel.n_targets,
            "M": _num(got[0]) if got else "",
            "sets": ";".join(member), "subgroups": ";".join(sub)}
        for i in sorted(m_reps):
            row[f"M_r{i}"] = _num(m_reps[i].get(lid, ""))
        rows.append(row)
    return rows


def _usable_m(panel, m_raw, min_coverage, allowed=None) -> dict:
    """{ID: M} for ligands scored on at least `min_coverage` of the panel's proteins AND
    of the panel's TARGETS.

    The target-level condition is not redundant. The gate counts proteins, but M is
    acutely sensitive to how many pockets of a protein were scored: `E[max of k uniforms]
    = k/(k+1)`, so losing one pocket of a two-pocket protein costs 1/6 of a percentile,
    and 78 % of the HCC panel's proteins have two pockets. A ligand that loses the SECOND
    pocket of all 211 two-pocket proteins still has protein coverage 1.000 while its M
    moves -7.6 sd; a realistic 5 % per-(ligand, target) failure rate applied to one group
    moves it -0.77 sd with minimum protein coverage still 0.970 -- about the size of the
    effect this test calls significant. Protein coverage alone cannot see either.

    `allowed` (the library's IDs) keeps anything the strata do not cover out of the
    controls: a ligand with no stratum cannot be matched on size, and an unmatched control
    is exactly what this design exists to avoid."""
    need_p = min_coverage * panel.n_proteins
    need_t = min_coverage * panel.n_targets
    return {lid: v[0] for lid, v in m_raw.items()
            if v[1] >= need_p and v[2] >= need_t
            and (allowed is None or lid in allowed)}


def _panel_of(run_dir, protein_of, pooled, replicates) -> Panel:
    """One Panel from the loaded tables, taking only this run folder's rows."""
    run_dir = Path(run_dir)
    key = str(run_dir)
    mine = {label: v for (k, label), v in pooled.items()
            if k == key and label in protein_of}
    reps = {}
    for (k, rep), by_label in replicates.items():
        if k != key:
            continue
        reps[rep] = {label: v for label, v in by_label.items() if label in protein_of}
    return Panel(run_dir.name or "panel", run_dir, protein_of, mine, reps)


def run(layout, args, progress) -> dict:
    inputs = _inputs(layout, args)
    (protein_of, jobs, meta), control = _all_jobs(layout, args)
    all_meta = dict(meta)
    all_jobs = list(jobs)
    if control is not None:
        all_jobs += control[2]
        all_meta.update(control[3])
    done = _finished(layout, args, all_meta)
    if done is not None:
        return done

    d = stage_dir(layout)
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").unlink(missing_ok=True)   # the outputs below are about to change
    (d / "failed.tsv").write_text("item\tstage\treason\n")
    logged: list = []

    # 1. the strata, from the library alone, written BEFORE any outcome is read.
    library = read_library(inputs["library"])
    strata = build_strata(library, args)
    _atomic_text(d / STRATA, json.dumps(
        {k: v for k, v in strata.items()
         if k not in ("of", "bins", "order", "fine_of")},
        indent=2) + "\n")
    if not strata["matches_preregistration"]:
        logged.append((STRATA, f"the strata edges computed from this library "
                               f"({strata['edges']}) differ from the pre-registered ones "
                               f"({PREREGISTERED_EDGES}, from {PREREGISTERED_LIBRARY}): "
                               f"the library changed, and the pre-registration must be "
                               f"amended explicitly, not silently"))

    # 2. the positive set.
    sets, unresolved = read_positive_sets(inputs["positive_set"], library)
    logged.extend(unresolved)

    # 3. the outcomes.
    failures: list = []
    dropped: list = []
    pooled, replicates = _load(all_jobs, all_meta, args, progress, failures, dropped)
    logged.extend(failures)
    logged.extend((item, why) for item, why, _n in dropped)
    dropped_rows = sum(n for _i, _w, n in dropped)

    panels = [_panel_of(layout.run_dir, protein_of, pooled, replicates)]
    if control is not None:
        cdir, cprotein, _cjobs, _cmeta = control
        panels.append(_panel_of(cdir, cprotein, pooled, replicates))
        if panels[1].name == panels[0].name:   # two run folders with the same basename
            panels[1].name = f"{panels[1].name}_control"

    for panel in panels:
        # a panel with nothing readable must not come out as a table of empty p-values:
        # every read failing looks exactly like "no effect" in the output otherwise.
        # `not panel.pooled` is not enough: a file that exists but whose every row is
        # unparseable gives {label: {}}, which is non-empty and sails through to a table of
        # blank statistics and blank p-values that reads as "no effect".
        if not any(panel.pooled.values()):
            _flush(d, logged)      # the read failures are the explanation; keep them
            raise RuntimeError(
                f"{panel.run_dir}: not one usable score was read from the "
                f"df_results_<label>.csv files of its {len(panel.protein_of)} targets "
                f"({len(panel.pooled)} file(s) read, {len(failures)} read failure(s), "
                f"{dropped_rows} unusable row(s); see {d / 'failed.tsv'}). "
                f"Has `polyphar dock` finished there?")

    verdict, control_detail = positive_control_verdict(
        inputs["positive_control"],
        display=getattr(args, "positive_control", POSITIVE_CONTROL))

    if verdict == "yes" and "source" not in control_detail:
        logged.append((str(inputs["positive_control"]),
                       "the positive control passes on its numbers but states no "
                       "`source`: nothing records how the self-rank, redocking-RMSD and "
                       "library-percentile figures were produced"))

    min_coverage = float(getattr(args, "min_coverage", MIN_COVERAGE))
    sizes = {r["ID"]: r["heavy_atoms"] for r in library}
    result_rows, drug_rows, surface, balance = [], [], [], []
    per_panel = {}
    for panel in panels:
        pct = percentile_matrix(panel.pooled)
        m_raw = m_statistic(pct, panel.protein_of)
        in_library = set(strata["of"])
        m_pooled = _usable_m(panel, m_raw, min_coverage, in_library)
        m_reps = {}
        for i, by_label in sorted(panel.replicates.items()):
            rp = percentile_matrix(by_label)
            m_reps[i] = _usable_m(panel, m_statistic(rp, panel.protein_of), min_coverage,
                                  in_library)
        per_panel[panel.name] = (panel, m_raw, m_pooled, m_reps)
        _write_csv(d / f"m_statistic_{_safe(panel.name)}.csv",
                   M_COLUMNS + [f"M_r{i}" for i in sorted(m_reps)],
                   _m_table(panel, library, strata, min_coverage, m_raw, m_reps, sets))

    primary = panels[0]
    # positives that resolved to a dockable record but did not clear the coverage rule.
    # Pre-registration section 9 promises every guideline drug dropped is NAMED.
    n_targets_of = {lid: v[2] for lid, v in per_panel[primary.name][1].items()}
    below = []
    for drug_set in sorted(sets):
        for lid in sorted(sets[drug_set]):
            if lid in per_panel[primary.name][2]:
                continue
            got = per_panel[primary.name][1].get(lid)
            where = f"{drug_set}:{lid} {sets[drug_set][lid]['name']}"
            if got is None:
                logged.append((where, f"dropped from the test: not scored on any target "
                                      f"of {primary.name}"))
            else:
                logged.append((where, f"dropped from the test: scored on {got[1]} of "
                                      f"{primary.n_proteins} proteins and {got[2]} of "
                                      f"{primary.n_targets} targets, below the "
                                      f"pre-registered --min-coverage {min_coverage:g}"))
            below.append(where)

    for panel in panels:
        _panel, _raw, m_pooled, m_reps = per_panel[panel.name]
        is_primary_panel = panel is primary
        ntg = {lid: v[2] for lid, v in _raw.items()}
        for drug_set in sorted(sets):
            members = sets[drug_set]
            family = (FAMILY_PRIMARY if drug_set == SET_HCC else FAMILY_EXPLORATORY)
            note = ""
            if drug_set != SET_HCC:
                note = (f"exploratory: n = {len(members)} in the positive set, of which "
                        f"n_positive on this row entered the test; no correction is "
                        f"applied and this arm makes no primary claim")
            result_rows.append(_arm_rows(panel, members, drug_set, family, m_pooled,
                                         m_reps, strata, args,
                                         verdict if family != FAMILY_EXPLORATORY
                                         else f"exploratory ({verdict})", note,
                                         sizes=sizes, targets=ntg, balance_out=balance))
            surface.extend(_surface_for(panel, members, m_pooled, drug_set))
            if is_primary_panel:
                drug_rows.extend(_per_drug(panel, members, m_pooled, strata, drug_set,
                                           library))
        if is_primary_panel and SET_HCC in sets:
            for sub in SUBGROUPS:
                members = {lid: v for lid, v in sets[SET_HCC].items()
                           if v["subgroup"] == sub}
                if not members:
                    continue
                result_rows.append(_arm_rows(
                    panel, members, f"{SET_HCC}_{sub}", FAMILY_SUBGROUP, m_pooled, m_reps,
                    strata, args, verdict,
                    "pre-planned subgroup: the HCC set is bimodal and the halves carry "
                    "opposite size biases", sizes=sizes, targets=ntg,
                    balance_out=balance))
        # the pre-registered SENSITIVITY analysis: the same arm on finer heavy-atom bins.
        if is_primary_panel and SET_HCC in sets:
            result_rows.append(_arm_rows(
                panel, sets[SET_HCC], SET_HCC, FAMILY_SENSITIVITY, m_pooled, m_reps,
                strata, args, verdict,
                "pre-registered SENSITIVITY analysis, not the primary test: heavy-atom "
                "deciles x rotatable-bond tertiles instead of quintiles x tertiles. The "
                "primary HA bins are wide (bin 4 is open-ended, 32-63 heavy atoms) and "
                "stratification removes the between-stratum size confound, not the "
                "within-stratum one. A primary result that does not survive here is a "
                "result that may be molecular size",
                sizes=sizes, targets=ntg, stratum_key="fine_of",
                arm_name=f"{SET_HCC} on {panel.name} (fine HA bins)",
                balance_out=balance))

    # The pre-registered prediction's OWN test (P1): HCC-on-HCC > HCC-on-AD, as a paired
    # within-stratum permutation on M_primary - M_control. Without it the prediction is
    # only ever read off two separate verdicts, which is the difference-of-significance
    # error -- and the pre-registration calls this the prediction that matters.
    if control is not None and SET_HCC in sets:
        m_a = per_panel[panels[0].name][2]
        m_b = per_panel[panels[1].name][2]
        diff = paired_difference(m_a, m_b)
        rep_diff = {}
        for i in sorted(set(per_panel[panels[0].name][3]) & set(per_panel[panels[1].name][3])):
            rep_diff[i] = paired_difference(per_panel[panels[0].name][3][i],
                                            per_panel[panels[1].name][3][i])
        result_rows.append(_arm_rows(
            panels[0], sets[SET_HCC], SET_HCC, FAMILY_CONTRAST, diff, rep_diff, strata,
            args, verdict,
            f"THE PRE-REGISTERED PREDICTION, tested directly: HCC-on-HCC > HCC-on-AD, as "
            f"a paired within-stratum permutation on M({panels[0].name}) - "
            f"M({panels[1].name}) over the {len(diff)} ligands carrying an M on both "
            f"panels. Concluding this prediction from 'one arm is significant and the "
            f"other is not' is the difference-of-significance error; this row is the "
            f"contrast itself. It carries no Holm correction: it is one test of one "
            f"pre-registered claim",
            sizes=sizes, targets=n_targets_of,
            arm_name=f"{SET_HCC}: {panels[0].name} minus {panels[1].name}",
            direction=DIRECTION_CONTRAST, balance_out=balance))

    if control is None:
        for row in result_rows:
            if row["family"] == FAMILY_PRIMARY:
                row["note"] = ((row["note"] + "; ") if row["note"] else "") + (
                    "the second primary test (the same drugs on another panel, the "
                    "pre-registered negative control) was NOT run: pass --control-panel "
                    "<run folder of the other disease>")

    _adjust(result_rows)
    # Holm is correct for the family of two, but it cuts the WRONG way on the arm that is
    # supposed to fail: a Holm-inflated non-significant negative control is weaker evidence
    # than its raw p. Say so on the row, not only in prose.
    if control is not None:
        for row in result_rows:
            if row["family"] == FAMILY_PRIMARY and row["panel"] != primary.name:
                row["note"] = ((row["note"] + "; ") if row["note"] else "") + (
                    f"negative control: this arm is the one we WANT to fail to reject, "
                    f"and Holm makes it harder to reject, which is backwards for such a "
                    f"test. Read the RAW p ({row['p_value']}) as the honest figure here; "
                    f"the adjusted p ({row['p_adjusted']}) is reported for completeness "
                    f"and is conservative in the unhelpful direction. The contrast row "
                    f"({FAMILY_CONTRAST} family) is what actually tests the prediction")

    _write_csv(d / RESULTS, RESULT_COLUMNS, result_rows)
    _write_csv(d / PER_DRUG, PER_DRUG_COLUMNS, drug_rows)
    _write_csv(d / BALANCE, BALANCE_COLUMNS, balance)
    write_no_pvalues(d / SURFACE, SURFACE_COLUMNS, surface)
    _atomic_text(d / SURFACE_MD, surface_markdown(surface))
    _atomic_text(d / SUMMARY_MD, _summary(panels, sets, strata, result_rows, verdict,
                                          unresolved, args, below, balance, min_coverage))

    _flush(d, logged)

    counts = {
        "panels": len(panels), "targets": primary.n_targets,
        "proteins": primary.n_proteins, "library": len(library),
        "positives": sum(len(v) for v in sets.values()),
        "positives_unresolved": len(unresolved),
        "positives_below_coverage": len(below),
        "arms": len(result_rows),
        "primary_tests": sum(1 for r in result_rows if r["family"] == FAMILY_PRIMARY),
        "contrast_tests": sum(1 for r in result_rows if r["family"] == FAMILY_CONTRAST),
        "balance_flagged": sum(1 for r in result_rows if r["balance_flag"]),
        "permutations": int(getattr(args, "permutations", N_PERMUTATIONS)),
        "replicates": len(primary.replicates),
        "strata": len(strata["sizes"]),
        "strata_edges_preregistered": bool(strata["matches_preregistration"]),
        "positive_control": verdict.split(":")[0],
        "failed": len(failures),
        "dropped_score_rows": dropped_rows,
    }
    write_manifest(d, STAGE,
                   args={"nprocs": getattr(args, "nprocs", 1), **_params(args)},
                   inputs=[p for p in inputs.values() if Path(p).is_file()],
                   versions={"statistic_version": 2,
                             "preregistration": PREREGISTRATION,
                             # the DOCUMENT, not a path string: a path records nothing, and
                             # an amendment made after the numbers exist left no trace.
                             "preregistration_sha256": preregistration_sha256(),
                             "preregistration_sha256_expected": PREREGISTRATION_SHA256,
                             "control_set": CONTROL_SET,
                             **{f"preregistered_{k}": v
                                for k, v in PREREGISTERED_EDGES.items()}},
                   counts=counts, timing={},
                   extra={"signature": _signature(layout, args, all_meta),
                          "positive_control": control_detail,
                          "outputs": {name: (d / name).stat().st_size
                                      for name in (STRATA, RESULTS, PER_DRUG, BALANCE,
                                                   SURFACE, SURFACE_MD, SUMMARY_MD)}})
    return counts


def _flush(stage_dir, logged) -> None:
    """Write everything noticed so far to failed.tsv (also on the way out of a failure:
    the read failures are what explains it)."""
    for item, reason in logged:
        append_failure(stage_dir, item, STAGE, reason)
    logged.clear()


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)) or "panel"


def _surface_for(panel, members, m_pooled, drug_set) -> list:
    """Per-ligand score lists for the descriptive grid, positives against the controls."""
    scores: dict = {}
    for label, by_id in panel.pooled.items():
        for lid, s in by_id.items():
            scores.setdefault(lid, []).append(s)
    pos = {lid: v for lid, v in scores.items() if lid in members and lid in m_pooled}
    ctl = {lid: v for lid, v in scores.items() if lid not in members and lid in m_pooled}
    return surface_rows(panel.name, drug_set, pos, ctl)


def _per_drug(panel, members, m_pooled, strata, drug_set, library) -> list:
    """Exploratory per-drug elevation against the controls of the drug's own stratum,
    with Benjamini-Hochberg FDR. Labelled exploratory on every row."""
    by_stratum: dict = {}
    for lid, m in m_pooled.items():
        if lid in members:
            continue
        by_stratum.setdefault(strata["of"][lid], []).append(m)
    size = {r["ID"]: r for r in library}
    rows = []
    for lid in sorted(members):
        if lid not in m_pooled:
            continue
        s = strata["of"][lid]
        controls = by_stratum.get(s, [])
        rows.append({"panel": panel.name, "drug_set": drug_set, "ID": lid,
                     "name": members[lid]["name"],
                     "heavy_atoms": size.get(lid, {}).get("heavy_atoms", ""),
                     "rotatable_bonds": size.get(lid, {}).get("rotatable_bonds", ""),
                     "stratum": s,
                     "M": _num(m_pooled[lid]), "n_control_in_stratum": len(controls),
                     "p_value": _num(single_drug_p(m_pooled[lid], controls)),
                     "q_value_bh": "", "status": "exploratory (BH FDR)"})
    # BH over the RESOLVABLE subset: one drug sitting in a stratum with no control (p = "")
    # used to wipe out every q-value in the table.
    usable = [r for r in rows if r["p_value"] != ""]
    if usable:
        for r, q in zip(usable, benjamini_hochberg([r["p_value"] for r in usable])):
            r["q_value_bh"] = round(q, 6)
    for r in rows:
        if r["p_value"] == "":
            r["status"] = (f"exploratory (no control in stratum {r['stratum']}: no p, and "
                           f"this drug is not in the BH family of {len(usable)})")
    return rows


def _adjust(rows) -> None:
    """Holm within each correction family; exploratory arms are left uncorrected.

    The primary family is exactly the two pre-registered tests (HCC drugs on this panel,
    HCC drugs on the control panel); the subgroup family is the targeted / cytotoxic
    split. Nothing else is corrected, and nothing else is in a family."""
    for family in (FAMILY_PRIMARY, FAMILY_SUBGROUP):
        members = [r for r in rows if r["family"] == family and r["p_value"] != ""]
        if not members:
            continue
        for r, p in zip(members, holm([r["p_value"] for r in members])):
            r["p_adjusted"] = round(p, 6)
            r["adjustment"] = f"holm({family}, m={len(members)})"
    for r in rows:
        if r["family"] == FAMILY_EXPLORATORY:
            r["adjustment"] = "none (exploratory; n is on this row)"
        elif r["family"] == FAMILY_CONTRAST:
            r["adjustment"] = ("none (one test of the one pre-registered prediction; "
                               "it is not in the Holm family of two)")
        elif r["family"] == FAMILY_SENSITIVITY:
            r["adjustment"] = "none (pre-registered sensitivity analysis, not a test)"


def _summary(panels, sets, strata, result_rows, verdict, unresolved, args,
             below=(), balance=(), min_coverage=MIN_COVERAGE) -> str:
    out = ["# `polyphar stats`: the cumulative weak-binding test", "",
           f"Pre-registered in `{PREREGISTRATION}` (sha256 "
           f"`{preregistration_sha256()[:16] or 'not on disk'}`); "
           "designed in `docs/superpowers/specs/2026-11-03-statistics-design.md`.", "",
           "## What the controls are", "",
           CONTROL_SET + ".", "",
           "## Positive control (pre-registered gate)", "",
           f"**{verdict}**", "",
           f"Pass thresholds: median self-recovery rank within the top "
           f"{SELF_RANK_MAX:g} of the panel over at least {MIN_SELF_RANK_N} ligands; "
           f"redocking success (RMSD <= {REDOCK_RMSD_MAX:g} A) >= "
           f"{REDOCK_SUCCESS_MIN:g} over at least {MIN_REDOCK_N} self pairs; at least "
           f"{LIBRARY_TOP_DECILE_MIN:g} of the co-crystal ligands in the top decile of "
           f"the production library at their own target, over at least "
           f"{MIN_LIBRARY_N} ligands; and at least "
           f"{MIN_PANEL_SCORED_FRACTION:g} of the panel's targets scored. "
           f"Pose-in-pocket is reported descriptively and is NOT a gate: its null on "
           f"random non-clashing centroids in the engine's own box is about 0.79. "
           f"If a gate fails, the primary result is reported as uninterpretable rather "
           f"than reinterpreted.", "",
           "## Panels", ""]
    for panel in panels:
        out.append(f"* **{panel.name}** ({panel.run_dir}): {panel.n_targets} targets / "
                   f"{panel.n_proteins} proteins, {len(panel.replicates)} replicate(s).")
    if len(panels) > 1:
        a, b = panels[0].accessions, panels[1].accessions
        shared = sorted(a & b)
        out += ["", f"Accession-level overlap between the panels: {len(shared)} of "
                    f"{len(a)} and {len(b)} ({', '.join(shared[:20])}"
                    f"{', ...' if len(shared) > 20 else ''}). The panels are not "
                    f"independent samples of protein space; shared targets compress the "
                    f"HCC-on-HCC minus HCC-on-AD difference."]
    else:
        out += ["", "**The negative-control panel was not supplied** (`--control-panel`), "
                    "so only one of the two pre-registered primary tests was run."]
    out += ["", "## Strata", "",
            f"Heavy-atom quintile edges {strata['edges']['heavy_atoms']}, "
            f"rotatable-bond tertile edges {strata['edges']['rotatable_bonds']} "
            f"({'match' if strata['matches_preregistration'] else 'DO NOT MATCH'} the "
            f"pre-registered edges). Bin sizes: "
            + ", ".join(f"{k}={v}" for k, v in sorted(strata["sizes"].items())) + ".",
            "", "## Positive set", ""]
    for name in sorted(sets):
        subs = {}
        for v in sets[name].values():
            subs[v["subgroup"] or "-"] = subs.get(v["subgroup"] or "-", 0) + 1
        out.append(f"* `{name}`: n = {len(sets[name])} dockable "
                   f"({', '.join(f'{k}: {v}' for k, v in sorted(subs.items()))}).")
    if unresolved:
        out += ["", "Rows that did not resolve to a dockable library record (named here "
                    "because the Methods must say which guideline drug was dropped):", ""]
        out += [f"* {item} - {why}" for item, why in unresolved]
    out += ["", "## Coverage", "",
            f"A ligand carries an `M` only if it was scored on at least "
            f"{min_coverage:g} of the primary panel's proteins **and** of its targets. "
            f"The target-level half of that rule is not redundant: `E[max of k uniforms] "
            f"= k/(k+1)`, so losing one pocket of a two-pocket protein costs 1/6 of a "
            f"percentile while protein coverage stays at 1.000, and most of the HCC "
            f"panel's proteins have two pockets. `n_targets` per positive and per control "
            f"is in `balance.csv`; per-ligand coverage is in `m_statistic_<panel>.csv`."]
    if below:
        out += ["", f"**{len(below)} positive(s) resolved to a dockable record but did "
                    f"not clear the coverage rule and are NOT in any arm** (also in "
                    f"`failed.tsv`):", ""]
        out += [f"* {item}" for item in below]
    else:
        out += ["", "Every resolved positive cleared the coverage rule."]

    out += ["", "## Arms", "",
            "| arm | family | n+ | n- | statistic | p | p adjusted | "
            "size balance (sd) | interpretable | replicates |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for r in result_rows:
        reps = ", ".join(str(r.get(f"stat_r{i}", "")) for i in range(3)
                         if r.get(f"stat_r{i}", "") != "")
        # the verdict must travel WITH the number: this table is what gets pasted into a
        # Methods draft, and it used to be byte-identical for pass, fail and absent.
        mark = (r["family"] in (FAMILY_PRIMARY, FAMILY_SUBGROUP, FAMILY_CONTRAST)
                and not str(r["interpretable"]).startswith("yes"))
        p = f'{r["p_value"]} (UNINTERPRETABLE)' if mark and r["p_value"] != "" \
            else r["p_value"]
        padj = f'{r["p_adjusted"]} (UNINTERPRETABLE)' if mark and r["p_adjusted"] != "" \
            else r["p_adjusted"]
        bal = r["balance_sd"]
        if r["balance_flag"]:
            bal = f"{bal} **FLAGGED**"
        out.append(f'| {r["arm"]} | {r["family"]} | {r["n_positive"]} | {r["n_control"]} '
                   f'| {r["statistic"]} | {p} | {padj} | {bal} '
                   f'| {str(r["interpretable"])[:40]} | {reps} |')

    out += ["", f"Permutations: {getattr(args, 'permutations', N_PERMUTATIONS)}, "
                f"seed {getattr(args, 'permutation_seed', SEED)}, one-sided ({DIRECTION}).",
            ""]

    # --- the pre-registered prediction, tested -------------------------------------------
    contrast = [r for r in result_rows if r["family"] == FAMILY_CONTRAST]
    prim = [r for r in result_rows if r["family"] == FAMILY_PRIMARY]
    out += ["## The pre-registered prediction: HCC-on-HCC > HCC-on-AD", ""]
    if contrast:
        c = contrast[0]
        out += [f"Tested directly, as a **paired within-stratum permutation on "
                f"M({panels[0].name}) - M({panels[1].name})** over the ligands carrying an "
                f"`M` on both panels: statistic **{c['statistic']}**, "
                f"p **{c['p_value']}** (one-sided, uncorrected).", "",
                "This row exists because the prediction cannot be read off the two primary "
                "arms. Each primary arm is tested against **its own** matched controls; "
                "concluding the prediction holds because one arm is significant and the "
                "other is not is the **difference-of-significance error**. The same "
                "ligand carries an `M` on both panels, so the difference is paired within "
                "ligand and every ligand-level nuisance cancels in it.", ""]
    if len(prim) == 2:
        a = next((r for r in prim if r["panel"] == panels[0].name), None)
        b = next((r for r in prim if r["panel"] != panels[0].name), None)
        if a and b and a["statistic"] != "" and b["statistic"] != "":
            bigger = a if float(a["statistic"]) >= float(b["statistic"]) else b
            out += [f"Descriptively, the larger arm statistic is **{bigger['arm']}** "
                    f"({a['arm']} {a['statistic']} vs {b['arm']} {b['statistic']}).", ""]
        if b:
            out += [f"**Which p to read on the negative control.** Holm over the family of "
                    f"two is correct, but it makes the negative control **harder to "
                    f"reject**, and this is the arm we want to fail to reject: a "
                    f"Holm-inflated non-significant control is weaker evidence than its "
                    f"raw p, not stronger. Raw p = **{b['p_value']}** is the honest figure "
                    f"to read for `{b['arm']}`; the adjusted p = {b['p_adjusted']} is "
                    f"reported for completeness and is conservative in the unhelpful "
                    f"direction.", ""]
    elif not contrast:
        out += ["**Not tested.** The control panel was not supplied, so there is no "
                "contrast and no negative-control arm. One significant arm on its own is "
                "not the pre-registered primary analysis.", ""]

    # --- the limitation the design does NOT remove ---------------------------------------
    flagged = [r for r in result_rows if r["balance_flag"]]
    out += ["## Limitation: stratification removes the BETWEEN-stratum size confound, "
            "not the WITHIN-stratum one", "",
            f"Matching on heavy-atom quintile x rotatable-bond tertile makes the "
            f"comparison size-matched **across** strata. It does not make the positives "
            f"and their controls the same size **inside** a stratum, and the primary bins "
            f"are wide: HA bin 4 is open-ended and spans 32-63 heavy atoms. With scores "
            f"that are a pure monotone function of heavy-atom count and carry no binding "
            f"signal at all, positives drawn uniformly inside the big strata come out "
            f"correctly null, but a within-stratum heavy-atom imbalance of about "
            f"**+0.6 sd already reaches P(p <= 0.05) = 0.61**, +0.8 sd reaches 0.99 and "
            f">= +1.15 sd reaches 1.00. **A within-stratum size imbalance of ~0.6 sd or "
            f"more can produce p < 0.05 from molecular size alone.**", "",
            f"`balance.csv` reports, per arm and per stratum, the positives' mean and "
            f"median heavy atoms against their matched controls' and the standardised "
            f"difference; the arm rows carry the weighted figure as `balance_sd` and a "
            f"`balance_flag` at |{BALANCE_FLAG_SD:g}| sd. The pre-registered sensitivity "
            f"analysis (the `{FAMILY_SENSITIVITY}` arm, heavy-atom **deciles** x "
            f"rotatable-bond tertiles) is the second reading: a primary result that does "
            f"not survive finer size bins is a result that may be molecular size.", ""]
    if flagged:
        out += ["**Flagged arms:**", ""]
        out += [f"* {r['arm']}: {r['balance_flag']}" for r in flagged]
        out += [""]
    else:
        out += [f"No arm is flagged: every arm's weighted within-stratum heavy-atom "
                f"imbalance is below {BALANCE_FLAG_SD:g} sd.", ""]

    out += ["## How to read the rest", "",
            "* `surface.csv` / `surface.md` are DESCRIPTIVE and carry no p-values.",
            "* `per_drug.csv` is exploratory with Benjamini-Hochberg FDR over the "
            "resolvable subset (a drug in a stratum with no control carries no p and is "
            "not in the BH family).",
            "* `replicate_sign_jitter_check` is **not corroboration**. corr(M_r0, M_r1) is "
            "about 0.9996, so it is True on essentially every arm, including arms with "
            "p = 0.45 and p = 0.80. It can only catch a catastrophic jitter regime in "
            "which replicates disagree about the sign at all; it says nothing about "
            "whether an effect is real.",
            "* An arm whose `interpretable` verdict is not `yes` has its p-values marked "
            "UNINTERPRETABLE in the table above. They are still printed, because hiding a "
            "number invites recomputing it elsewhere; the label is the point.", ""]
    return "\n".join(out)
