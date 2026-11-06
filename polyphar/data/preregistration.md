# Pre-registration: the cumulative weak-binding test

**Date:** 2026-11-04
**Status:** pre-registration. Written **before `dock` has run** and therefore before any
value of the statistic exists, and committed in the same commit as the `polyphar stats`
stage that pins its sha256 as `stats.PREREGISTRATION_SHA256`. That pin, not this sentence,
is what makes the claim checkable: `polyphar stats` records the document's sha256 in every
`08_stats/manifest.json`, and `tests/test_stage_stats.py` fails until an edit to this file
is accompanied by a new constant and an entry in **§ Amendments** below. The file's mtime
proves nothing (it lives on a Windows drive); the git history and the hash do.
**Design:** `docs/superpowers/specs/2026-11-03-statistics-design.md`
**Implementation:** `polyphar stats` (`polyphar/stages/stats.py`, outputs in `08_stats/`)

This document fixes the analysis. Its only value is that it exists before the numbers do.
Anything not written here is exploratory and must be labelled as such wherever it is
reported. Any change to this document after a result is known must be made as a visible
amendment, dated, saying what changed and why — not as an edit.

---

## 1. The claim

> Approved HCC drugs bind more of the HCC panel, at a useful level, than comparable
> molecules do.

"Comparable" carries the whole argument: approved drugs are larger than the library median
and Vina-family scoring is roughly additive over atoms, so an unmatched comparison comes
out positive whether or not the claim is true. The design below exists so that **ligand
size and pocket generosity cannot produce the result**.

## 2. The statistic

**Step 1 — within-target percentile.** For each target (one pocket of one protein), every
ligand scored on it is ranked by `score_min` from `07_docking/df_results_<label>.csv` — the
best of that ligand's replicate draws at that target (see the 2026-11-06 amendment).

A ligand with fewer draws than its target's best-attested count is **dropped** from that
target, counted and named: min-of-k improves monotonically with k, so ranking a 3-draw
ligand against a 2-draw one would reward the extra draw rather than the molecule. The
dropped ligand is then subject to the coverage rule below like any other gap.

* Direction: **more negative = better binding = percentile closer to 1**.
* Ties: **mid-rank**.
* `percentile = 1 - (rank - 0.5) / N`, N = ligands scored on that target. The percentiles
  of every target therefore have mean exactly 0.5.

This **removes** pocket generosity rather than adjusting for it. No descriptors, no fitted
model, no calibration run.

**Step 2 — the statistic.** For ligand *i*:

```
M_i = mean over PROTEINS of ( max over that protein's targets of percentile_t(i) )
```

Per protein, best pocket — **not** a mean over all targets, which would weight a protein
contributing two pockets twice. The HCC panel is 271 proteins / 482 targets.

**Coverage rule.** A ligand enters the test only if it is scored on at least **98 %** of
the panel's **proteins** *and* at least **98 %** of the panel's **targets**
(`--min-coverage 0.98`). Ligands below that are excluded, counted
(`counts.positives_below_coverage`) and, when they are positives, **named** in
`failed.tsv` and in `summary.md`. No ligand is given an imputed score.

Both halves of that rule are load-bearing.

* **Why 0.98 and not 0.9.** Under informative missingness — a ligand failing to dock where
  it fits worst, which is the realistic direction — `M` is biased upward by up to
  `0.5 x (1 - coverage)`: **+0.011 at 0.98, +0.026 at 0.95, +0.051 at 0.90**. The effect
  this test is powered to detect is of that order, so a 0.9 gate admits a bias the size of
  the result.
* **Why targets and not only proteins.** `E[max of k uniforms] = k/(k+1)`, so losing one
  pocket of a two-pocket protein costs 1/6 of a percentile, and about 78 % of the HCC
  panel's proteins contribute two pockets. A ligand that loses the *second* pocket of all
  211 two-pocket proteins keeps protein coverage **1.000** while its `M` moves **−7.6 sd**;
  a 5 % per-(ligand, target) failure rate applied to one group moves it −0.77 sd with
  minimum protein coverage still 0.970. A protein-level gate cannot see either. The
  per-arm distribution of `n_targets` for positives and controls is reported in
  `balance.csv`.

## 3. Strata (the matching scheme)

**Quintile bins of heavy-atom count crossed with tertile bins of rotatable-bond count**
(15 cells). Nearest-neighbour matching was considered and rejected: strata are easier to
describe and harder to tune after the fact.

* The rotatable-bond count is the `ligands` stage's **pinned, non-standard** definition —
  single, non-ring (decided by bridge detection on the heavy-atom bond graph),
  non-terminal heavy-atom bonds, **amide C–N NOT excluded**. RDKit's / Open Babel's
  conventional Lipinski count differs for **1,316 of the 2,083** e-Drug3D records and must
  not be substituted.
* A bin is the half-open interval `[edge[k-1], edge[k])`, i.e. `bin = #{e in edges : e <= value}`.

### The bin edges, computed now

Computed from **e-Drug3D_2083**, over the **1,960 records that survive the dockability
filter** (rotatable bonds <= 15 and MW <= 900), on 2026-11-04, before any docking score
for the new panels existed:

| axis | quantiles | **edges** |
|---|---|---|
| heavy atoms | 0.20, 0.40, 0.60, 0.80 | **17, 22, 26, 32** |
| rotatable bonds | 1/3, 2/3 | **4, 7** |

Realised bin sizes (the counts are unequal because both variables are integers; this is
expected and is recorded, not corrected):

| | rotB bin 0 (<4) | 1 (4–6) | 2 (>=7) | total |
|---|---|---|---|---|
| **HA bin 0 (<17)** | 233 | 93 | 15 | 341 |
| **1 (17–21)** | 186 | 175 | 59 | 420 |
| **2 (22–25)** | 108 | 158 | 97 | 363 |
| **3 (26–31)** | 57 | 170 | 190 | 417 |
| **4 (>=32)** | 23 | 97 | 299 | 419 |
| **total** | 607 | 693 | 660 | **1,960** |

All 15 cells are populated. `polyphar stats` recomputes the edges from the run folder's own
`07_docking/library/ligands.csv`, writes them to `08_stats/strata.json` **before any
docking result is read**, and compares them against the values above; a mismatch is
recorded (`counts.strata_edges_preregistered = false`) and logged, never silently adopted.

### What the controls are

**The library is e-Drug3D**: approved and withdrawn small-molecule drugs, not a generic
screening set. The size-matched comparison set is therefore **already "other approved
drugs of the same size"**. This is the contrast a reviewer asks for first, it is
load-bearing, and it is easy to miss, so it is stated in `results.csv`, `surface.md` and
`summary.md` of every run.

### 3b. What stratification does NOT remove, and the pre-registered sensitivity analysis

Stratification removes the **between-stratum** size confound. It does **not** remove the
**within-stratum** one, and the primary bins are wide — HA bin 4 is open-ended and spans
**32–63** heavy atoms.

Quantified on the real 1,960-record library, with a 271-protein / 482-target panel whose
scores are a **pure monotone function of heavy-atom count** and carry no binding signal at
all:

| positive set | result |
|---|---|
| drawn uniformly within the big strata | correctly null (median p 0.57–0.59) |
| the 25 largest molecules | statistic **+0.10**, **p = 0.0001** |
| within-stratum HA imbalance **+0.62 sd** | P(p ≤ 0.05) = **0.61** |
| **+0.79 sd** | **0.99** |
| **≥ +1.15 sd** | **1.00** |

The residual size-driven component of `M` is about **0.05** in every stratum (slope
`dM/dHA` 0.0073–0.0410 per atom × within-stratum sd(HA) 7.1 down to 1.0 — the two factors
trade off so the product is flat across all 15 cells), which is 2–3× the effect this test
detects at 67 % power.

> **Stated limitation: a within-stratum heavy-atom imbalance of about 0.6 sd or more can
> produce p < 0.05 from molecular size alone, with no binding signal present.**

For the plausible real HCC set (sorafenib 32, lenvatinib 30, regorafenib 33,
cabozantinib 37, apatinib 40, doxorubicin 39, sunitinib 29, erlotinib 29) the
within-stratum position is near the middle, and in HA4-RB2 *below* it (positives mean HA
35.5 against the stratum's 38.8); under the pure-size generator that set gives statistic
**−0.0041, p = 0.59** — null, as the design promises. The real analysis is therefore
probably safe, **but by accident of which drugs these are, not by construction.** Two
things are pre-registered because of that, and the test itself is **not** changed:

1. **A balance table on every arm** (`08_stats/balance.csv`): per stratum, the positives'
   mean and median heavy atoms against their matched controls', the controls' sd, the
   standardised difference, and the mean `n_targets` of each group. The arm rows of
   `results.csv` carry the weighted standardised difference as `balance_sd` and a
   `balance_flag` raised at **|0.6| sd** (`stats.BALANCE_FLAG_SD`). `summary.md` states the
   limitation and lists the flagged arms.
2. **A sensitivity analysis, pre-registered as a sensitivity analysis and not as the
   primary:** the primary HCC arm repeated on **heavy-atom deciles × rotatable-bond
   tertiles** (`stats.HA_QUANTILES_FINE`), reported as the `sensitivity` family and
   carrying no correction. A primary result that does not survive finer size bins is a
   result that may be molecular size.

## 4. The test

* **Statistic:** mean `M` over the positive set minus mean `M` over the matched controls,
  formed stratum by stratum and combined with the positives' own stratum weights
  (`w_s = n_positive_s / n_positive`). The controls are therefore averaged at the
  positives' size composition. Strata holding no positive or no control contribute
  nothing.
* **Null:** permute the drug / non-drug label **within strata**.
* **Permutations:** **10,000**.
* **Direction:** **one-sided, positives > controls.** The claim is directional.
* **p-value:** `(1 + #{permuted statistic >= observed}) / (1 + 10,000)`. The add-one
  convention: p is never reported as 0 (the design document says "proportion of
  permutations at least as extreme"; the difference is 1e-4 and only makes p larger).
* **Seed:** 42. The whole permutation runs in one process from one seeded generator, so
  `--nprocs` cannot change a result.

## 5. The arms and the multiplicity correction

| arm | status | correction |
|---|---|---|
| HCC drugs on the **HCC** panel | **primary** | Holm, family of 2 |
| HCC drugs on the **AD** panel | **primary (negative control)** | Holm, family of 2 — **but read the raw p**, see below |
| **`M` on the HCC panel − `M` on the AD panel** (paired) | **the pre-registered prediction's own test** (`contrast`) | none; one test of one claim |
| HCC targeted / HCC cytotoxic on the HCC panel | pre-planned subgroup | Holm within the subgroup family |
| HCC drugs on the HCC panel, **heavy-atom deciles** | pre-registered **sensitivity** (§ 3b) | none; not a test |
| AD drugs on either panel | **exploratory**, n about 5 | none; **n is reported wherever the arm appears** |
| per-drug elevation | **exploratory** | Benjamini–Hochberg FDR, labelled exploratory |

**Two primary tests, Holm-corrected.** The AD-drug arms carry no correction because they
make no primary claim: their set is tiny (guideline AD therapy intersected with dockable
e-Drug3D is roughly donepezil, rivastigmine, galantamine, memantine, tacrine, possibly
brexpiprazole), the drugs are known multi-target by design, and they sit in different size
strata from the HCC agents.

### Pre-registered prediction, with direction — and how it is tested

> **HCC-on-HCC > HCC-on-AD.**

Arm 2 holds the drug set fixed and changes only the panel, so nothing about the AD drugs'
own properties can contaminate it. The **accession-level overlap between the two panels**
is reported alongside the result: the panels are not independent samples of protein space,
and shared targets compress the difference.

**This prediction is tested directly, as its own statistic.** The two primary arms are each
tested against **their own** matched controls, so concluding that the prediction is met
because one arm is significant and the other is not is the **difference-of-significance
error** — and this is the prediction that matters. The test is therefore a **paired
within-stratum permutation on `M_HCC_panel − M_AD_panel`**, over the ligands carrying an
`M` on both panels, reported as the `contrast` family of `results.csv`:

* the same ligand carries an `M` on both panels, so the difference is **paired within
  ligand** and every ligand-level nuisance (size, flexibility, general stickiness) cancels
  in it;
* the **same** within-stratum label permutation distributes the contrast exactly as it
  distributes each arm's own statistic, so nothing new is assumed;
* scale is not an obstacle: `sd(M)` is essentially invariant to panel size (60 → 500
  proteins) and to the pockets-per-protein mix (0 % → 100 % two-pocket proteins, sd 0.1979
  → 0.2004) once there is shared ligand-level signal, which docking guarantees;
* it carries **no Holm correction**: it is one test of one pre-registered claim.

`summary.md` also states, descriptively, which arm's statistic is larger.

**Which p to read on the negative control.** Holm over the family of two is correct, but it
makes arm 2 — the arm we *want* to fail to reject — **harder** to reject. A Holm-inflated
non-significant negative control is **weaker** evidence than its raw p, not stronger. Both
are reported on that row and in `summary.md`, and the **raw p is named as the honest figure
to read** for that arm, with the adjusted p given for completeness.

If the negative-control panel is not supplied, only one of the two primary tests has been
run, there is no contrast, and every output says so. A single-panel result is not the
pre-registered primary analysis. Pointing `--control-panel` at the primary run folder is
refused: a panel run against itself produces two identical p-values and a negative control
that cannot fail.

## 6. Positive control: the pass thresholds, and the consequence of failing

The receptors are 63 % AlphaFold models and the engine carries an uncharacterised
+0.33 kcal/mol offset from smina. Nothing otherwise anchors a score to real binding, so a
percentile means nothing until the following readouts are produced by redocking the
co-crystallised ligands labelled from the panel's own experimental pockets.

**How many ligands to expect.** An earlier draft of this section said **687**. That is the
count of labelled pockets over **all 8,005 experimental net1 pockets**, which is not what
is redocked: only a **panel** pocket can have a self-rank, so the stage labels the
**panel's targets**. Measured label rate on panel-style pockets is **7 of 52 = 13.5 %**;
with 482 targets of which about 37 % are experimental chains (~178), the realistic
expectation is **about 24 labelled ligands**, one per physical ligand (a ligand filling a
pocket and its nested subpocket is counted once — it was being counted twice, which
inflated n without adding information). Minimum-n rules below exist because of that.

| readout | **pass threshold** | **minimum n** |
|---|---|---|
| **self-recovery rank** — each labelled ligand docked against the whole panel, the rank of its own pocket | **median self-rank within the top 5 %** of the panel | **10 ligands** |
| **redocking RMSD** — for the self pairs, the symmetry-insensitive RMSD between replicate 0's best pose and the ligand's deposited coordinates, in place, no superposition | **>= 50 % within 2.0 Å** | **10 self pairs** |
| **library percentile at home** — the co-crystal ligand's score percentile among the production library at its own target, same receptor, box and engine | **>= 50 % in the top decile** (chance 0.10) | **10 ligands** |
| **panel coverage** — panel targets that produced a score, over the panel's targets | **>= 0.98** | — |

**Why pose-in-pocket is no longer a gate.** It was pre-registered at **>= 70 %** on the
framing "the centroid is inside the cavity rather than on a surface patch at the box edge".
That framing presupposed a *padded* box. `dock.pocket_box` is the **unpadded bounding box
of the pocket's own lining atoms**, and the engine keeps the ligand inside it, so there is
no edge to park at and the criterion nearly restates a constraint the engine already
enforces. Measured on real net2 pockets, random **non-clashing** centroids confined to that
box pass **about 79 %** of the time (0.30 / 0.59 / 0.76 at insets of 2 / 4 / 6 Å) against
**about 85 %** for the true deposited poses: **the gate sat below its own null**, and an
engine placing ligands at random among non-clashing positions in its own box passed it. The
asymmetry criterion itself is sound — a centroid in van der Waals contact with real lining
passes only 17.5 % — but the engine cannot produce the geometry it rejects. It is therefore
**reported descriptively, for every (ligand, target) pair, and gates nothing.**

**Why RMSD <= 2 Å replaces it.** The stage already holds the deposited ligand's
coordinates (it is what the redocked SDF is built from), so the self pairs admit a real
**redocking-success rate** with a published base rate — Vina-family cognate redocking is
usually reported at 60–80 % — and a null of **0.0 %**: random non-clashing *rigid
placements* of the same ligand in the same box passed 0 of 1,871 draws. The threshold is
set at **0.50**, below published cognate performance, because this is a sanity check and
not a benchmark to be optimised against.

**Why a third readout.** `M` is built from **within-target, across-ligand** percentiles.
The self-rank is **within-ligand, across-target**. Those are rows and columns of the same
matrix and one does not imply the other, so the axis `M` rests on was never validated. The
third readout costs nothing — `07_docking/df_results_<home>.csv` and
`08_control/docking/df_results_<home>.csv` are the same receptor, box and engine, differing
only in the ligand set — and it states exactly what `M` needs: a known binder reaches the
top decile of the production library at the pocket it was crystallised in.

**Why a minimum n, and why `unknown` below it.** Without one, a single labelled ligand
ranking first on a 20-target panel reads `median 0.025` and passes. A readout over fewer
than 10 observations is reported as **`unknown`, not `yes` and not `no`**: it is not
evidence either way, and the spec already treats `unknown` as not-interpretable.

**Why a panel-coverage rule.** A ligand missing from any scored target is dropped (§ 2's
logic, applied per ligand). An incomplete **target**, however, is excluded from *every*
ligand's denominator. Simulated on a ligand genuinely ranked 30th of 482 (true 0.0613,
correctly failing): 100 targets failing **at random** gives 0.0589 and still fails; the
same 100 failing because they were **the best-scoring ones** gives **0.0013 and passes**.
The correlation is plausible in exactly that direction — the most generous pockets are the
largest, the largest have the biggest boxes, and the biggest boxes are likeliest to hit
`--unit-timeout` or OOM. At coverage 0.98 the worst case (every unscored target having
beaten the home pocket) can move the median by at most 0.0207 on a 482-target panel,
against 0.21 at the 0.79 coverage of that simulation. **The residual is stated, not
removed:** a true 0.0613 could still be read as 0.0406 under fully adversarial attrition at
0.98 coverage. The worst-case median is therefore reported alongside the estimate as
`self_rank.median_percentile_adversarial`, descriptively.

**The consequence, fixed now:**

> **If either control fails, the primary result is reported as UNINTERPRETABLE. It is not
> reinterpreted, not reported with a caveat, and not replaced by a different statistic
> chosen after seeing it.**

And:

> **An absent positive control is not a pass.** Until the redocking has been run, the
> primary result is reported as uninterpretable-unknown.

This is enforced in code, not left to prose: `polyphar stats` reads the verdict from
`inputs/positive_control.json`
(`{"median_self_rank_percentile": ..., "redock_rmsd_success_fraction": ...,
"library_top_decile_fraction": ..., "n_redocked": ..., "n_redock_rmsd": ...,
"n_library_percentile": ..., "panel_targets_scored": ..., "panel_targets_total": ...}`) and
writes an `interpretable` verdict onto **every primary and subgroup row** of `results.csv`
and into `summary.md` — `yes`, `no: <which threshold failed> — the primary result is
reported as UNINTERPRETABLE, not reinterpreted`, or `unknown: positive control not run` /
`unknown: <it does not support a verdict>`. The thresholds live in `stats.THRESHOLDS`
(`SELF_RANK_MAX`, `REDOCK_RMSD_MAX`, `REDOCK_SUCCESS_MIN`, `LIBRARY_TOP_DECILE_MIN`,
`MIN_SELF_RANK_N`, `MIN_REDOCK_N`, `MIN_LIBRARY_N`, `MIN_PANEL_SCORED_FRACTION`); changing
any of them is an amendment to this document.

**The gate validates its input.** Every readout is a fraction and must be a real number in
`[0, 1]`. A boolean, a string that parses to a number outside that range, a non-finite
value, a negative self-rank or a fraction of 99 is **not** a pass and **not** a near-miss:
the verdict is `unknown: <path> is out of range`, because an invalid positive control must
land in `unknown` and never in `yes`. The descriptive `pose_in_pocket_fraction` is
range-checked on the same terms although it gates nothing. A `source` key describing how
the numbers were produced is expected; its absence does not change the verdict but is
written to `failed.tsv`.

**A non-finite score cannot produce a passing rank.** `NaN` compares false against
everything, so a `NaN` *home* score would be counted as beaten by nothing and tied with
itself — percentile **0.0**, better than the best achievable `0.5/N` — and an all-`NaN`
score matrix would read as a perfect control. The control reads its own `df_results`
through the same reader the main analysis uses, which rejects non-finite values and names
them, and the ranking function refuses a non-finite score outright.

Where the verdict is not `yes`, every primary, subgroup and contrast p-value in the
`summary.md` arms table is printed with the suffix `(UNINTERPRETABLE)`. The numbers are
still printed — hiding a number invites recomputing it somewhere without the label — but
the label travels with them, because that table is what gets pasted into a Methods draft.

### 6b. What a pass licenses, and what it does not

The label's biases are stated wherever the control is reported, in the `caveat` key of the
result itself and in `08_control/summary.md`, because the numbers are quoted by people who
will not read this document:

1. **The home pocket is the holo conformation, induced-fit around this very ligand.**
   DoGSiteScorer runs on the pdbfixer structure, i.e. the holo backbone with the heterogens
   stripped, so the control self-docks into a receptor whose side chains are already
   arranged for the ligand it is recovering. **This is the easiest benchmark in the field
   and it is the largest bias in this section.** It says nothing about cross-docking or
   about the ~63 % of the panel that is AlphaFold.
2. The absence of a co-crystal ligand does not mean a pocket is undruggable (most such
   structures were solved apo), and the label skews toward **large** pockets and
   **well-studied** proteins.
3. **"Drug-like ligand" overstates the criterion**, which is only a 12–60 heavy-atom
   non-solvent HETATM residue in contact with the lining. Across 1,142 net2 chains the
   commonest labelled codes are **PSU 49** (pseudouridine), **7PE 22** (a PEG oligomer the
   `SKIP` list misses), **NAP 20**, **OMG 19**, **FAD 18**, **ATP 12**, **ADP 10** — about
   14 % modified RNA residues and about 10 % nucleotide cofactors, neither of which is a
   drug.
4. **The self-rank cannot detect the defect its own rationale names.** It is a relative
   within-ligand comparison and is therefore invariant to any monotone transform applied
   uniformly across targets — **including the +0.33 kcal/mol smina offset, which cancels
   exactly.**

A pass licenses exactly this: *within this panel and this engine, a molecule scores better
at the site it is known to occupy than at ~95 % of other panel sites, recovers that site's
deposited pose, and outranks the production library there — so a within-target percentile
carries site-discriminative content.* Nothing about absolute affinity, and nothing that
calibrates the `score_threshold x min_targets` surface of § 8.

It is a sanity check, not an optimisation target.

## 7. Replicate spread (no analytic standard error)

`M` is recomputed **separately per replicate** (r0, r1, r2) from each replicate's own raw
poses, and the arm statistic is reported per replicate beside the permutation p. If sign
and rough magnitude hold across all three, the result does not depend on replicate jitter.

**No analytic standard error is derived for `M`.** In particular `mean |delta| / sqrt(N)`
is not one: for a normal, `E|X - Y| = 1.128 sigma`, so that figure is biased by ~13 % and
is not the standard deviation of anything. It must not reappear.

## 8. The descriptive surface

The `score_threshold x min_targets` grid is reported as a figure/table, **descriptively,
with no p-values attached**. It is how the field states these results and it should be
shown; it is not the test, and it contributes no tests to any correction. `polyphar stats`
makes this structural, as an **allowlist**: `surface.csv` and `surface.md` may carry only
the fields of `stats.SURFACE_COLUMNS`, and **any** other field name — in the column list or
smuggled in through a row dict — is refused, whether or not it looks like a p-value. A
denylist was tried and is not tight enough: `p-value` (the commonest spelling),
`adj.P.Val` (R/limma), `Pr(>|t|)`, `pFDR`, `log10p`, `prob`, `confidence` and a Cyrillic
`р_value` all defeated it. The `PVALUE_LIKE` pattern is kept as a second line, applied to
the allowlist itself and to NFKC-normalised, homoglyph-folded names, so an expected column
can never be renamed into a p-value without the write failing first.

## 9. What is deliberately NOT pre-registered

* Which drugs are in the positive set. That is clinical-guideline reading (EASL / NCCN /
  ESMO), intersected with dockable e-Drug3D records after the stage-4 filter, with no
  approval-date cut and with antibodies and metal complexes excluded. It is supplied to
  the stage as `inputs/positive_sets.csv`, and both counts are recorded: the guideline set
  size and how many survived the dockability intersection. **Every guideline drug dropped
  is named** in `failed.tsv` and in `summary.md` — whether it was dropped for
  undockability (it did not resolve to a library record) or for coverage (it resolved but
  was not scored on enough of the panel, § 2).
* The optional larger negative drug set (approved drugs with no oncology indication at
  all, n in the hundreds). It is a strengthening to be used only if the primary result is
  marginal, and it would be a separate, labelled analysis.

## 10. Ranking on PocketMatch

Nothing in `polyphar stats` ranks on, reads or reports PocketMatch similarity. If any
ranking is ever added, it must use **`Pmax`, not `Pmin`**, or carry an explicit
minimum-size / size-ratio guard: `Pmin == 1` is near-meaningless because containment
dominates it, and about 49,600 artificial `1.0` rows would sit at the top of any
distribution built on it. A "perfect cross-protein match" must never be reported from
`Pmin` alone.

---

## Amendments

Every entry here must be dated, must say what changed and why, and must be accompanied by
a new `stats.PREREGISTRATION_SHA256`; the test suite fails until both are done together.

### 2026-11-06 — the estimator, before any result exists

`dock` has produced no live draw for any of the seven subtypes (net1's earlier 109 units
were set aside as `superseded_mispaired_20261105T162341` when the receptor/pocket pairing
defect was found), and no `08_stats` or `08_control` directory exists anywhere. No value of
the statistic existed when this was written.

* **§ 2 step 1**: ligands are ranked by **`score_min`**, not `score_mean`.

  Vina-family docking is a stochastic search over a *deterministic* scoring function. Given
  a pose the score is fixed; what varies between replicates is whether the search found the
  good pose. Every draw is therefore an upper bound on the global minimum of the scoring
  function in that box and never below it, so `E[draw] >` the quantity of interest: **the
  mean is biased by construction**, by the expected failure of the search, while the
  minimum over draws converges to it from above and is the consistent estimator.

  This is not only a precision argument. The bias is not uniform — a flexible ligand's
  search is harder, so its draws sit further above the true minimum than a rigid ligand's.
  Ranking on the mean therefore penalises exactly the molecules with the most
  conformational space. The size-matched strata (§ 2) already defend against this, since
  positives and controls are matched on heavy atoms and rotatable bonds and so on search
  difficulty too; ranking on the minimum removes the bias at its source rather than
  relying on that defence.

* **§ 2 step 1**: ligands with fewer draws than their target's best-attested count are
  dropped, counted and named. min-of-k is a function of k, so an uneven draw count is a
  bias under this estimator where it was only a nuisance under the mean.

* The replicate spread (`score_sd`, `score_se`) is **reported and never ranked on**. It is
  a sampling-adequacy diagnostic — whether the search converged — not an error bar on
  affinity. Nothing in this document treats it as one.

* Not changed: the replicate count (3), the within-target percentile, `M`, the strata, the
  permutation test, the coverage rule, the control gate.

### 2026-11-04 — revision before any result exists

Two independent reviews of the implementation (one statistical, one on the guarantees and
the wiring) were run **before `dock` had produced a single score**, so no value of the
statistic existed when the following was written. Nothing here was chosen after seeing a
number.

* **§ 2** coverage rule: 0.90 → **0.98**, and extended from proteins to **proteins and
  targets**, with the bias and `E[max of k uniforms]` arithmetic that forces both.
* **§ 3b** added: the within-stratum size limitation stated plainly and quantified, a
  per-arm/per-stratum **balance table**, and a **sensitivity analysis** on heavy-atom
  deciles — pre-registered as a sensitivity analysis, not as the primary. The test itself
  is unchanged.
* **§ 5** added: the prediction *HCC-on-HCC > HCC-on-AD* is now **tested**, as a paired
  within-stratum permutation on the difference, instead of being read off two separate
  significance verdicts (the difference-of-significance error). Also recorded: Holm makes
  the negative control harder to reject, so its **raw** p is the honest figure for that
  arm.
* **§ 6** added: the positive-control gate validates its input and an out-of-range value
  lands in `unknown`, never `yes`; the `UNINTERPRETABLE` label now travels with every
  p-value in the arms table.
* **§ 8** the no-p-values guarantee is an **allowlist**, covering `surface.md` as well as
  `surface.csv`.
* **§ 9** "every guideline drug dropped is named" now covers coverage drops, not only
  undockability.
* This document's sha256 is pinned in code (header).

### 2026-11-04 (second) — the positive control, after two reviews of the `control` stage

Written **before `dock` has run**: no v2 run folder exists on this machine, no co-crystal
ligand has been redocked, and no value of any readout in § 6 exists. Two independent
reviews of the `polyphar control` implementation found that the control **as pre-registered
could not fail**. Everything below is § 6 and § 6b only; §§ 1–5 and 7–10 are unchanged, and
the primary test is untouched.

* **§ 6 the expected n: 687 → about 24.** 687 counts labelled pockets over all 8,005
  experimental net1 pockets; the stage labels the **panel's** targets, because only a panel
  pocket can have a self-rank. Measured label rate on panel-style pockets 7/52 = 13.5 %,
  over ~178 experimental targets of 482. The reason for the change is an arithmetic
  correction, not a result.
* **§ 6 pose-in-pocket is no longer a gate.** `POSE_IN_POCKET_MIN = 0.70` sat **below its
  own null**: the docking box is the pocket's unpadded lining bounding box, the engine
  keeps the ligand inside it, and random non-clashing centroids in that box pass ~79 %
  against ~85 % for true deposited poses. It is now reported descriptively.
* **§ 6 a redocking-RMSD gate replaces it:** `>= 50 %` of the self pairs within
  **2.0 Å**, symmetry-insensitive, in place. Field-standard criterion, published base rate
  60–80 %, measured null 0 of 1,871 random non-clashing rigid placements.
* **§ 6 a third gated readout:** the co-crystal ligand's percentile among the production
  library **at its own target**, `>= 50 %` in the top decile (chance 0.10). The first two
  readouts validate the transpose of the axis `M` is built on; this one validates that
  axis.
* **§ 6 minimum n per readout (10), returning `unknown` below it.** Without one, a single
  ligand ranking first on a 20-target panel passed.
* **§ 6 a panel-coverage gate (>= 0.98)**, with `panel_targets_scored` and
  `panel_targets_total` now in the result JSON rather than only in the manifest, and the
  adversarial worst-case median reported beside the estimate. Score-correlated **target**
  attrition turned a genuine failure (true self-rank 0.0613) into a comfortable pass
  (0.0013) when the 100 best-scoring of 482 targets dropped out.
* **§ 6 a non-finite score is refused** rather than ranked. A `NaN` home score produced
  percentile 0.0 — better than best-of-N — and an all-`NaN` matrix produced the verdict
  `yes`.
* **§ 6b added**, stating what a pass licenses and the four biases that travel with the
  number, chief among them that **the home pocket is the holo conformation induced-fit
  around this very ligand**, which is the easiest benchmark in the field and was not
  previously stated anywhere.
* The thresholds now live in one `stats.THRESHOLDS` dict rather than two constants.
