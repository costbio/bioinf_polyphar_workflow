# Statistics: testing the cumulative weak-binding claim

**Date:** 2026-11-03
**Status:** design, not implemented
**Scope:** stage 6 of the panel redesign
(`docs/superpowers/specs/2026-11-03-panel-selection-redesign-design.md`)

## The claim, written as something that can fail

> Approved HCC drugs bind more of the HCC panel, at a useful level, than comparable
> molecules do.

"Comparable" is doing the work. Approved drugs are bigger than library median (31 vs 25
heavy atoms), and Vina-family scoring is roughly additive over atoms, so a drug beats a
threshold on more targets *for size reasons alone*. Any test that does not neutralise
that will come out positive whether or not the claim is true. The same applies to the
pockets: a large generous cavity gives every molecule a good score, so a panel that
happens to be large-pocket-heavy inflates every count.

So the design has exactly one job: build a comparison in which **ligand size and pocket
generosity cannot produce the result**, and then see whether anything is left.

## The decision that makes the rest simple: rank, don't residualise

The previous draft had two bad options — fit `E[score | geometry, size]` and residualise,
or match on size and throw data away. Residualisation needs a calibration we do not have
(the existing one mismatches production on engine, ligand geometry, sampling design and
ligand range), and with ~6 collinear geometry descriptors the geometry term is barely
identified anyway.

There is a third route that needs no model and no calibration run:

**Step 1 — rank within each target.** For every one of the 482 targets, rank all 1,960
ligands by mean score and convert to a percentile in [0, 1]. A generous pocket scores
everything well, but percentiles are uniform on every target *by construction*. Pocket
geometry is now gone — not adjusted for, **removed**. No descriptors, no fitted model, no
calibration.

**Step 2 — compare against size-matched molecules.** Within-target percentile still
carries the ligand-size effect: a big molecule ranks high on *every* target. So a drug's
summary statistic is compared not against the whole library but against library
molecules in the same size stratum.

This disposes of problem 2 in the old draft entirely, and it disposes of the "geometry
descriptors are entangled with pocket similarity, so over-adjusting removes signal"
worry — there is no adjustment to over-apply.

### The statistic: mean percentile, not a threshold count

`min_targets` — "how many targets does this drug beat cutoff X on" — has two defects. It
discards the magnitude of every score, and it is a threshold, so replicate jitter near
the cutoff flips targets in and out in a way that averaging the score does not fix. That
is why the old draft needed a power statement for threshold flipping and never had one.

Use instead, for ligand *i*:

```
M_i = mean over the 271 proteins of  max_t in protein  percentile_t(i)
```

**Per protein, best of its pockets** (decided, user 2026-11-03) — not a mean over all 482
targets, which would weight the proteins that contributed two pockets twice. The statistic
is then "how much of the proteome does this molecule touch", which is the claim as stated.

A molecule that scores moderately well on many targets has a high `M`. That is a literal
transcription of "cumulative weak binding across a set of related pockets". It is
continuous, uses every score, and **has no threshold to flip across** — so the power
question that blocked the old statistic does not arise.

Report the `threshold x min_targets` surface as a **figure, descriptively, with no
p-values attached**. It is how the field states these results and it should be shown; it
is just not the test.

### The test: permutation within size strata

- Stratify the library on heavy-atom count (quintiles) crossed with rotatable-bond count
  (tertiles) — **decided (user, 2026-11-03)**, over nearest-neighbour matching, because
  strata are easier to describe in the paper and harder to tune after the fact. Bin edges
  are computed on the library once and recorded before any outcome is looked at.
- Statistic: mean `M` over the positive set minus mean `M` over matched controls.
- Null: permute the drug/non-drug label **within strata**, 10,000 times.
- p = proportion of permutations at least as extreme. One-sided; the claim is directional.

Exact, no distributional assumption, and the stratification is visible rather than buried
in a model's functional form.

## The negative control

The audit's F2 noted there was no drug set expected *not* to be polypharmacological. Two
controls are available, and they are not equally good.

### The control that is already built in

**The library is e-Drug3D** — approved and withdrawn small-molecule drugs, not a generic
screening set. So the size-matched comparison set is *already* "other approved drugs of
the same size". Nothing extra is needed to rule out "approved drugs just score well", and
it is the contrast a reviewer asks for first. State this explicitly where the matching is
described; it is easy to miss and it is load-bearing.

### Panel specificity: HCC drugs against the AD panel

The remaining question the matched control cannot answer is whether the signal is about
*this panel* or about any panel. That needs the same drugs run against a different one:

| arm | status |
|---|---|
| HCC drugs on HCC panel | primary |
| **HCC drugs on AD panel** | **negative control — the one that matters** |
| AD drugs on AD panel | exploratory, see below |
| AD drugs on HCC panel | exploratory, see below |

Arm 2 holds the drug set fixed and changes only the panel, so nothing about the AD drugs'
own properties can contaminate it. Pre-registered prediction, with direction:
**HCC-on-HCC > HCC-on-AD.**

Report the accession-level overlap between the two panels alongside the result: the panels
are not independent samples of protein space and may share targets, which would compress
the difference.

### Why the AD-drug arms are exploratory, not primary

Three problems, in order of severity.

**1. The set is tiny.** Guideline AD therapy intersected with dockable e-Drug3D is roughly
donepezil, rivastigmine, galantamine, memantine, tacrine, possibly brexpiprazole — **n of
about 4 to 6**. The antibodies (lecanemab, donanemab) are excluded by dockability.
Benzgalantamine (2024) and the dextromethorphan-bupropion combination approved for AD
agitation on 2026-04-30 post-date e-Drug3D, though both components of the latter are in
it. A null from five molecules is weak evidence, and it is a null the design *predicts*,
which makes it close to unfalsifiable.

**2. The drugs are known multi-target, by design.** Donepezil binds sigma-1 at low
nanomolar affinity as an agonist, and its memory effect in mice is blocked by a sigma-1
antagonist, so sigma-1 is plausibly mechanism rather than an incidental off-target.
Galantamine is an allosteric nicotinic modulator in addition to AChE. The multi-target-
directed-ligand literature is anchored specifically on AD, where broad target coverage is
a deliberate design goal. These are close to the worst available choice for a
"monogamous drug" control.

**3. The size bias runs the opposite way.** Memantine is 12 heavy atoms, rivastigmine 18,
galantamine 21, donepezil 28, against a library median of 25 and HCC targeted agents at
31. The two drug sets occupy different strata, so AD-on-AD and HCC-on-HCC are not matched
on matching quality.

Problems 2 and 3 both push the AD arms toward *false positives*, so a null there would
still be informative. But with n about 5 there is not enough power to rely on it. Report
the arms with n stated and label them exploratory.

### An optional strengthening, if more power is wanted

A larger negative drug set can be built from e-Drug3D itself: approved drugs with **no
oncology indication at all** (n in the hundreds), size-matched the same way, run against
the HCC panel. No new pipeline run, no new docking beyond what the full library already
gets. Many non-oncology drugs do hit kinases, so this is conservative — which is the safe
direction. Worth doing if the primary result is marginal; not needed if it is clean.

## Positive control: does the score mean anything at all?

Unchanged from the earlier analysis, and still a prerequisite. The receptors are 63%
AlphaFold models and the engine carries an uncharacterised +0.33 kcal/mol offset from
smina. Nothing currently anchors score to real binding.

**Redock the co-crystallised ligands.** `rg_calibrate.py` already labelled 687
experimental pockets by the drug-like ligand sitting in them, read from HETATM records in
the extracted chains, with coordinates matching the pocket files (`pockets.py`
`_author_atoms` maps pdbfixer output back within 0.001 A).

Two readouts:

1. **Self-recovery rank.** Dock each labelled ligand against the whole panel. Its own
   pocket should rank far above chance. Report the median self-rank and the fraction in
   the top 1% / 5%. This is the number that says whether a percentile means anything.
2. **Pose-in-pocket.** For replicate 0's poses, what fraction of best-scoring poses have
   their centroid inside the DoGSiteScorer cavity rather than on a surface patch at the
   box edge? This is the assumption the whole pocket framing rests on and it has never
   been checked.

The label's biases must be stated where the result is: absence of a co-crystal ligand
does not mean a pocket is undruggable (most such structures were solved apo), and the
label skews toward large pockets (a 12-60 heavy-atom ligand needs room) and toward
well-studied proteins. It is a sanity check, not an optimisation target.

## The positive set

Replace the current exact string match on one `indications` value plus four hand-added
names (one of which, ramucirumab, is an antibody and cannot be docked) with:

- current clinical-guideline systemic therapies (EASL / NCCN / ESMO) for HCC, and the
  corresponding guideline set for Alzheimer's (the latter feeds the exploratory arms
  only, and its expected n of about 5 should be confirmed against e-Drug3D early, since
  it is small enough to change whether those arms are worth running);
- intersected with **dockable** small molecules present in e-Drug3D after the stage-4
  filter (rotatable bonds <= 15, MW <= 900);
- no approval-date cut;
- antibodies and metal complexes excluded;
- **targeted agents and cytotoxics reported separately.** The set is bimodal and the
  halves carry opposite size biases; pooling them hides that.

Both counts recorded — guideline set size, and how many survived the dockability
intersection. If a guideline drug is dropped for undockability, name it.

## Replicate jitter, now that the statistic is continuous

Mean percentile inherits jitter from the scores, but mildly: averaging over 482 targets
suppresses per-score noise by ~sqrt(482), and ranking is monotone so small score
perturbations usually do not move a rank at all.

Still, do not estimate it analytically. **Compute `M` separately per replicate** (r0, r1,
r2) and report the spread across replicates alongside the permutation p. If the sign and
rough magnitude hold in all three, the result does not depend on jitter. If they do not,
no SE formula would have saved it.

This also retires the old blocker: the SE derivation from mean |delta| / sqrt(N) was
wrong (for a normal, E|X-Y| ~ 1.13 sigma) and replicate independence was unjustified
while the batch-composition root cause stays unsettled. The per-replicate spread assumes
neither.

## Multiple testing

With one pre-registered statistic per (drug set, panel) pair the problem nearly
evaporates:

| test | what it is |
|---|---|
| HCC drugs on HCC panel | primary |
| HCC drugs on AD panel | primary (negative control) |
| targeted / cytotoxic split | pre-planned subgroup, Holm within family |
| AD drugs on either panel | exploratory, n about 5, reported with n |

**Two primary tests**, Holm-corrected. The AD-drug arms carry no correction because they
make no primary claim; their n is reported wherever they appear. Per-drug results (which drugs individually
show elevation) are **exploratory with Benjamini-Hochberg FDR**, labelled as such.

The `threshold x min_targets` surface carries no p-values at all, so it contributes no
tests.

## What must be fixed in `parse_scores` before any of this

`similarity.py:714` drops pairs with `Pmax == 1.0` **and** `Pmin == 1.0`. Those are
cross-protein pairs with perfect pocket similarity — the strongest possible evidence for
the claim, silently discarded. Even though the similarity test is deferred, the
descriptive output is wrong until this is fixed.

## Pre-registration

The design's value is almost entirely in being fixed before the numbers exist. Write
these down, commit them, and date the commit **before `dock` runs**:

1. the statistic (`M` = mean within-target percentile);
2. the matching scheme (strata definition or k-NN on (HA, rotB));
3. the permutation count and the one-sided direction;
4. the two primary tests and the Holm correction;
5. the negative-control prediction, with direction (HCC-on-HCC > HCC-on-AD);
6. the positive-control thresholds that would count as a pass (set them *now* — e.g.
   median self-rank in the top 5%, pose-in-pocket >= 70% — and if the controls fail,
   the primary result is reported as uninterpretable rather than quietly reinterpreted).

Point 6 is the one most likely to be skipped and the one that most determines whether
the paper is believable.

## Decisions taken (user, 2026-11-03)

1. **Matching: quintile x tertile strata**, not k-NN.
2. **`M` is per protein, best of its pockets**, not a mean over all 482 targets.
3. **The AD-drug arms are exploratory, not primary.** The panel-specificity control is
   HCC drugs against the AD panel, which holds the drugs fixed. See "The negative
   control" above for the three reasons.

## Implementation order

1. `parse_scores` fix (independent, do now).
2. Positive set rebuild — needs guideline reading, not code.
3. Pre-registration document, committed before `dock` runs.
4. `polyphar stats`: percentile matrix, `M`, matching, permutation, the four tests.
5. Positive control: redock labelled ligands, self-rank and pose-in-pocket.
