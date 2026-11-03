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
M_i = mean over the 482 targets of percentile_t(i)
```

A molecule that scores moderately well on many targets has a high `M`. That is a literal
transcription of "cumulative weak binding across a set of related pockets". It is
continuous, uses every score, and **has no threshold to flip across** — so the power
question that blocked the old statistic does not arise.

Report the `threshold x min_targets` surface as a **figure, descriptively, with no
p-values attached**. It is how the field states these results and it should be shown; it
is just not the test.

### The test: permutation within size strata

- Stratify the library on heavy-atom count (quintiles) crossed with rotatable-bond count
  (tertiles), or nearest-neighbour matching on (HA, rotB) at ~5 controls per drug —
  decide once on the real library, before looking at outcomes, and record the choice.
- Statistic: mean `M` over the positive set minus mean `M` over matched controls.
- Null: permute the drug/non-drug label **within strata**, 10,000 times.
- p = proportion of permutations at least as extreme. One-sided; the claim is directional.

Exact, no distributional assumption, and the stratification is visible rather than buried
in a model's functional form.

## The negative control the old draft was missing

The audit's F2 noted there was no set of drugs expected *not* to be polypharmacological.
There is one, and it costs nothing: **the Alzheimer's drugs, tested against the HCC
panel, and the HCC drugs against the AD panel.**

Both sets are approved small molecules, both skew large relative to library median, both
get the same preparation and the same engine. If within-target ranking and size matching
work, the cross-indication arm should come out null. If the HCC drugs light up the AD
panel as brightly as their own, the statistic is measuring "approved drug" or "big
molecule", not disease-specific polypharmacology — and that is the single most
informative number the study can produce.

This should be stated as a **pre-registered prediction with a direction**, not computed
and then interpreted: HCC-on-HCC > HCC-on-AD, and AD-on-AD > AD-on-HCC.

A caveat to record rather than hide: the two panels are not independent samples of
protein space, and some targets may be shared. Report the panel overlap (accession-level)
alongside the result.

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
  corresponding guideline set for Alzheimer's;
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
| AD drugs on AD panel | primary |
| HCC drugs on AD panel | negative control |
| AD drugs on HCC panel | negative control |
| targeted / cytotoxic split | pre-planned subgroup, reported, Holm within family |

Four primary-family tests, Holm-corrected. Per-drug results (which drugs individually
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
4. the four primary tests and the Holm correction;
5. the two negative-control predictions, with direction;
6. the positive-control thresholds that would count as a pass (set them *now* — e.g.
   median self-rank in the top 5%, pose-in-pocket >= 70% — and if the controls fail,
   the primary result is reported as uninterpretable rather than quietly reinterpreted).

Point 6 is the one most likely to be skipped and the one that most determines whether
the paper is believable.

## Open choices for the user

1. **Strata or k-NN matching?** Strata are simpler to describe and to defend; k-NN
   matches more tightly but needs a distance choice. Recommendation: **quintile x tertile
   strata**, for describability.
2. **Does `M` average over all 482 targets, or over the 271 proteins (best pocket each)?**
   Averaging over targets weights proteins that contributed 2 pockets twice.
   Recommendation: **per protein, best of its pockets**, so the statistic is
   "how much of the proteome does this drug touch", which is the claim as stated.
3. **Are the AD drugs really a fair negative control**, given that some AD drugs are
   known promiscuous (donepezil, galantamine)? If so, the negative control is
   conservative, which is the safe direction — but it should be said.

## Implementation order

1. `parse_scores` fix (independent, do now).
2. Positive set rebuild — needs guideline reading, not code.
3. Pre-registration document, committed before `dock` runs.
4. `polyphar stats`: percentile matrix, `M`, matching, permutation, the four tests.
5. Positive control: redock labelled ligands, self-rank and pose-in-pocket.
