# Panel selection and GPU docking redesign

**Date:** 2026-11-03 (third rewrite; see git history for superseded versions)
**Status:** design — stages 1-3 implemented, stage 5 designed, stages 4 and 6 open
**Supersedes:** hub-and-neighbours target selection in `polyphar/stages/select.py`

This file is rewritten rather than amended whenever the numbers change. An audit found
that its worst defect had been carrying several generations of figures side by side.
**Every number here comes from `verify_panel.py` run against the shipped modules on
`disease_subtypes_v2/Liver_cancer/net1`, unless stated otherwise.**

## Intent

The claim is **drug-centric and cumulative**: an approved HCC drug may owe part of its
efficacy to weak but real binding across a set of *related* pockets, beyond its primary
target. The old pipeline selected targets as "highest-degree node of the pocket-similarity
graph, plus its neighbours", which served that claim badly — the hub is a topology
artifact and the target count tracked how many structures were processed rather than
biology.

Two parts replace it:

1. A **panel** defined by explicit pocket filters, docked in full.
2. Pocket similarity used **post hoc, on the docking results** — not to build the panel.

## The central methodological decision

**Pocket similarity never selects or groups targets.** It is used only after docking:

> Are the targets a given drug hits more similar to one another than a random set of the
> same size drawn from the same panel?

This replaced an earlier design that pre-grouped pockets into clusters, dropped (user,
2026-11-03) after it collapsed on the dockable set — see "Why the cluster layer was
dropped". The post hoc form is better on three counts: no circularity, since similarity
touches neither panel construction nor group definition; no arbitrary thresholds, since
the full `Pocketsim` distribution is a continuous measure; and a clean permutation null
over random same-size target sets.

## Panel definition

### How the filters were chosen — and chosen wrongly twice

Recorded because the failures are instructive and must not be repeated.

**First attempt: `segments <= 3`,** a count of contiguous runs in the lining-residue
numbering. It is provenance-biased — pdbfixer renumbers experimental chains 1..N with
missing loops modelled in, while AlphaFold pockets carry UniProt positions. Measured, it
passed AF pockets at 39.4% against 17.5% for experimental ones, a 21.9-point gap.

**Second attempt: `radius of gyration <= 8.5 A`,** adopted because it closed that gap to
1.5 points. It was never checked against whether it selects *druggable* pockets, which is
its job. Calibrating against the 687 of 8,005 experimental net1 pockets holding a
co-crystallised drug-like ligand:

| | ligand-occupied (n=687) | no ligand (n=7,318) |
|---|---|---|
| median Rg | **10.48** | 7.34 |

`Rg <= 8.5` kept **8%** of pockets known to bind a drug and 69% of the rest —
**anti-enriching roughly eightfold**. Real binding sites have *large* lining-residue
envelopes; a shallow surface dent has a small one. Covering 90% of ligand-bound pockets
would need `Rg <= 13.97`, which filters nothing.

**What the field actually does.** Druggability prediction is a trained-model problem, not
a geometry-threshold problem: DoGSiteScorer's `drugScore` is an SVM over size, shape and
hydrophobicity descriptors with ~88% reported accuracy, and fpocket, SiteMap, PockDrug and
DrugPred all learn a score from the same descriptor families. We had stacked four
hand-chosen geometric thresholds on top of a validated classifier.

**Third attempt, and the current one.** Filters chosen on physical grounds and stated
plainly, rather than tuned to net1 — one subtype with a known-bad history.

### The filters

| filter | what it does |
|---|---|
| `drugScore >= 0.5` | DoGSiteScorer's trained druggability SVM |
| `300 <= volume <= 1500 A^3` | the cited druggable range is ~300-1000, stricter criteria near 1000; 1500 is a deliberate allowance. **The ceiling is what rejects the oversized interdomain cavities** (~2,600 A^3 against ~270 typical); 2 of the 100 worst survive |
| subpocket parent `<= 1500 A^3` | without it, a subpocket carved from a 25,840 A^3 cleft passes the ceiling and becomes its protein's representative |
| AF pocket `pLDDT >= 70` | a sanity floor only — pLDDT does **not** detect the oversized cavities (the worst average 84.8 against 83.3 overall). It is local per-residue confidence, not domain packing; PAE would be the right metric |
| `surf/vol <= 0.45` | enclosure. Ligand-bound pockets sit at a median 0.410 against 0.510 for all pockets. It does **not** reject the oversized cavities (0.539, barely distinguishable) |
| TM policy C | drop multi-pass TM proteins; drop AF pockets overlapping a TM segment by > 10% |
| Uni-Dock grid legality | a guard, not a selector: <=128 points per axis and <=531,441 total at 0.375 A |

**Volume and `surf/vol` are both kept because they do different jobs** — volume rejects
oversized cavities and cannot judge shape; `surf/vol` judges shape and cannot separate the
oversized ones.

### Two honest caveats

- **`drugScore >= 0.5` does no work.** `DrugScore.csv`'s `Selected` column already means
  exactly that, so the step is inherited from `pockets` and constrains nothing further.
  Raise it above 0.5 for it to filter.
- **`surf/vol <= 0.45` has no literature value.** It is placed near the observed
  ligand-bound centre. Unlike Rg the direction is physically unambiguous, but it is still
  a descriptive choice and Methods must say so.

### The funnel

| step | pockets remaining |
|---|---|
| selected (drugScore >= 0.5) | 11,599 |
| volume 300-1500 | 4,705 |
| subpocket parent | 3,621 |
| AF pLDDT >= 70 | 3,458 |
| surf/vol <= 0.45 | 1,647 |
| TM policy C | 1,435 |
| Uni-Dock grid legal | 1,367 |

**Panel: 158 proteins, 1,367 eligible pockets, 32% AlphaFold** (50 AF, 108 experimental).
One row per protein, the eligible pocket with the highest drugScore, ties to the smallest
pocket name. Every rejected pocket carries a reason in `failed.tsv`, and rejections plus
eligible equals selected.

The AlphaFold share falling from 53-68% to **32%** materially reduces the interface-pocket
exposure below, which remains the panel's main uncontrolled risk.

### The docking box is unchanged

`pocket_box` (`dock.py`) takes the unpadded bounding box of the `_res_P_*.pdb` lining
residues. Deriving it from the DoGSiteScorer cavity (`_gps_`) was tested and is **worse**:
the cavity spans 18.1 A median unpadded, but Vina needs the whole ligand inside, and
cavity + 4 A already exceeds the current box. The lining-residue envelope is effectively
"cavity + ~3 A". Leave it alone.

## Ligand library

e-Drug3D, 2,083 records. Its inclusion limit is MW <= 2000, so it holds peptides far
outside Vina's domain — afamelanotide at 64 rotatable bonds, ganirelix 58, degarelix 53,
sugammadex at 1,994 Da. Including them stalled a calibration run for hours.

**Dockability filter: exclude `rotatable bonds > 15` OR `MW > 900`** — 123 of 2,083
excluded (5.9%), keeping 1,960. No HCC or AD positive-set drug is affected.

**The rotatable-bond definition is non-standard and must be pinned in Methods.**
`ligsize.py` counts single, non-ring, non-terminal heavy-atom bonds via bridge detection
and does **not** exclude amide C-N bonds; the conventional RDKit/Lipinski definition does,
and so does `obabel`. Swapping implementations silently changes which ligands are
excluded.

### Conformer preparation (stage 4, not yet implemented)

Vina-family docking samples torsions but **not ring conformations or stereochemistry**, so
the input conformer's ring pucker and stereo centres are fixed for the whole run.
e-Drug3D's coordinates are of unknown provenance and RDKit reports "tagged as 2D, but at
least one Z coordinate is not zero" for every record.

**Minimise in place with MMFF94**, preserving the library's geometry and stereochemistry.
Regenerate with ETKDG only for records that fail sanitisation, and log each one. One
conformer per ligand suffices. Preparation is deterministic, done once, and the prepared
library is hashed into the dock fingerprint.

This changes results relative to the smina baseline, so the agreement figures below
**must be re-measured after it**.

## Docking engine: Uni-Dock 1.2.0 on GPU

Validated 2026-11-03 on a Colab L4 against a finished smina run with identical receptor
PDBQTs, ligands and boxes.

- **Throughput: median 0.176 s/ligand**, flat across box volumes from 3,679 to
  20,112 A^3. CPU smina was ~200 core-seconds/ligand on 64 cores.
- **Agreement with smina: n=580, Pearson r=0.962, Spearman rho=0.951**, Uni-Dock
  +0.326 kcal/mol less negative. A shift of the score scale, not a redesign.
- **Determinism is broken** (upstream issue #133, open since 2024-05): 95-98% of ligands
  score differently across identical same-seed repeats, mean |delta| 0.07-0.13, max 1.202.

| scope | GPU-hours | units (L4 @ 1.94/hr measured) |
|---|---|---|
| 158-protein panel, 1,960 ligands, 1 replicate | 15.1 | 29 |
| all 7 subtypes, 1 replicate | 106 | 206 |
| **all 7 subtypes, 5 replicates** | 530 | **1,028** |

Against a 1,600-unit balance (~64%), units expiring after 90 days. The 7.0x subtype
multiplier is an overestimate — AD subtypes have 78/88/132/21 input accessions against
net1's 461 — so this errs safe.

**Scores are means over 5 replicates with a reported standard error** (user, 2026-11-03).

**Unresolved, and blocking the statistic.** The SE must be derived from a proper standard
deviation, not from mean |delta| / sqrt(N) (for a normal, E|X-Y| ~ 1.13 sigma). Replicate
independence is unjustified while the root cause is unsettled — `--max_gpu_memory` is an
engine flag `unidocktools` rejects, so the batch-composition hypothesis is untested. And
**the statistic is a threshold count**: `min_targets` counts targets below a cutoff, so
jitter near the cutoff flips classifications in a way averaging the score does not fix. No
power statement exists.

### Engine practical notes

- Install: `micromamba create -c conda-forge unidock` plus `pip install unidock-tools` —
  31 seconds measured, no Docker daemon needed.
- `unidock_tools` writes **one SDF per input ligand**, preserves input SD fields
  (`<ID>`, `<name>`, `<cas>`) and adds **`<docking_score>`** — not the engine's `ENERGY=`.
  Key on `<ID>`; RDKit rewrites the title.
- Use `--search_mode balance`. Do not mirror smina's exhaustiveness 12.
- Filter absurd scores (upstream issue #205 reports +300 kcal/mol on a Meeko ligand).
- **Never use Uni-Dock before 1.2.0**: earlier versions silently corrupted ~20% of results
  for boxes >= 30 A.

## Statistics (stage 6, open)

The test is panel-wide. Three things must be right.

**1. The positive set.** Currently an exact string match on one `indications` value plus
four hand-added names, one of which (ramucirumab) is an antibody that cannot be docked.
Replace with current clinical-guideline systemic therapies (EASL/NCCN/ESMO) intersected
with dockable small molecules present in e-Drug3D, reporting both counts. No approval-date
cut. Exclude antibodies and metal complexes. Report targeted agents and cytotoxics
**separately** — the set is bimodal and the halves carry opposite size biases. Alzheimer's
is in scope.

**2. Ligand size and pocket geometry.** Library median is 25 heavy atoms against 31 for
HCC drugs (apatinib 88th percentile, sorafenib 75th, fluorouracil 2nd). Vina scoring is
roughly additive over atoms, so targeted agents clear any cutoff at an elevated rate for
size reasons alone; the naive Fisher test would be significant against a random panel.
Add rotatable bonds alongside heavy-atom count.

Two routes, **not yet chosen**: residualise on a fitted
`E[score | pocket geometry, ligand size]`, or match. Residualisation currently rests on a
calibration that mismatches production on four axes — wrong engine (smina), wrong ligand
geometry (unprepared), wrong sampling design, wrong ligand range. With 48 pockets and ~6
collinear descriptors the geometry term is weakly identified regardless. **Either re-run
the calibration under production conditions, or choose matching and say so.** Note that
geometry descriptors are entangled with pocket similarity, so over-adjusting can remove
signal.

**3. Multiple testing.** The legacy `score_threshold` x `min_targets` grid over three
subtypes is thousands of tests and its best p-value is not a p-value. Either pre-register
one combination **before** looking at new results, or correct across the grid and report
the whole surface.

### The post hoc similarity test

For each drug, take its hit set and summarise pairwise pocket similarity within it from
`Pocketsim.csv`. Compare against permutation: random same-size target sets drawn from the
panel, optionally matched on pocket volume. This needs no new docking.

One inherited filter must be examined first: `parse_scores` (`similarity.py:714`) drops
pairs with `Pmax == 1.0` as a "legacy filter". Since same-accession pairs are excluded
separately, this discards *cross-protein* pairs with perfect similarity — the strongest
possible evidence for the claim. Count them among the 1,367 eligible pockets before
inheriting the rule.

### Positive control — a prerequisite

The study would otherwise report cumulative sub-threshold binding on receptors that are
32% AlphaFold models, with a new engine carrying an uncharacterised +0.33 kcal/mol offset,
and **nothing anchoring score to real binding**. The subtype gene lists do not contain the
known targets of the approved HCC drugs, so that control is unavailable.

**Redock co-crystallised ligands and check recovery.** The machinery exists: the Rg
calibration labelled 687 experimental pockets by the drug-like ligand sitting in them,
read from HETATM records in the already-downloaded extracted chains. Coordinates match the
pocket files exactly, since `pockets.py` `_author_atoms` maps pdbfixer output back with
0.001 A tolerance.

**Caveat on that label:** absence of a co-crystal ligand does not mean a pocket is not
druggable — most such structures were simply solved apo. The label is also biased toward
large pockets, because a 12-60 heavy-atom ligand needs room, and toward well-studied
targets, because crystallographers solve complexes for sites they already care about. It
is a sanity check, not an optimisation target.

## Why the cluster layer was dropped

An earlier design pre-grouped eligible pockets into Louvain communities at Pmax >= 0.7 and
tested within them. On the dockable set that layer collapsed: 73 pockets / 84 edges became
23 / 15; seven clusters became four, all at the 3-protein minimum; coverage fell from
31/209 (15%) to 10/183 (5%). The mitotic and replication groups vanished.

Two findings made it untenable even at its former size:

- **8 of 10 cluster proteins would have been docked at a pocket that is not the clustered
  pocket.** Clusters were defined over all eligible pockets while docking uses one pocket
  per protein chosen by drugScore, so the statistic would have measured binding to the
  wrong site.
- The 7 former clusters held 46 protein-slots over 31 distinct proteins (CENPI in four),
  so the same score would have been reused across tests proposed as the multiple-testing
  unit, with no correction scheme.

**The disappearance was informative.** What vanished were obligate-complex subunits
(kinetochore: NDC80, SPC24, SPC25, SKA1, CENPI, CENPK; pre-replication complex: CDC6,
CDT1, ORC1, ORC6). For a monomeric AlphaFold model of a complex subunit, the largest
compact high-drugScore groove is plausibly the **missing partner's binding surface**.
Subunits of one complex sharing exposed interface geometry is exactly what that artifact
predicts, so reading their co-clustering as validation inverted the inference.

**There is no control for interface pockets.** The AlphaFold controls target only
oversized cavities. At 32% AF this is reduced but not eliminated, and must be stated in
Methods.

Two factual errors in the earlier cluster reading, withdrawn: NDC80 has no RWD domain
(SPC24 and SPC25 do; NDC80 has a calponin-homology domain plus coiled coil), and the
homeobox proteins did not cluster separately.

## Goals, non-goals, anti-goals

**Goals.** Panel defined by explicit published filters with no reference to network
topology. Panel size a property of the filters. AlphaFold models retained with artifact
control that works. Similarity used post hoc on results. Statistics corrected for ligand
size, pocket geometry and multiple testing, with a positive control. No DoGSiteScorer
rerun.

**Non-goals.** Hub or centrality selection. Clustering of the similarity network — tested
twice, dropped; do not reintroduce. Size-conditional edge thresholds. Radius of gyration
or sequence-segment compactness — measured to be wrong. PAE-based AlphaFold filtering
(record as a future option). Changes to PocketMatch or DoGSiteScorer invocation.
Maybridge. Reconciling 2023 target counts.

**Anti-goals.** A panel whose size changes with pocket coverage. Any claim of biological
grouping arising from selection. A significant p-value explainable by ligand size or
pocket geometry. Silently dropping proteins.

## Architecture

### Stage ownership and order

`ids -> uniprot -> structures -> fix -> pockets -> annotate -> select -> similarity -> dock`

| stage | change |
|---|---|
| `uniprot` | TM segment positions into a **side file** `01_uniprot/tm.tsv` — re-running `uniprot` refreshes the PDB xref, which feeds `structures`, creating new chains and models that `pockets` would have to score. Unrecoverable inside the DS2 licence window |
| **`annotate`** (new) | `04_pockets/annotation/pocket_annotation.csv`: lining residues (range-encoded), box dimensions, parent volume, AF pocket pLDDT, plus Rg/segments/span kept as descriptive columns only. Derivable from existing outputs, so DoGSiteScorer is never re-run |
| `select` | panel filters; `final_targets.csv` (schema unchanged, so `dock` is unaffected), `eligible_pockets.csv`, `panel_funnel.tsv` |
| `similarity` | **runs after `select`**, on the panel's eligible pockets |
| `dock` | Uni-Dock, replicates, remote GPU — see `2026-11-03-dock-stage-gpu-design.md` |

**`similarity` moved after `select`** because it no longer feeds selection. On net1 that
is ~934,000 pairs instead of 67,169,845. The stage directory numbers now read
`04_pockets -> 06_targets -> 05_similarity`; they are directory names, not semantics, and
renumbering would invalidate every existing run folder.

`Is_Hub` and `Similarity` in `final_targets.csv` are meaningless now; `dock` reads
neither. They are retained empty for schema compatibility.

### Containers

Two images. `polyphar-core` (CPU, stages 00-06) needs a Debian/Ubuntu base with **i386
multiarch**, because PocketMatch uses the upstream 32-bit binaries at commit 4648bff (the
64-bit build overflows). DoGSiteScorer is never in the image — host-licensed, fed through
the `04_pockets/queue/` file queue. smina stays, reference-only.

`polyphar-dock` (CUDA 12.x **runtime** base, stage 07) carries Uni-Dock,
`unidock-tools`, RDKit for conformer prep, OpenBabel for receptor PDBQT.

**On Colab, do not use the image.** VMs are ephemeral, so Apptainer would pull and flatten
3-4 GB every session against a measured 31-second micromamba install. Three runtime paths:
Colab Pro+ (primary) via micromamba and a **pinned lockfile**; TRUBA via an Apptainer
`.sif` from the image, on a persistent filesystem; RunPod/vast.ai via the Docker image
directly. Reproducibility rests on the lockfile, which the image is also built from —
without it the Colab path drifts as conda-forge moves and replicate scores stop being
comparable across sessions.

### Remote GPU execution

POSEIDON has no GPU. Data movement is not the constraint: ~100 MB in and ~400 MB out per
subtype, against 7.9 GB `02_structures` and 9.8 GB `04_pockets`.

**The Colab operating pattern is load-bearing.** `colab exec` reliably hangs the local
client when driven from POSEIDON (three times out of three, once with
`RuntimeError: Connection was lost`) **but the kernel completes server-side and the VM
survives**. So: upload, short exec that launches with `nohup`, poll every ~2 minutes with
a tiny exec that tails the log (the poll doubles as a liveness ping), download, then
`colab stop` — a GPU session burns units while it merely exists.

The Colab CLI needs **Python >= 3.12**; `colab upload` requires both local and remote
paths; the default `--auth` is `oauth2` despite the README, and `--auth adc` works from
`gcloud auth application-default login --no-launch-browser`.

**This is designed but not sized.** Production is ~1,028 units over 530 GPU-hours on a
client that hangs, under a ~24-hour session ceiling — roughly 22 attended session-days.
TRUBA is the fallback and has never been sized: no partition, no queue wait, no validated
Apptainer build.

## Open, and blocking

1. **Residualisation or matching**, and whether the calibration is re-run under production
   conditions.
2. **The SE derivation and a power statement** for threshold flipping near the
   `min_targets` cutoff.
3. **The `Pmax == 1.0` legacy filter** — count what it discards before inheriting it.
4. **Positive-set source**, if a curated list exists.
5. **Interface-pocket exposure** at 32% AlphaFold.
6. **Thresholds validated on net1 only.** `surf/vol <= 0.45` and the volume range have not
   been checked on a second subtype, and net1 is the subtype whose 2023 gene list was
   wrong. Only net1 has v2 pocket output; the DS2 licence expires 2026-11-16.

## Implementation status

| stage | state |
|---|---|
| 1. annotate | **done** — `polyphar/stages/annotate.py`, verified on real data |
| 2. panel in select | **done** — `polyphar/panel.py` + rewritten `select.py` |
| 3. similarity after select | **done** — runs on `eligible_pockets.csv` |
| 4. ligand preparation | open, small |
| 5. dock rework | **designed** (`2026-11-03-dock-stage-gpu-design.md`), not built; the large one |
| 6. statistics | open, blocked on the decisions above |

Suite: 407 passed. `test_stage_dock.py` needs openbabel, absent on the dev laptop.
