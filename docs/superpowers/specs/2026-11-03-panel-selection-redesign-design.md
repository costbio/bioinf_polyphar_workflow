# Panel selection and GPU docking redesign

**Date:** 2026-11-03 (rewritten after audit; see git history for the superseded version)
**Status:** design, approved for staged implementation
**Supersedes:** hub-and-neighbours target selection in `polyphar/stages/select.py`

This file was rewritten rather than amended. The previous version carried three
generations of numbers side by side and contradicted itself; an audit found that to be
its worst defect. **Every number here is from the same run of `reconcile.py` against
`disease_subtypes_v2/Liver_cancer/net1` unless stated otherwise.**

## Intent

The claim is **drug-centric and cumulative**: an approved HCC drug may owe part of its
efficacy to weak but real binding across a set of *related* pockets, beyond its primary
target. The old pipeline selected targets as "highest-degree node of the pocket-similarity
graph, plus its neighbours", which served that claim badly — the hub is a topology
artifact and the target count tracked pocket coverage rather than biology.

The redesign has two parts:

1. A **panel** defined by explicit pocket filters, docked in full.
2. Pocket similarity used **post hoc, on the docking results** — not to build the panel.

## The central methodological decision

**Pocket similarity is never used to select or group targets.** It is used only to ask,
after docking:

> Are the targets a given drug hits more similar to one another than a random set of the
> same size drawn from the same panel?

This replaces an earlier design that pre-grouped pockets into clusters and tested within
them. That design was dropped (user, 2026-11-03) after `reconcile.py` showed it collapsed
on the dockable set — see "Why the cluster layer was dropped". The post hoc form is better
on three counts:

- **No circularity.** Similarity touches neither panel construction nor group definition,
  so it cannot manufacture the structure it then detects.
- **No arbitrary thresholds.** No Pmax cut, no minimum group size, no community-detection
  resolution. The full `Pocketsim` distribution is used as a continuous measure.
- **A clean null.** Permutation over random same-size target sets from the panel,
  optionally matched on pocket volume. No distributional assumptions.

It also contains the interface-artifact problem (below): if a drug's hit set is enriched
for mutually similar pockets, that holds whether or not those pockets are partner-binding
surfaces. The artifact question then affects interpretation of *which* proteins appear,
not whether the effect is real.

## Panel definition

A pocket is eligible when all hold. The funnel is from `reconcile2.py` on net1:

| filter | pockets surviving |
|---|---|
| DoGSiteScorer 2 selected, drugScore >= 0.5 | 11,599 |
| volume 300-2000 A^3 | 4,992 |
| subpocket rejected if its parent exceeds 2000 A^3 | 3,780 |
| AlphaFold pockets: mean pocket pLDDT >= 70 | 3,590 |
| TM policy C: drop multi-pass (>= 2 TM segments) | 2,676 |
| TM policy C: drop pockets overlapping a TM segment > 10% | 2,672 |
| **radius of gyration <= 8.5 A** | **537** |
| Uni-Dock grid legal at 0.375 A (<=128/axis, <=531,441 total) | 537 |

**Panel: 164 proteins, 537 eligible pockets, 57% AlphaFold** (93 AF, 71 experimental).
One row per protein, the eligible pocket with the highest drugScore (ties: smallest
pocket name).

### Compactness is measured in space, not in sequence

This replaces an earlier `segments <= 3` rule (a segment break being a gap of more than 5
residues in the sorted lining-residue numbers) **and** a separate 30 A cap on the docking
box's longest side. Both are gone.

`segments` measures compactness in **sequence**, so it depends on the residue-numbering
coordinate system — and experimental pockets are numbered by pdbfixer (1..N, with missing
loops modelled in) while AlphaFold pockets carry UniProt positions. Measured on net1 that
made the filter badly provenance-dependent:

| filter, at matched stringency (645 pockets) | AF pass | experimental pass | gap |
|---|---|---|---|
| `segments <= 3` | 39.4% | 17.5% | **21.9 points** |
| **Rg <= 8.60 A** | 25.1% | 23.7% | **1.5 points** |
| box longest side <= 23.3 A | 21.8% | 25.1% | 3.3 points |
| asphericity <= 0.217 | 16.0% | 27.6% | 11.7 points |
| Rg / n_atoms^(1/3) <= 1.514 | 14.0% | 28.5% | 14.5 points |

Note the direction: the sequence rule favoured **AlphaFold** by 22 points, and the gap
survived controlling for pocket size (widest at 15-30 lining residues, where ordinary drug
pockets sit). Radius of gyration is identical across both classes by construction.

Adopting it (user, 2026-11-03) also:
- **Subsumes the box cap.** Over the eligible set the box's longest side is median 21.2 A
  and **maximum 28.8 A** — nothing exceeds 30 A, because a compact pocket bounds its own
  box. One criterion replaces two.
- **Removes the SIFTS dependency entirely** (see below).
- Drops the AlphaFold share from 68% to **57%**, easing the interface-artifact exposure.
- Leaves panel size essentially unchanged (164 against 183).

The cutoff is **8.5 A**, chosen as a round value in the physically plausible range for a
drug-binding site (Rg typically 5-8 A), not inherited from the stringency of the rule it
replaces — the earlier 8.60 A was calibrated that way and that anchor was weak. Eligible
pockets run median Rg 8.08. **Validate the cutoff on a second subtype before the
production run.**

### The docking box itself is unchanged

`pocket_box` (`dock.py`) takes the unpadded bounding box of the `_res_P_*.pdb` lining
residues. Deriving it from the DoGSiteScorer cavity (`_gps_`) was tested and is **worse**:
the cavity spans 18.1 A median unpadded, but Vina needs the whole ligand inside, and
cavity + 4 A already exceeds the current box. The lining-residue envelope is effectively
"cavity + ~3 A". Leave it alone.

### Other filter notes

- The **upper** volume bound, the Rg criterion and the subpocket-parent rule are the
  load-bearing additions. None existed before, because drugScore rewards large pockets
  and has no ceiling.
- TM policy C costs **4 pockets** over policy A. Earlier spikes used policy A; the
  discrepancy was real but harmless.
- Grid legality removes nothing once Rg is applied.
- The TM pocket-overlap check runs on AlphaFold models only, where PDB residue number
  equals UniProt position. For experimental chains only the protein-level multi-pass rule
  applies, on the assumption that a crystallised construct of a single-pass protein is the
  soluble domain. State that in Methods.

## Ligand library

Source: e-Drug3D (2,083 records). Its inclusion limit is MW <= 2000, so it contains
peptides far outside Vina's domain — afamelanotide at 64 rotatable bonds, ganirelix 58,
degarelix 53, colistimethate 47, daptomycin 39, sugammadex at 1,994 Da.

**Dockability filter: exclude `rotatable bonds > 15` OR `MW > 900`** — 123 of 2,083
excluded (5.9%), keeping 1,960. No HCC or AD positive-set drug is affected.

Two things the builder must pin:
- **The rotatable-bond definition is not standard and must be stated in Methods.**
  `ligsize.py` counts single, non-ring, non-terminal heavy-atom bonds via bridge
  detection and does **not** exclude amide C-N bonds; the conventional RDKit/Lipinski
  definition does, and `obabel` does. Swapping implementations silently changes which
  ligands are excluded. Pin one and publish the exclusion list.
- There is no class flag in the SDF (only `status`, `name`, `cas`, `ID`), so the filter
  must be computed.

### Conformer preparation

Vina-family docking samples torsions but **not ring conformations or stereochemistry**,
so the input conformer's ring pucker and stereo centres are fixed for the whole run.
e-Drug3D's coordinates are of unknown provenance, and RDKit reports "tagged as 2D, but at
least one Z coordinate is not zero" for every record.

**Minimise in place with MMFF94**, preserving the library's geometry and stereochemistry.
Regenerate with ETKDG only for records that fail sanitisation, and log each one. Fix the
2D/3D flag. One conformer per ligand suffices. Preparation is deterministic, done once,
and the prepared library is hashed into the dock fingerprint.

This changes results relative to the smina baseline, so the smina-versus-Uni-Dock
agreement figures below **must be re-measured after it**.

## Docking engine: Uni-Dock 1.2.0 on GPU

Validated 2026-11-03 on a Colab L4 against a finished smina run with identical receptor
PDBQTs, ligands and boxes.

- **Throughput: median 0.176 s/ligand**, flat across box volumes from 3,679 to
  20,112 A^3. CPU smina was ~200 core-seconds/ligand on 64 cores.
- **Agreement with smina: n=580, Pearson r=0.962, Spearman rho=0.951**, Uni-Dock
  +0.326 kcal/mol less negative. A shift of the score scale, not a redesign.
- **Determinism is broken** (upstream issue #133, open since 2024-05): 95-98% of ligands
  score differently across identical same-seed repeats, mean |delta| 0.07-0.13, max 1.202
  kcal/mol, against a ~6.4 kcal/mol score range.

| scope | GPU-hours | units (L4 @ 1.94/hr measured) |
|---|---|---|
| 164-protein panel, 1,960 ligands, 1 replicate | 15.7 | 30 |
| all 7 subtypes, 1 replicate | 110.0 | 213 |
| **all 7 subtypes, 5 replicates** | 550 | **1,067** |

Against a 1,600-unit balance (~67%), units expiring after 90 days. Caveat: the 7.0x
subtype multiplier is an overestimate — AD subtypes have 78/88/132/21 input accessions
against net1's 461 — so this errs safe.

**Scores are means over 5 replicates with a reported standard error** (user, 2026-11-03).

Unresolved, and it must be before the statistic is trusted: the SE must be derived from a
proper standard deviation, not from mean |delta| / sqrt(N) (for a normal,
E|X-Y| ~ 1.13 sigma). Replicate independence is unjustified while the root cause is
unsettled — `--max_gpu_memory` is an engine flag `unidocktools` rejects, so the
batch-composition hypothesis is untested. And **the statistic is a threshold count**:
`min_targets` counts targets below a cutoff, so jitter near the cutoff flips
classifications in a way averaging the score does not fix. The variance of the count is
driven by how many ligand-target pairs sit within ~1 sigma of the threshold. No power
statement exists.

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

## Statistics

The test is panel-wide. Three things must be right.

**1. The positive set.** Currently an exact string match on one `indications` value plus
four hand-added names, one of which (ramucirumab) is an antibody that cannot be docked.
Replace with current clinical-guideline systemic therapies (EASL/NCCN/ESMO) intersected
with dockable small molecules present in e-Drug3D, reporting both counts. No approval-date
cut. Exclude antibodies and metal complexes. Report targeted agents and cytotoxics
**separately** — the set is bimodal and the halves carry opposite size biases.
Alzheimer's is in scope and short enough to enumerate by hand.

**2. Ligand size and pocket geometry.** Library median is 25 heavy atoms against 31 for
HCC drugs (apatinib 88th percentile, sorafenib 75th, fluorouracil 2nd). SMINA/Vina scoring
is roughly additive over atoms, so the targeted agents clear any cutoff at an elevated rate
for size reasons alone; the naive Fisher test would be significant against a random panel.
Add rotatable bonds alongside heavy-atom count.

Two routes, and the choice is **not yet made**: residualise on a fitted
`E[score | pocket geometry, ligand size]`, or match. Residualisation currently rests on a
calibration that mismatches production on four axes — wrong engine (smina), wrong ligand
geometry (unprepared), wrong sampling design (the document described `build_calib.py`'s
maximin sample, which `build_calib2.py` supersedes with a box-volume-stratified one), and
wrong ligand range. With 48 pockets and ~6 collinear descriptors the geometry term is
weakly identified regardless. **Either re-run the calibration under production conditions
or choose matching and say so.** Note that geometry descriptors are entangled with pocket
similarity, so over-adjusting can remove signal.

**3. Multiple testing.** The legacy `score_threshold` x `min_targets` grid over three
subtypes is thousands of tests and its best p-value is not a p-value. Either pre-register
one combination **before** looking at new results, or correct across the grid and report
the whole surface.

### The post hoc similarity test

For each drug, take its hit set (targets beating the threshold) and compute a summary of
pairwise pocket similarity within that set from `Pocketsim.csv`. Compare against
permutation: random same-size target sets drawn from the panel, optionally matched on
pocket volume. This is the test of the study's actual claim and needs no new docking.

One inherited filter must be examined first: `parse_scores` (`similarity.py:714`) drops
pairs with `Pmax == 1.0` as a "legacy filter". Since same-accession pairs are excluded
separately, this discards *cross-protein* pairs with perfect similarity — the strongest
possible evidence for the claim. Count them among the 537 eligible pockets before
inheriting the rule.

### Positive control — a prerequisite, not optional

The study would otherwise report cumulative sub-threshold binding on receptors that are
57% AlphaFold models, with a new engine carrying an uncharacterised +0.33 kcal/mol offset,
and **nothing anchoring score to real binding**. The subtype gene lists do not contain the
known targets of the approved HCC drugs, so that control is unavailable.

**Redock co-crystallised ligands and check recovery.** `lig_name` is empty for all 11,599
pockets, so the ligands must come from HETATM records in the already-downloaded mmCIFs —
one parsing pass, no new downloads. This was previously filed as "later if wanted"; it is
promoted to a prerequisite.

## Why the cluster layer was dropped

An earlier design pre-grouped eligible pockets into Louvain communities at Pmax >= 0.7 and
tested within them. On the **dockable** set (`reconcile.py`) that layer collapses:

| | pre-cap (209-protein panel) | dockable (183) |
|---|---|---|
| cluster graph | 73 pockets, 84 edges | **23 pockets, 15 edges** |
| clusters of >= 3 proteins | 7 | **4** |
| panel coverage | 31/209 (15%) | **10/183 (5%)** |

All four survivors are exactly 3 proteins, sitting on the minimum threshold. The
mitotic/replication groups vanish entirely.

Two further findings made the layer untenable even at its former size:

- **8 of 10 cluster proteins would have been docked at a pocket that is not the clustered
  pocket.** Clusters were defined over all eligible pockets while docking uses one pocket
  per protein chosen by drugScore. The statistic would have measured binding to the wrong
  site, severing the mechanism the study rests on.
- The 7 former clusters held 46 protein-slots over 31 distinct proteins (CENPI in four),
  so the same score would have been reused across tests proposed as the multiple-testing
  unit, with no correction scheme.

**The disappearing clusters are themselves informative.** The box cap removes large
pockets, and the clusters that vanished were the obligate-complex subunits (kinetochore:
NDC80, SPC24, SPC25, SKA1, CENPI, CENPK; pre-replication complex: CDC6, CDT1, ORC1,
ORC6). For a monomeric AlphaFold model of a complex subunit, the largest compact
high-drugScore groove is plausibly the **missing partner's binding surface**. Subunits of
one complex sharing exposed interface geometry is exactly what that artifact predicts, so
the earlier reading of their co-clustering as evidence of validity inverted the inference.
The cap appears to act as an accidental interface-artifact filter.

**There is no control for interface pockets.** The AlphaFold controls target only
oversized pockets. At 57% AF this remains an open exposure for the panel, and must be
stated in Methods.

Two factual errors in the earlier cluster reading, withdrawn: NDC80 has no RWD domain
(SPC24 and SPC25 do; NDC80 has a calponin-homology domain plus coiled coil), and the
homeobox proteins did not cluster separately.

## Goals, non-goals, anti-goals

**Goals.** Panel defined by explicit published filters, with no reference to network
topology or an anchor. Panel size a property of the filters. AlphaFold models retained
with artifact control that works. Similarity used post hoc on results. Statistics
corrected for ligand size, pocket geometry and multiple testing, with a positive control.
No DoGSiteScorer rerun.

**Non-goals.** Hub or centrality selection in any form. Clustering or community detection
of the similarity network — tested twice and dropped; do not reintroduce. Size-conditional
edge thresholds. PAE-based AlphaFold filtering (record as a future option). Changes to
PocketMatch or DoGSiteScorer invocation. Maybridge. Reconciling 2023 target counts.

**Anti-goals.** A panel whose size changes with pocket coverage. Any claim of biological
grouping arising from selection. A significant p-value explainable by ligand size or
pocket geometry. Silently dropping proteins — every exclusion is recorded with its reason.

## Architecture

### Stage ownership

| stage | change | why there |
|---|---|---|
| `uniprot` | TM segment positions into a **side file** `01_uniprot/tm.tsv` | re-running `uniprot` refreshes the PDB xref, which feeds `structures`, creating new chains and models that `pockets` would have to score — unrecoverable inside the DS2 licence window |
| new annotation pass | `04_pockets/pocket_annotation.csv`: lining residues, segments, span, box dims, parent volume, AF pocket pLDDT | derivable from existing outputs, so DoGSiteScorer is never re-run; a sibling file keeps DS output pristine |
| `similarity` | apply pocket filters when building the pocket table | cleans the reported network at source |
| `select` | panel algorithm; same `final_targets.csv` schema | `dock` stays compatible |
| `dock` | Uni-Dock, replicates, remote GPU | see below |
| analysis notebook | the three statistical fixes plus the post hoc similarity test | where the inference lives |

`Is_Hub` and `Similarity` in `final_targets.csv` are meaningless now; `dock` reads
neither. Drop or blank them and document the choice.

### Containers

Two images. `polyphar-core` (CPU, stages 00-06) needs a Debian/Ubuntu base with **i386
multiarch**, because PocketMatch uses the upstream 32-bit binaries at commit 4648bff (the
64-bit build overflows). DoGSiteScorer is never in the image — host-licensed, fed through
the `04_pockets/queue/` file queue. smina stays, reference-only.

`polyphar-dock` (CUDA 12.x **runtime** base, stage 07) carries Uni-Dock, `unidock-tools`,
RDKit for conformer prep, OpenBabel for receptor PDBQT.

**On Colab, do not use the image.** VMs are ephemeral, so Apptainer would pull and flatten
3-4 GB every session against a measured 31-second micromamba install. Three runtime paths:

| host | mechanism |
|---|---|
| **Colab Pro+** (primary) | micromamba + a **pinned lockfile** |
| TRUBA | Apptainer `.sif` from the image (persistent filesystem) |
| RunPod / vast.ai | the Docker image directly |

Reproducibility rests on the **lockfile**, which the image is also built from. Without it
the Colab path drifts as conda-forge moves and replicate scores stop being comparable
across sessions. For the TRUBA path the image must be Apptainer-safe: read-only rootfs, no
fixed UID, every writable path in the bind-mounted run folder, no baked host paths,
working under `apptainer run --nv`, built from a registry (`docker://`).

### Remote GPU execution

POSEIDON has no GPU, so stage 07 runs elsewhere. Data movement is not the constraint — a
docking run needs only the selected targets' receptors and pockets, ~100 MB in and ~400 MB
out per subtype, against 7.9 GB `02_structures` and 9.8 GB `04_pockets`.

**The Colab operating pattern is load-bearing.** `colab exec` reliably hangs the local
client when driven from POSEIDON (three times out of three, once with
`RuntimeError: Connection was lost`) **but the kernel completes server-side and the VM
survives**. So: upload, short exec that launches with `nohup`, poll every ~2 minutes with
a tiny exec that tails the log (the poll doubles as a liveness ping), download, then
`colab stop` — a GPU session burns units while it merely exists.

The Colab CLI needs **Python >= 3.12**; `colab upload` requires both local and remote
paths; the default `--auth` is `oauth2` despite the README, and `--auth adc` works from
`gcloud auth application-default login --no-launch-browser`.

**This is designed but not sized.** Production is 1,191 units over 614 GPU-hours on a
client that hangs, under a ~24-hour session ceiling — roughly 26 attended session-days.
TRUBA is the fallback and has never been sized: no partition, no queue wait, and no
Apptainer build validated. Size one of them before committing.

### Dock stage rework — a stage-sized piece of work

`dock.py` is built around a deterministic CLI engine and must change substantially:

- A chunk's `.done` job key includes the smina version and assumes determinism. With a
  non-deterministic engine and N replicates it needs a replicate index, and a `.done`
  marker no longer certifies *which* scores.
- `RESULT_COLUMNS = ["Pose", "SMINA_Score", "Ligand", "ID", "Pocket"]` has no replicate or
  SE dimension.
- `tests/test_stage_dock.py::test_cpu_and_chunk_invariance` **cannot pass** with Uni-Dock
  and must be retired, along with `test_reproduces_2023_targets`.
- `run`'s preflight checks `smina --version` and will exit 2 without it.
- There is **no execution model for a remote GPU**. DoGSiteScorer was solved with an
  explicit host queue and worker; this needs an equivalent — a queue protocol,
  cross-host manifest reconciliation, and how the polling loop drives `polyphar dock`.

If a `--targets` selector is added: pruning must be **additive** (out-of-scope labels are
"not mine", not "obsolete"), because `_orphans()`/`_obsolete()` will otherwise move a
previous selection's results. Note `_stage_fingerprint` already hashes `plan.targets`, so
a selector changes the fingerprint on its own — resume will not short-circuit.

### Provenance the dock manifest must record

Image digest or lockfile hash; Uni-Dock and `unidock-tools` versions; CUDA version; **GPU
model**; `search_mode`; `--batch_size`; **replicate index and seed**; prepared-library
hash; per-score mean and SE across replicates.

### SIFTS mapping — dropped

An earlier version made this "the largest implementation risk": map experimental pocket
residues to UniProt through SIFTS, to remove a cross-class bias in `segments <= 3`.
Two measurements retired it.

**The predicted bias ran the wrong way.** The reasoning was that pdbfixer models in
missing loops, so experimental numbering is contiguous, so experimental pockets get fewer
segments and the filter favours them. Measured: experimental pockets have **more**
segments (median 5, mean 5.63) than AlphaFold ones (median 4, mean 4.33), and pass
`segments <= 3` at 17.5% against 39.4%. Restoring the real gaps by mapping to UniProt
would have pushed experimental counts higher still and **widened** the gap.

**And it was unverifiable where it mattered.** `author_segments()` falls back to a
seqres-derived range when PDB bounds are unusable, and that applies to **71% of eligible
experimental pockets and 65% of experimental panel rows** — so the pdbfixer → author →
UniProt chain breaks for roughly 213 of the 327 pockets the fix would have acted on.

Switching to Rg removes the need entirely: it is coordinate-system independent, so there
is nothing to map.

Recorded for whoever revisits this: the cached SIFTS flat file **is** present at
`02_structures/sifts/pdb_chain_uniprot.tsv.gz` and **does** carry `SP_BEG`/`SP_END`
(the offset `Protein_info.csv` discards after `structures.py:182-191`), and
`pockets.py` `_coord_key()`/`_author_atoms()` already maps pdbfixer output back to author
numbering. So if a UniProt mapping is ever needed for another purpose, it is far more
tractable than the earlier text claimed — the blocker is only the discarded `SP_BEG`.

## Open, and blocking the statistic

1. Residualisation or matching, and whether the calibration is re-run under production
   conditions.
2. The SE derivation and a power statement for threshold flipping near the `min_targets`
   cutoff.
3. The `Pmax == 1.0` legacy filter — count what it discards before inheriting it.
4. Positive-set source, if a curated list exists.
5. `MITOTIC`'s provenance — a hardcoded 23-gene set in seven spikes, found nowhere else,
   whose survival count was printed in the threshold sweeps that chose the filters. Any
   enrichment claim depending on it is not independent evidence. Document or drop.
6. Interface-pocket exposure at 68% AF.

## Staged implementation plan

Each stage ends with verification before the next begins.

1. **Annotation pass** — `pocket_annotation.csv` from existing outputs: lining residues,
   radius of gyration, box dimensions, parent-pocket volume, AF pocket pLDDT. No SIFTS.
   Gates everything.
2. **Panel in `select`** — filters, one row per protein, funnel counts, `failed.tsv`
   reasons. Verify: reproduces 164 proteins / 537 pockets on net1.
3. **Similarity-side filtering** — pocket filters before PocketMatch; re-check the ETA
   model (537 eligible pockets give ~144k pairs against 67.2M).
4. **Ligand preparation** — MMFF94 in place, deterministic, hashed.
5. **Dock stage rework** — Uni-Dock, replicates, remote GPU execution model.
6. **Statistics** — positive set, size/geometry adjustment, multiple testing, the post hoc
   similarity test, the redocking positive control.

Stage 1 is the natural starting point and depends on no open question above.
