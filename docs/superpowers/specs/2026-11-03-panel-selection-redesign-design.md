# Panel and cluster selection redesign

**Date:** 2026-11-03
**Status:** design, approved for planning
**Supersedes:** the hub-and-neighbours target selection in `polyphar/stages/select.py`

## Intent and why

The study's scientific claim is **drug-centric and cumulative**: an approved HCC drug may
owe part of its efficacy to weak but real binding across a set of *related* pockets, in
addition to its primary target — which for many drugs is not clearly identified. The
secondary claim is that the same machinery runs forward, proposing multi-target scaffolds
for discovery.

The current pipeline does not serve that claim. It is **target-centric and anchored on a
single pocket**: it takes the highest-degree node of the Pmax >= 0.7 pocket-similarity
graph as a "hub" and returns the hub plus its neighbours. The hub is a topology artifact
with no relationship to any drug's primary target, and the resulting target count is not
a stable quantity.

This redesign has three parts:

1. A **panel** defined by explicit pocket filters — the published universe of targets.
2. A **cluster** layer: communities of mutually similar pockets *within* the panel, which
   is where pocket similarity earns its place and which supplies the unit for the
   cumulative-effect statistic.
3. A **docking selector**, so which clusters actually get docked is a scheduling decision
   made at the last moment rather than baked into selection.

Note on terminology: these are **clusters**, not families. They are network communities;
homology is not established, and "family" would claim more than has been shown.

## Evidence this rests on

All figures from HCC net1, `disease_subtypes_v2/Liver_cancer/net1`, 2026-11-02/03.
Reproduction scripts were throwaway and live in `~/scratch_clust` on POSEIDON.

### The hub is not a distinguished node

- Degrees: 74, 68, 66, 63, 63, 63, 61, 61, 60, ... with 38 nodes at >= 50 (mean 8.9,
  median 4.0). The top node beats the runner-up by 9% on a smooth distribution.
- The hub is `P48448_AL3B2_HUMAN_AF-P48448-F1_P_1`: the P_1 pocket of an **AlphaFold
  model** of ALDH3B2, which has no experimental structure in the run.
- The top-*betweenness* node is a different protein entirely (5HT2C, degree 23), which is
  why the superseded `Network analysis/` notebooks disagree with the script.
- 75 targets == hub degree + 1. In 2023 the same rule gave 8-11 because the pocket stage
  had output for only ~38% of models. Higher coverage mechanically raises hub degree, so
  2023 and v2 target counts are **not comparable**.
- Among the 74 neighbours only 684 of 2,701 possible pairs are edges (25% density), and
  similarity to the hub spans 0.701-0.841, median 0.751.

### Clustering: rejected on the unfiltered network, reinstated on the filtered one

This is the most important lesson in this document, and the sequence matters.

**On the unfiltered network, clustering fails.** Louvain at resolutions 0.5-2.0 on the
Pmax >= 0.7 graph (2,594 pocket nodes, 11,585 edges, 109 components): the best community
for the mitotic/replication group is 291 pockets / 123 proteins at **5% density**, purity
**0.19**, mean internal Pmax 0.760 — identical to the network at large. Against a 6.4%
base rate that is ~3x enrichment: real but far too weak to define targets. Neither
druggability filters (internal Pmax stayed flat at 0.748-0.786 across drop-AF,
volume >= 300/600/1000 and drugScore >= 0.6/0.8) nor a size-conditional edge threshold
(purity unchanged at 0.19, fewer mitotic genes recovered) changed this.

**That network was the wrong object to cluster.** It contained only 2,594 of 11,599
pockets, biased hard to large ones (median volume 1,009 A^3 against 269.7 overall), so
every clustering test was run on a pocket set dominated by the oversized clefts later
identified as the artifacts. Filtered variants were *induced subgraphs* of that same
biased graph, not networks rebuilt from filtered pockets.

**On the filter-eligible pocket set, clustering works.** Rebuilding from the 645 eligible
compact pockets (209 proteins), pairs extracted from the existing `PocketMatch_score.txt`
so PocketMatch is not re-run: 205,127 comparable pairs, Pmax median 0.275, p99 0.577,
p99.9 0.663, max 0.851. At **Pmax >= 0.7**, 7 clusters of >= 3 proteins:

| proteins | pockets | density | internal Pmax | members |
|---|---|---|---|---|
| 6 | 9 | 0.39 | 0.751 | CDC6, CENPI, EME1, ORC1, TERT, UBE2C |
| 10 | 11 | 0.40 | 0.751 | CDT1, CENPK, DIAP3, EME1, NDC80, ORC1, PAX6, PNMA5, SKA1, USH1G |
| 9 | 9 | 0.25 | 0.748 | BUB1B, CDC6, CENPI, CENPK, EFHC2, EME1, ORC6, S100P, UBE2C |
| 12 | 12 | 0.20 | 0.716 | HXC9, MMP1, NDC80, SPC24, SPC25, PAX6, SIX1, SOLH1, PKHG4, TFDP3, TNNI3, TOP2A |
| 3 | 4 | 0.50 | 0.724 | BRDT, CENPI, S100P |
| 3 | 3 | 0.67 | 0.725 | PKH4B, SPC24, TNNC1 |
| 3 | 3 | 0.67 | 0.708 | BRDT, CENPI, EFHC2 |

**Two claims previously made here were wrong and are withdrawn** (audit, 2026-11-03):

- It said NDC80, SPC24 and SPC25 share "related RWD domains". SPC24 and SPC25 have RWD
  domains; **NDC80/HEC1 does not** (N-terminal calponin-homology domain plus coiled coil).
- It said "the homeobox proteins cluster separately, as DNA-binding grooves should". They
  do not: HXC9, PAX6 and SIX1 sit in cluster 4 **together with** MMP1, TOP2A, TNNI3,
  SPC24, SPC25, SOLH1, PKHG4 and TFDP3, and PAX6 is also in cluster 2. Cluster 4 is the
  least coherent group in the table (density 0.20, lowest internal Pmax 0.716).

What survives: CDC6 and ORC1 are both AAA+ ATPases of the pre-replication complex with
homologous nucleotide-binding folds, and CDT1 sits in the adjacent cluster.

**And the coherence that remains may be an artifact, not a validation.** The AlphaFold
controls here target only *oversized* pockets (segments, volume ceiling, subpocket rule).
There is **no control for interface pockets** — grooves that exist only because an AF
model is a monomer with its obligate partner subunit absent. The clusters are dominated by
exactly that class: NDC80, SPC24, SPC25, SKA1, CENPI, CENPK, BUB1B (kinetochore), CDC6,
CDT1, ORC1, ORC6 (pre-replication complex), mostly coiled-coil or protein-interaction
subunits with no small-molecule pharmacology. For a monomeric model of an obligate complex
subunit, the largest compact high-drugScore groove is plausibly the missing partner's
binding surface. Subunits of one complex sharing exposed interface geometry is what the
artifact predicts, so **citing their co-clustering as evidence of validity inverts the
inference.** This must be resolved before any cluster is presented as biologically
meaningful.

Coverage is the limitation: clusters of >= 3 proteins cover **31/209 panel proteins
(15%)** at Pmax >= 0.7, or 73/209 (35%) at >= 0.65. Pmax >= 0.7 is the decided threshold
(user, 2026-11-03): the tighter clusters are the defensible ones, and 0.65 reintroduces a
29-protein lump.

### pLDDT does not detect the AlphaFold artifacts

- Top-100-degree AF nodes: mean pocket pLDDT **84.8**, *above* the 83.3 all-AF average.
  Pockets >= 3000 A^3 average 80.8 against 83.7 for those < 1000 A^3.
- Mechanism: pLDDT is per-residue *local* confidence. Domains can each be predicted at
  90+ while packed wrongly relative to one another, leaving a large spurious interdomain
  void of uniformly high pLDDT. The metric that captures this is **PAE**, not pLDDT.

### The artifacts are oversized pockets

| AF pocket group | lining residues | segments | span | volume |
|---|---|---|---|---|
| all selected AF pockets | 16 | 3 | 113 | 363 A^3 |
| in the flat >= 0.7 network | 46 | 5 | 229 | 1,279 A^3 |
| top-100 by degree | **90** | **7** | 256 | **2,809 A^3** |
| volume >= 3000 A^3 | 137 | 6 | 377 | 5,054 A^3 |

A real small-molecule site is 15-25 residues in 2-4 sequence segments. The hubs are
clefts and voids. `segments <= 3` removes 83-98% of the artifact hubs; the `volume`
ceiling and the segment cap together leave 2 of 100 surviving.

### Subpockets circumvent the volume ceiling

DS2 emits `P_n` and `P_n_m`. **303 of 791 eligible pockets are subpockets, and 146 have a
parent exceeding the 2000 A^3 cap** (median parent 2,093 A^3, maximum 25,840 A^3).
Excluding subpockets whose parent exceeds the ceiling: 791 -> 645 pockets,
241 -> **209 proteins**. Without this rule the ceiling is defeated wholesale.

### Per-pocket quality filters cannot shrink a per-protein panel

drugScore 0.5 -> 0.8 discards 77% of pockets but only 12% of the panel (424 -> 373),
because a protein is kept if *any* of its pockets passes. The DoGSiteScorer drugScore
cutoff the 2023 method rests on therefore selects pockets, not targets.

### Transmembrane proteins: the artifact is multi-pass, not single-pass

- Of 424 proteins with pockets, 87 have >= 1 TM segment and 48 are multi-pass.
- Dropping *pockets* that overlap a TM segment by > 10% removes only **2 proteins**: the
  single-pass proteins' selected pockets are in soluble domains.
- Dropping every TM protein costs 21 proteins for no measured benefit, including CAH9
  (carbonic anhydrase IX, a validated cancer target whose druggable site is the
  extracellular catalytic domain), INSRR, EPHA5 and PTK7.

### Drug-side statistics are confounded by ligand size

- eDrug3D heavy-atom counts: library median 25, HCC-drug median 31. Percentiles: apatinib
  88%, doxorubicin 87%, cabozantinib 85%, regorafenib 78%, sorafenib 75%, lenvatinib 68%,
  sunitinib/erlotinib 64%.
- SMINA scoring is approximately additive over atoms, so the targeted agents clear any
  score cutoff at an elevated rate **for size reasons alone**. As written, the Fisher test
  would be significant against a randomly chosen panel.
- The positive set is also **bimodal** (fluorouracil 2nd percentile, gemcitabine 22nd), so
  cytotoxics are systematically missed while kinase inhibitors are inflated.

## Goals

1. Define the panel by **explicit, published pocket filters**, with no reference to
   network topology, centrality or an anchor.
2. Make panel size a property of the filters, not of how many structures were processed.
3. Keep AlphaFold models, with artifact control that actually works.
4. Give pocket similarity a real role: **clusters within the panel**, as the unit of the
   cumulative-effect claim.
5. Make the docking scope a run-time selector, so cost is incurred incrementally.
6. Make the TM policy configurable with a defensible default.
7. Clean the reported similarity network at source, by filtering pockets before
   `similarity` rather than post-processing edges.
8. Fix the drug-side statistics: positive set, ligand-size adjustment, pocket-geometry
   adjustment, multiple testing — with a **built-in negative control**.
9. Require **no DoGSiteScorer rerun** (the DS2 evaluation license expires 2026-11-16).

## Non-goals

- Hub or centrality-based selection in any form.
- Clustering the **unfiltered** network. It was tested at four resolutions and does not
  work; do not reintroduce it.
- A size-conditional edge threshold. Superseded by pocket-level size filters, which are
  simpler and correct the same confound at source.
- PAE-based AlphaFold filtering. Theoretically the right metric for interdomain voids, but
  the segment and volume caps already remove the artifacts. Record as a future option.
- Changes to PocketMatch, the `similarity` computation itself, or `pockets` /
  DoGSiteScorer invocation.
- The Maybridge library.
- Reconciling 2023 target counts. They are not comparable; report that, do not reconcile.

## Anti-goals

- A panel whose size changes when pocket coverage changes.
- Any selection rule requiring a defence of an arbitrary anchor protein.
- Calling clusters "families", or claiming homology that has not been shown.
- A significant p-value a reviewer can explain by ligand size or pocket geometry.
- Silently dropping proteins. Every exclusion must be recorded with its reason.

## Decisions already made

By the user, 2026-11-03:

- **Panel by filters, not by hub.** Hub selection is removed, not parameterised.
- **Keep AlphaFold models** in the panel.
- **TM policy C is the default**: drop multi-pass TM proteins, and drop pockets that
  overlap a TM segment. **Policy B** (drop every protein with >= 1 TM segment) must be
  available as an option and reported as a sensitivity analysis.
- **Clusters, not families**, in naming and in prose.
- **Cluster threshold Pmax >= 0.7.**
- **Docking scope is a `dock` argument**, decided per run, not fixed at selection.
- **A negative control arm is part of the experiment**, with controls matched on pocket
  volume and residue count at minimum.
- **Do not reuse 2023 docking scores** for calibration. Measure on current data.
- **`pockets` runs first for all remaining subtypes**, before any of this is implemented.
  The pocket data is the only deadline-bound, irreplaceable part.
- **Sensitivity panels are not docked separately.** Both proposed ones are subsets of the
  primary panel, so they are obtained by filtering the primary docking results.
- **Experimental pocket residues are mapped through SIFTS** before segment/span are
  measured, fixing the coordinate-system bias below.
- **The primary-target positive control is dropped.** The subtype gene lists do not
  contain the known targets of the approved HCC drugs, and that is accepted rather than
  worked around.
- **Alzheimer's is in scope** for the positive-set rebuild.
- **Positive set** = current clinical-guideline systemic therapies, intersected with
  dockable small molecules actually present in eDrug3D, with both counts reported. No
  approval-date cut.
- **Library-level dockability filter: exclude ligands with `rotatable bonds > 15` OR
  `MW > 900`.** Applies to the real analysis, not only the calibration.
- **The docking engine changes from smina to Uni-Dock 1.2.0 on GPU** (validated
  2026-11-03; see "Docking engine" below). The 2023-continuity requirement is dropped and
  a change of ligand protonation is accepted, which is what made this possible.
- **Scores are means over 5 replicates with a reported standard error** (user, 2026-11-03;
  raised from 3), replacing the old "`--nprocs` does not change results" invariant.
- **Box-size constraints stay in place** (user, 2026-11-03): a pocket must have its
  longest box side <= 30 A *and* be grid-legal for Uni-Dock at 0.375 A spacing. Losing
  proteins to this is accepted. **The panel is therefore 183 proteins, not 209.**
- **Ligand conformers are prepared before docking** (user, 2026-11-03). See "Ligand
  conformer preparation" below.
- **The pose-in-pocket validation is run before the production screen.** It is the one
  test that can still overturn a design decision.
- **Colab Pro+ is the GPU host** for now, driven from POSEIDON. TRUBA (Apptainer +
  SLURM) is the fallback if attended sessions become impractical.

Earlier, still binding (see `CLAUDE.md` and the modernization plan): DS2 is the primary
pocket tool; DoGSiteScorer runs on the host only; stages run manually per subtype.

**My calls, not the user's** — a sensible builder could change these:

- Filter thresholds: `volume 300-2000 A^3`, `segments <= 3` (break at
  `next - prev > 5`), `drugScore >= 0.5`, AF pockets `mean pocket pLDDT >= 70`, TM pocket
  overlap `> 10%` rejects, and subpockets rejected when the parent exceeds the ceiling.
  On net1 this gives **209 proteins / 645 eligible pockets**.
- Keeping `pLDDT >= 70` at all, given it does not detect artifacts: it is a cheap sanity
  floor (86% of pockets pass) and worth a line in Methods.
- TM segments defined by UniProt `ft_transmem` only, parsing `TRANSMEM a..b` — not
  intramembrane or beta-barrel features. This produced the 87/48 split.
- Lining residues = distinct residue numbers in the DS2 `_res_P_*.pdb` file, not the
  descriptor `sumAA` column.
- Pocket pLDDT = mean over lining **residues**, not atoms (atom weighting favours large
  residues), read back from the `02_structures/` AF model because DS2 zeroes the
  B-factor column in its pocket files.
- **TM policy is not part of the similarity-side filter**, so switching to policy B never
  forces a PocketMatch re-run.
- **Panel membership does not require network membership.** Filters define the panel; a
  pocket that `similarity` dropped for its own reasons (cabbage-core rejection, the
  32-character name guard) still belongs to the panel, and `select` may therefore emit a
  target absent from the network it also writes.
- **Exclusion reporting stays per-stage.** An end-to-end funnel from the input gene list
  is a separate reporting script, out of scope.
- 209 is a documented expected value, **not a hard acceptance gate**. Recalibrate on the
  first real run and say what deviation is tolerable before stopping.

## Constraints

- **The DS2 evaluation license expires 2026-11-16.** Nothing here may require re-running
  `pockets`. All new per-pocket annotation must be derivable from existing `04_pockets/`
  output, `02_structures/` models and UniProt.
- **TM data must go in a side file**, `01_uniprot/tm.tsv` from a TM-only fetch keyed by
  accession — *not* new `uniprot.csv` columns. Re-running `uniprot` refreshes the `PDB`
  cross-reference column from a current release; `structures` takes its whole work list
  from that column and rebuilds `Protein_info.csv` every run, so new cross-references
  produce new chains and models that `pockets` would have to score, and the new
  `Protein_info.csv` hash invalidates the `similarity` fingerprint. Inside the license
  window that is unrecoverable.
- Stage fingerprints do not include the code version, and resume treats existing outputs
  as done. Use a fresh run folder or delete the affected stage directories.
- `dock` consumes only `Protein`, `Pocket`, `Receptor_Path`, `Pocket_Path` from
  `final_targets.csv` (`dock.py:388`), and has no prerequisite-manifest guard.
- Target machine: 64 cores; dock calibrated at 57 s/ligand. Run folders are on an
  external disk, so per-file I/O is slow (scanning 11,599 pocket files takes ~4 minutes).
- Seven subtypes: HCC net1/net2/net3, AD Cluster1-4. `Alzhemier` stays misspelled.

## Design

### Where each piece lives

| Stage | Change | Why there |
|---|---|---|
| `uniprot` | TM segment **positions** and subcellular location into `01_uniprot/tm.tsv` via a TM-only fetch | side file avoids the PDB-refresh cascade above |
| new annotation pass | `04_pockets/pocket_annotation.csv` keyed by `Pocket Filename`: lining-residue count, segment count, sequence span, parent-pocket volume, and (AF only) mean/min pocket pLDDT | derivable from existing outputs, so DoGSiteScorer is never re-run; a sibling file keeps DS output pristine |
| `similarity` | apply the pocket filters when building the pocket table | cleans the reported network at source |
| `select` | panel algorithm + cluster detection; same `final_targets.csv` schema plus a `Cluster` column; new `clusters.csv` | `dock` stays compatible |
| `dock` | `--targets` selector | docking scope becomes a run-time decision |
| `docking_screen_analysis.ipynb` | statistical fixes | where the inference lives |

Annotation is a separate idempotent pass rather than part of `pockets`, because `pockets`
cannot be re-run inside the license window. Whether it is a `polyphar annotate`
subcommand or a cached lazy computation is the builder's call; it must be cached.

### Pocket eligibility

A pocket is eligible when all hold:

1. `drugScore >= 0.5` (DS2; existing `Selected` column)
2. `300 <= volume <= 2000` A^3
3. `segments <= 3`
4. if a subpocket (`P_n_m`): its parent `P_n` volume must also be `<= 2000` A^3
5. if an AlphaFold model: `mean pocket pLDDT >= 70`
6. TM policy does not reject it

The **upper** volume bound, the segment cap and the subpocket rule are the load-bearing
additions. None exists today, because drugScore rewards large pockets and has no ceiling.

### TM policy

- **C (default).** Reject every pocket of a multi-pass protein (>= 2 TM segments). For
  remaining proteins with >= 1 TM segment, reject a pocket whose lining residues overlap
  a TM segment by more than 10%.
- **B (option).** Reject every pocket of any protein with >= 1 TM segment.
- **A (option, sensitivity only).** Reject multi-pass proteins only.

Overlap is in UniProt coordinates. For AlphaFold models the PDB residue number equals the
UniProt position, so this is direct. For experimental chains the design **assumes** the
crystallised construct of a single-pass protein is the soluble domain with the TM helix
removed, and applies only the protein-level multi-pass rule. State that in Methods.

### Coordinate systems: experimental pockets must be mapped through SIFTS

`pockets.py` numbers experimental pocket residues as in the pdbfixer output (renumbered
1..N), and `fix.py` runs `findMissingResidues` + `addMissingAtoms`, so missing loops are
modelled in and numbering is contiguous. AlphaFold pockets skip pdbfixer and are numbered
in UniProt positions.

Gap-filling therefore **systematically lowers segment counts for experimental chains**
relative to AF ones — a cross-class bias on `segments <= 3`, which gates both panel and
cluster membership, in a panel that is ~75% AF.

**Decided (user, 2026-11-03): map experimental pocket residues to UniProt positions
through SIFTS before measuring segments and span.** This also makes the TM pocket-overlap
check available for experimental chains, so the soluble-construct assumption can be
dropped.

**This section was mis-scoped and is rewritten** (audit, 2026-11-03). The task is more
tractable than claimed, and the real blocker is different from the one stated.

Already done, and not to be rebuilt:
- **pdbfixer → author numbering exists.** `pockets.py` has `_coord_key()` and
  `_author_atoms()` (L550-592): a heavy-atom coordinate-key map from the pdbfixer output
  back to author residue numbers in the `02_structures` extracted chain, with 0.001 A
  tolerance. `filter_segments` already uses it for the fusion-chain rule.
- **SIFTS does not need re-fetching.** The flat file is cached in the run folder at
  `02_structures/sifts/pdb_chain_uniprot.tsv.gz` (`structures.py:47-52`). Reading it in
  the annotation pass touches neither `uniprot.csv` nor `Protein_info.csv`, so the
  re-fetch hazard described earlier is moot.

**The actual blocker:** `Protein_info.csv`'s `SIFTS segments` stores **author bounds
only** (`PDB_BEG:PDB_END`). `SP_BEG`/`SP_END` are read in `structures.py:182-191` to
compute coverage and then discarded, so the author→UniProt offset is not recoverable from
`Protein_info.csv`. Either read it from the cached flat file, or add `SP_BEG`/`SP_END` to
`INFO_COLUMNS`.

Two correctness caveats:
- The flat file is **segment-range-level**. Within a segment the mapping is affine only
  with no insertion codes, no microheterogeneity and monotonic author numbering.
  Residue-level fidelity needs the per-entry SIFTS XML, a different download. The builder
  must state which is used.
- `author_segments()` falls back to a **seqres**-derived range when PDB bounds are
  unusable (`SIFTS segments source = seqres`), and `filter_segments` already counts those
  chains as `segment_unverified`. For them the UniProt offset is unrecoverable by any
  route. **Measure that fraction for net1 first** — `segments <= 3` gates both panel and
  cluster membership, so an unmappable fraction is a silent hole in the cross-class bias
  fix this section exists to deliver.

### Panel algorithm

1. Read the DrugScore table, the annotation sibling file and the TM side file.
2. Keep eligible pockets.
3. Group by UniProt accession; **one row per protein**, choosing the eligible pocket with
   the highest drugScore (ties: smallest pocket name, for determinism).
4. Apply `--exclude` as today (default O15350/p73).
5. Write `final_targets.csv` with the existing columns plus `Cluster`.

`Is_Hub` and `Similarity` are now meaningless. `dock` reads neither; the builder should
drop or blank them and document the choice.

### Cluster detection

1. Build the graph over **eligible pockets** (not one-per-protein), edges where
   `Pmax >= 0.7`, excluding same-accession pairs. Pairs come from the existing
   `PocketMatch_score.txt`, so PocketMatch is not re-run.
2. Louvain with a fixed seed, resolution 1.0 (`networkx.algorithms.community`,
   networkx >= 3.0 is already available).
3. Keep clusters spanning **>= 3 distinct proteins**.
4. Write `clusters.csv`: cluster id, protein count, pocket count, density, mean internal
   Pmax, members. Tag each panel row with its `Cluster` (blank when not in one).

A protein may appear in several clusters via different pockets — legitimate, but the
statistic must **count proteins, not pockets**, one pocket per protein per cluster.

`centrality.csv` and `unique_pairs.csv` remain descriptive outputs.

### `dock --targets` selector

`--targets all | clusters | cluster:<id> | file:<path>`, default `all`.

Two hazards the builder must handle:

- **Pruning must be additive.** `dock` currently moves outputs of labels absent from
  `final_targets.csv` into `07_docking/obsolete/<timestamp>/` (`dock.py:44`) — moved, not
  deleted, so nothing is destroyed, but docking `cluster:1` then `cluster:2` would shunt
  cluster 1's results out of the active set and the manifest. Out-of-scope labels must be
  treated as "not mine", not "obsolete", and pruning must consider only the selected
  scope.
- ~~The selector must enter `_stage_fingerprint`'s `opts`.~~ **Withdrawn — this was a
  misdiagnosis** (audit, 2026-11-03). `_stage_fingerprint` (`dock.py:625-643`) hashes
  `{"final_targets_sha256", "targets": rows, "library_md5", **opts}`, where `rows`
  iterates `plan.targets`. A selector that reduces `plan.targets` **already changes the
  fingerprint**, so the stated failure mode does not exist.

  **The real hazard is the inverse and is worse:** because the fingerprint *does* change,
  resume will not short-circuit, and `_orphans()` (`dock.py:657-673`) will classify the
  previous selection's `df_results_<label>_fixed.csv` and `<label>/` directories as
  orphans for `_obsolete()` to move. That is the same additive-pruning problem below,
  reached from the other direction.

### The docking box: leave `pocket_box` alone, clamp the tail

Investigated 2026-11-03 after the first calibration run ran ~3.7x slower than the
benchmark predicted. Conclusions, in order of importance:

**The current box derivation is sound and should not change.** `pocket_box`
(`dock.py`) takes the **unpadded** bounding box of all atoms in the `_res_P_*.pdb` file,
which holds the complete lining residues. Over the 631 eligible net1 pockets that gives a
median box of **9,550 A^3** with a median longest side of **24.6 A** — squarely inside
normal Vina practice.

**Deriving the box from the cavity would be worse.** DoGSiteScorer writes three files per
pocket: `_res_` (complete lining residues), `_atms_` (lining atoms) and `_gps_` (cavity
grid points). The cavity spans 18.1 A (median, unpadded), but Vina needs the whole ligand
inside the box, so a cavity box needs 4-8 A of padding per side:

| box source | median | p90 | max | median longest side |
|---|---|---|---|---|
| current `_res_` envelope, unpadded | **9,550 A^3** | 24,691 | 57,067 | 24.6 A |
| cavity `_gps_`, unpadded | 2,970 A^3 | 13,166 | 38,970 | 18.1 A |
| cavity + 4 A per side | 11,391 A^3 | 31,777 | 74,168 | 26.1 A |
| cavity + 6 A per side | 18,653 A^3 | 45,814 | 97,700 | 30.1 A |

Cavity + 4 A is already *larger* than the current box, and cavity + 6 A is double it. The
lining-residue envelope is effectively "cavity + ~3 A": side chains extend past the cavity
by about the margin a ligand needs. Switching to `_gps_` costs 0.8x-0.5x throughput for no
gain in validity.

**Worth doing: clamp the tail.** 18% of cavities exceed 25 A in their longest dimension
and 9% exceed 30 A; those drive the p90 of 24,691 A^3 and carry a real validity risk — a
box with a 45 A side can swallow a neighbouring groove, so a good score need not involve
the selected pocket at all, which would sever the link between pocket similarity and
binding. Clamp any box dimension at ~30 A about the pocket centroid and **flag** pockets
whose true extent exceeds the clamp rather than silently truncating them.

**Worth doing: give the ETA model a box-volume term.** `benchmark` produces a single
seconds-per-ligand figure and cannot predict a 6x range in box volume, so every `dock`
ETA is wrong. This matters before committing server time to seven subtypes.

**Worth doing: carry box volume as a covariate in the score model.** More search space
means more chance of a spuriously good pose, and box volume could correlate with cluster
membership, so it is a confounder for the cluster-versus-control test.

### Docking engine: Uni-Dock 1.2.0 on GPU

Validated 2026-11-03 on a Colab L4 against a finished smina run, using identical
receptor PDBQTs, identical ligands and identical boxes. Logs: `unidock_tests.log` and
`unidock_results.json` on POSEIDON; bundle `unidock_test.tgz`.

**Why change engines at all.** smina 2020.12.10 is based on AutoDock Vina 1.1.2 and is
the end of its line (last SourceForge upload 2021-12-01). On 64 CPU cores it measured
~200 core-seconds per ligand, which puts the full 209-protein panel at ~15 days per
subtype and ~61 days for all seven — infeasible, and it was forcing the design down to
clusters-plus-controls only. The user dropped the 2023-continuity requirement and
accepted a protonation change, removing both objections.

**Measured throughput: median 0.176 s/ligand on one L4.**

| target | box | s/ligand |
|---|---|---|
| 7UNG_W_P_11 | 3,679 A^3 | 0.180 |
| 6LK0_A_P_2_1 | 7,897 A^3 | 0.174 |
| 2JT0_A_P_1 | 13,031 A^3 | 0.177 |
| 2KMU_A_P_1 | 20,112 A^3 | 0.157 |

| scope | GPU-hours | compute units (L4 @ 1.94/hr, measured) |
|---|---|---|
| 183-protein panel, 1,960 ligands, 1 replicate | 17.5 | 34 |
| all 7 subtypes, 1 replicate | 122.8 | 238 |
| **all 7 subtypes, 5 replicates** | 614 | **1,191** |

Against a 1,600-unit balance, so the production screen consumes about **74%** of it.
That is feasible but leaves little margin for reruns, and units expire after 90 days.
If the margin becomes a problem, the fallbacks in order are: drop to 3 replicates
(714 units), run 5 replicates only for ligand-target pairs within ~0.3 kcal/mol of a
decision threshold, or move production to TRUBA. The validation run itself cost 0.91
units, so exploratory work is effectively free.

Panel size note: 183 is the count after **both** the 30 A longest-side cap and
Uni-Dock grid-legality at 0.375 A spacing (492 eligible pockets). Grid-legality alone
would give 203 and the cap alone 183; the user accepted losing the proteins.

### Ligand conformer preparation

Vina-family docking samples torsions but **not ring conformations or stereochemistry**,
so the input conformer's ring pucker and stereo centres are fixed for the whole run. A
strained or wrong input geometry therefore cannot be recovered by any amount of search,
and e-Drug3D's own 3D coordinates are of unknown provenance and quality. RDKit also
reports "molecule is tagged as 2D, but at least one Z coordinate is not zero" for every
record in the library, i.e. the dimensionality flag disagrees with the coordinates.

**Decision (user, 2026-11-03): prepare ligand conformers before docking.**

Recommended approach, to be confirmed by the builder:
- **Minimise in place with MMFF94** rather than regenerating from SMILES, so the
  library's curated geometry and stereochemistry are preserved. Regenerate with ETKDG
  only for records that fail sanitisation or minimisation, and log every such record.
- Fix the 2D/3D flag so downstream tools do not mis-handle the coordinates.
- One conformer per ligand is sufficient for Vina-type search, which explores torsions
  itself; multiple ring conformers would only matter for macrocycles, which the
  `rot > 15 or MW > 900` filter already removes.
- `unidocktools ligandprep` exists and may do this; compare it against an explicit RDKit
  step before choosing, and record which was used.
- Preparation is **deterministic and done once**, with the prepared library hashed into
  the dock stage fingerprint, so it never becomes a source of run-to-run variation.

This is a methodological change from the 2023 pipeline and from the smina baseline, both
of which docked e-Drug3D's shipped coordinates. The smina-versus-Uni-Dock agreement
figures above were measured on *unprepared* ligands and will need re-measuring after
this change.

**Throughput is flat across box volume** (0.180 s/ligand at 3,679 A^3 against 0.157 at
20,112). On CPU smina box volume was the dominant cost driver; on GPU it is noise. So the
*performance* argument for capping box size is gone and only the pose-validity argument
remains — see the box section above, and the open test below.

**Agreement with smina is strong**, better than expected: n=580, Pearson r=0.962,
**Spearman rho=0.951**, with Uni-Dock **+0.326 kcal/mol less negative** (medians -5.73
against -6.06). Uni-Dock uses the published Vina weights but a reimplemented CUDA
evaluator and a different sampling regime (exhaustiveness 384 at `balance`, not 12), so
this is a shift of the score scale, not a redesign of it. Any threshold carried over from
smina needs re-centring by about +0.33 kcal/mol. Caveat: n=580 is enough to establish the
relationship, not to calibrate the offset precisely or to show it is constant across the
score range.

**Determinism: issue #133 is confirmed and quantified.**

| target | ligands differing across 3 same-seed repeats | mean abs delta | max abs delta |
|---|---|---|---|
| 7UNG_W_P_11 | 138/145 (95%) | 0.069 | 0.626 |
| 6LK0_A_P_2_1 | 142/145 (98%) | 0.134 | 1.202 |
| 2JT0_A_P_1 | 142/145 (98%) | 0.094 | 0.819 |

Essentially every ligand's score moves between identical runs, by ~0.1 kcal/mol on
average against a ~6.4 kcal/mol score range, with excursions to 1.2. This matters more
here than for most users, because the study's claim concerns *weak* cumulative binding and
jitter of that size is comparable to the signal.

**Decision: run 5 replicates and average** (user, 2026-11-03), reporting a per-score
standard error and propagating it into the cluster-versus-control test.

**The SE reasoning used here was wrong in kind** (audit, 2026-11-03) and must be redone.
It divided a *mean absolute delta across repeats* by sqrt(N); that statistic is not a
standard deviation (for a normal, E|X-Y| ~ 1.13 sigma), and the quantity wanted is
sigma/sqrt(N). It also assumes replicate independence, which the document itself cannot
justify while the root cause is unsettled — if the jitter is batch-composition dependent
rather than floating-point, replicates in the same batch configuration are correlated and
sqrt(N) does not apply. And averaging does not control the tail: max |delta| was
1.202 kcal/mol.

**A deeper problem the averaging does not touch:** the statistic is a *threshold count*
(`min_targets` counts targets below a score cutoff), not a score. Jitter near a hard
cutoff flips classifications, and averaging the score does not fix that in any simple way
— the variance of the count is driven by how many ligand-target pairs sit within ~1 sigma
of the threshold. No power statement exists, so "0.06 is enough" is unsupported, and the
document elsewhere concedes jitter of this size "is comparable to the signal".
The pipeline's reproducibility guarantee is restated: not "results do not depend on
`--nprocs`" but **"scores are means over N replicates with a reported SE"**. This is
better science than a single deterministic-looking number, and 3 replicates for all seven
subtypes is affordable (814 of 1,600 units).

Root cause is **not** settled. `--max_gpu_memory` is an engine flag that `unidocktools`
rejects as unrecognised, so the VRAM-driven batch-composition hypothesis could not be
tested through the wrapper. With only 150 ligands everything fits one batch (default
`--batch_size` 18000) and the jitter persisted anyway, which points at floating-point
non-associativity rather than batch composition — but settling it needs the `unidock`
engine binary called directly.

**Practical notes for the builder.**
- Install: `micromamba create -p <prefix> -c conda-forge unidock` plus
  `pip install unidock-tools` — 31 seconds, no Docker daemon required.
- `unidock_tools` writes **one SDF per input ligand** (`<stem>_<i>.sdf`) in the save dir
  plus a `results.csv`. It preserves the input SD fields (`<ID>`, `<name>`, `<cas>`,
  `<status>`) and adds **`<docking_score>`** — not the engine's `ENERGY=` tag. Key on the
  `<ID>` SD field; RDKit rewrites the molecule title.
- Use `--search_mode balance`. Do not try to mirror smina's exhaustiveness 12; the engine
  assumes high exhaustiveness with capped steps.
- Add a sanity filter for absurd scores (open issue #205 reports +300 kcal/mol on a
  Meeko-prepared ligand). None appeared in 580 pairs here.
- RDKit emits "molecule is tagged as 2D, but at least one Z coordinate is not zero" on
  these ligands. Appears benign; confirm before production.

### GPU host and the Colab operating pattern

Colab Pro+ is the host, driven from POSEIDON (which has no GPU). The operating pattern
below is **load-bearing** and should be reused for the production run.

`colab exec` **reliably hangs the local client** when driven from POSEIDON — observed
three times out of three, once with an explicit `RuntimeError: Connection was lost` —
**but the kernel completes the work server-side and the VM survives** (session stays
allocated, status returns to IDLE). So never hold a long `colab exec`. Instead:

1. `colab upload <local> <remote>` (both paths are required)
2. a short exec that unpacks and launches the job with `nohup ... &`
3. poll every ~2 minutes with a tiny exec that tails the log; the poll doubles as a
   liveness ping so Colab does not reclaim the VM as idle
4. `colab download` the results, then **`colab stop -s <name>`** — a GPU session burns
   compute units while it merely exists

Scripts on POSEIDON: `launch.py`, `poll.py`, `driver.py`. The Colab CLI needs
**Python >= 3.12** (conda env `colabcli`, symlinked into `~/.local/bin`, which is on the
login-shell PATH). Despite its README the default `--auth` is `oauth2`; `--auth adc` works
from `gcloud auth application-default login --no-launch-browser`, which is the right
choice for a headless host. GPU allocation succeeded for L4 on request.

Data movement is not a constraint: a docking run needs only the selected targets'
receptor and pocket files — about 100 MB in and 400 MB out per subtype, against the
7.9 GB `02_structures` and 9.8 GB `04_pockets` of the full stage outputs.

### The experiment: clusters versus matched controls

Dock each cluster plus a **matched set of non-cluster panel proteins of the same size**,
and test:

> Does an approved HCC drug accumulate more sub-threshold binding within a cluster than
> across an equal number of unrelated druggable pockets?

This has a built-in negative control, tests the cumulative-effect hypothesis directly,
and holds ligand size constant because the same drug is compared across two target sets.

Control construction:

- Exact match on **AlphaFold vs experimental** — pocket geometry differs systematically
  between the two and the panel is ~75% AF.
- Nearest-neighbour match on standardised `volume`, `depth`, `lid/hull` and
  `hydrophobicity`. Report achieved balance so match quality is visible.
- Controls must not be in **any** cluster.
- **Every pairwise Pmax within a control set below 0.5**, so "unrelated" is true by
  construction rather than assumed.
- Matched in number, per cluster; the statistic is computed per cluster.

Note the donor pool is ~178 proteins for ~31 controls, so exact matching on several
continuous variables will fail. Use nearest-neighbour on a few standardised variables or
coarse bins.

### Calibration run (in progress 2026-11-03)

`volume` and residue count are not the only things that make a pocket score well with
anything. `depth`, enclosure (`lid/hull`) and `hydrophobicity` are at least as important,
and which of them actually matter is an empirical question. The user ruled out reusing
2023 scores, so it is being measured on current data.

Design: **48 eligible net1 pockets** chosen by maximin (space-filling) design over
standardised `volume`, `depth`, `lid/hull`, `hydrophobicity`, `siteAtms`, balanced 24 AF /
24 experimental, against a **fixed 250-ligand e-Drug3D sample stratified by heavy-atom
count** (4-112 heavy atoms, median 24). 12,000 dockings, ~3 h on 64 cores. Run folder
`disease_subtypes_v2/Liver_cancer/net1_calib`, self-contained (only the chosen receptors
and pocket files copied in). Same smina settings, seed 42 and pdbfixer pH 7.0 prep as the
real runs, or the model will not transfer.

The selected 48 span the full eligible range on every descriptor (volume 310-1972,
depth 10.4-53.6, lid/hull 0.03-0.57, hydrophobicity 0.58-0.86, siteAtms 43-482).

Because the ligand sample is stratified by size, this one run calibrates **both** the
pocket-geometry term and the ligand-size term.

**The payoff is residualisation rather than matching.** If descriptors predict a pocket's
baseline generosity, test on

> `score - E[score | pocket geometry, ligand size]`

which uses all 209 panel proteins instead of discarding to ~31 matched controls, and
dissolves the small-donor-pool problem. Matching then becomes a robustness check and the
control arm still sits on top.

Caveats: fitting on net1 pockets and applying to net1 results is mildly circular — hold
out calibration pockets to check the fit, and re-validate on another subtype once its
pockets exist. The calibration pockets are themselves panel members; harmless for a
geometry model, but document it.

### Library-level dockability filter

e-Drug3D's inclusion limit is MW <= 2000, so it contains peptides and other molecules far
outside the validated domain of Vina-family scoring. Library profile: median 24 heavy
atoms / 352 Da / **5 rotatable bonds**; p99 91 / 1,324 Da / 35 rot; maximum 128 / 1,994 Da
/ **64 rot**. The extreme tail is GnRH analogues (afamelanotide 64 rot, ganirelix 58,
degarelix 53, cetrorelix 50, goserelin/histrelin/triptorelin/nafarelin/leuprolide 42-44),
cyclic peptide antibiotics (colistimethate 47, daptomycin 39, bacitracin 35), bleomycin
(41), sugammadex (a cyclodextrin, 1,994 Da / 40 rot) and long-chain surfactants
(colfosceril palmitate 40, polidocanol 37, nonoxynol-9 35).

**This was found the hard way.** The first calibration run stalled: the ligand sample was
stratified across the library's *full* heavy-atom range (4-112), which guaranteed the
largest molecules were included, and single smina calls on those chunks ran for over two
hours at 100% CPU while the small-ligand chunks finished in minutes.

**Filter (user, 2026-11-03): exclude `rotatable bonds > 15` OR `MW > 900`.** Excludes 123
of 2,083 (5.9%), keeping 1,960.

Notes for the builder:

- There is **no class flag to filter on**. The SDF carries only `status`, `name`, `cas`
  and `ID` — no MW, no molecule class. The filter must be computed. Pure-stdlib
  computation is possible (`docs/.../2026-11-03-panel-selection-spikes/ligsize.py`:
  element masses for MW; rotatable bonds as single, non-ring, non-terminal bonds between
  heavy atoms, with non-ring found by bridge detection on the heavy-atom bond graph).
  `obabel` is already in the image if preferred.
- The filter catches **flexibility, not peptide-ness**. It removes ~106 legitimate
  flexible oral drugs (aliskiren 20 rot / 552 Da, dabigatran etexilate 20 / 628,
  vilanterol 16 / 486) alongside the ~17 true peptides. That is correct — torsion count is
  what degrades Vina search and scoring — but the exclusion list must be reported.
- **No HCC or AD positive-set drug is affected.** The HCC agents top out at apatinib
  (40 heavy atoms) and are comparatively rigid.
- This is not only a runtime fix. These ligands would otherwise dominate any best-score
  ranking on atom count alone, which is the extreme tail of the ligand-size confound.

### Drug-side statistical fixes

In `Gene_Disease_Association/docking_screen_analysis.ipynb`:

1. **Rebuild the positive set.** Currently an exact string match on one `indications`
   value (`'hepatocellular carcinoma (HCC); renal cell carcinoma (RCC)'`) plus four
   hand-added names, one of which (RAMUCIRUMAB) is a monoclonal antibody that cannot be
   docked.

   Replace with: **current clinical-guideline systemic therapies (EASL / NCCN / ESMO),
   intersected with dockable small molecules actually present in eDrug3D, reporting both
   counts** — "N guideline drugs, of which N are dockable and present". Stating the
   intersection explicitly removes the need for any approval-date cut, which matters
   because eDrug3D is a fixed snapshot and HCC therapy moved repeatedly (sorafenib 2007,
   regorafenib 2017, lenvatinib 2018, cabozantinib 2019, then the immunotherapy
   combinations).

   Exclude non-dockable entities: antibodies (ramucirumab, and the atezolizumab +
   bevacizumab and durvalumab + tremelimumab combinations) and metal complexes
   (cisplatin, oxaliplatin). Record drugs absent from eDrug3D entirely (donafenib,
   tivantinib, brivanib). Report targeted agents and cytotoxics **separately** — the set
   is bimodal and the halves carry opposite size biases.

   **Alzheimer's is in scope** and gets the same treatment; its approved list is short
   enough to enumerate and verify by hand (donepezil, rivastigmine, galantamine,
   memantine; lecanemab and donanemab are antibodies and undockable).

   If the user supplies a curated list from their own literature work, it supersedes this.
2. **Adjust for ligand size and pocket geometry** using the calibration model above.
   Heavy-atom count is the dominant ligand term; add **rotatable bonds**, since smina
   penalises flexibility and the HCC kinase inhibitors are comparatively rigid.
3. **Handle multiple testing.** The `score_threshold` x `min_targets` grid over three
   subtypes is thousands of tests and its best p-value is not a p-value. The
   cluster-versus-control test is a cleaner design that may remove this problem rather
   than needing a correction for it; if the grid is kept, pre-register one combination
   **before** looking at new results, or correct across the grid and report the whole
   surface.

`min_targets` remains the cumulative-effect operator and is the right statistic — now
computed **within a cluster**, which is both more interpretable and a natural
multiple-testing unit.

### The primary-target positive control is dropped — and its absence helps

A natural positive control would be: for each approved HCC drug, does the set of pockets
it scores well against include **its known primary target**? The subtype gene lists do not
contain the targets of sorafenib, lenvatinib, regorafenib or cabozantinib, so the check
cannot be run within the seven-subtype structure. Decided (user, 2026-11-03): accept this
rather than build a purpose-made gene list.

This should be stated as a **strength, not a gap**. If a drug's known targets are absent
from the panel by construction, then any binding it shows there is necessarily away from
its primary mechanism — which is precisely the cumulative off-primary-target effect the
study hypothesises. The signal is cleaner for their absence.

What is genuinely lost is a basic check that the docking is meaningful at all. A free
substitute would have been re-docking co-crystallised ligands, but `lig_name` is empty for
all 11,599 selected pockets, so DS2 as configured did not record them. They could be
recovered from HETATM records in the already-downloaded structures if a sanity check is
wanted later. Not blocking.

### Reporting requirements

Per subtype: counts at each filter step so the funnel from pockets to panel is
reconstructible; every excluded protein with its reason in `failed.tsv`; and the panel's
AlphaFold share, which is the main methodological exposure and must not be buried.

**Sensitivity analyses cost no extra docking.** Both the experimental-structures-only
panel and TM policy B are *subsets* of the primary panel, so they are produced by
filtering the primary docking results. Note that under the final filter definition the
experimental-only subset is roughly 25% of 209, about **50 proteins** — the earlier "178"
figure came from a filter combination without the volume ceiling, segment cap or
subpocket rule. A sensitivity analysis on ~50 targets is correspondingly weak; say so.

This holds only because both are subsets. Any future sensitivity analysis defined by
*loosening* a filter, admitting proteins outside the panel, would need its own docking.

## Container design: split into two images

The current single `polyphar:1.0` image (1.53 GB, CPU-only) cannot carry the new docking
stage: Uni-Dock needs a CUDA runtime, POSEIDON has no GPU at all, and the GPU hosts we
will actually use (Colab, TRUBA) have **no Docker daemon** and must run the image through
Apptainer. **Decision: two images.**

### Why split rather than one fat image

- A CUDA runtime base adds 2-3 GB to stages that never touch a GPU, and every stage from
  `ids` to `select` runs on POSEIDON, which has no GPU.
- The GPU image must convert cleanly to Apptainer. Conversion downloads and flattens the
  whole image, so size is a direct cost on Colab's disk and session time.
- The two have different release cadences: stages 00-06 are stable, the docking stage is
  the part that is changing.
- Neither host can build the other's image.

The shared `polyphar` Python package is installed into both, from the same source, so
there is one implementation and two runtimes.

### `polyphar-core:<ver>` — CPU, stages 00-06

Base must be Debian/Ubuntu with **i386 multiarch enabled**: PocketMatch uses the upstream
prebuilt 32-bit binaries at commit 4648bff, because the 64-bit build overflows. That
constraint rules out a plain slim-python or Alpine base.

Contains: the `polyphar` package; PocketMatch 32-bit binaries plus `libc6:i386`;
openmm/pdbfixer, OpenBabel, RDKit, Biopython, pandas, networkx; and the client side of
the DoGSiteScorer host queue.

Explicitly **not** in the image:
- **DoGSiteScorer.** Host-licensed, device-locked; it runs on the host and is fed through
  the file queue in `04_pockets/queue/`. Never bake it in, never mount it in, never use a
  privileged container or the Docker socket.
- **Uni-Dock / CUDA.**

smina stays, marked reference-only: it is what the CPU baseline and the
smina-versus-Uni-Dock agreement figures were measured with, and keeping it preserves the
ability to reproduce them.

### `polyphar-dock:<ver>` — GPU, stage 07

Base: `nvidia/cuda:12.x-runtime-ubuntu22.04` — **runtime, not devel**, since nothing is
compiled at run time.

Contains: the `polyphar` package (dock stage and shared code); Uni-Dock 1.2.0 from
conda-forge plus `unidock-tools` from pip; RDKit for ligand conformer preparation;
OpenBabel for receptor PDBQT preparation.

**On Colab, do not use the image.** Colab Pro+ is the primary GPU host (user decision),
and its VMs are ephemeral — every session starts clean. Converting and flattening a 3-4 GB
image through Apptainer on each fresh VM would cost minutes of session time and disk every
go, whereas `micromamba create -c conda-forge unidock` plus `pip install unidock-tools`
took a **measured 31 seconds** in the 2026-11-03 validation. So the stage has three
runtime paths for the same code:

| host | mechanism | why |
|---|---|---|
| **Colab Pro+** (primary) | micromamba + a **pinned lockfile** | ephemeral VMs; 31 s beats a multi-GB pull per session |
| TRUBA | Apptainer `.sif` built from the image | persistent filesystem, so the flatten cost is paid once |
| RunPod / vast.ai | the Docker image directly | an instance *is* a container from your image |

Reproducibility therefore rests on a **pinned conda lockfile**, not on the image, and the
image is built from that same spec so all three paths resolve identical versions. The
dock manifest records which path was used along with the resolved versions. Without the
lockfile the Colab path would silently drift as conda-forge moves, which is the one thing
that would make replicate scores incomparable across sessions.

**Apptainer compatibility is a hard requirement for the TRUBA path, and it constrains the
build:**
- Apptainer mounts the image **read-only** and runs as the invoking user, so nothing may
  be written inside the image and no fixed UID may be assumed. The current Docker wrapper
  passes `--user $(id -u):$(id -g)`; that has no Apptainer equivalent and must not be
  depended on.
- Every writable path must live in the bind-mounted run folder. No `$HOME` writes, no
  caches under `/opt`, no temp files outside the mount.
- No host-specific paths baked in (notably not `/content`, which is Colab-only).
- Must work under `apptainer run --nv`, which is how GPU passthrough is granted.
- Conversion path is registry-based: `apptainer build polyphar-dock.sif docker://<registry>/polyphar-dock:<ver>`.
  `docker-daemon://` will not work on a host without a daemon.

### Orchestration across the two

1. POSEIDON runs `polyphar-core` for stages 00-06 through `bin/polyphar-docker`, with
   DoGSiteScorer on the host, exactly as now.
2. A new **dock bundle** step exports only what docking needs: the selected targets'
   receptors and pocket files, the prepared ligand library, and `final_targets.csv`.
   Measured at ~100 MB for a full panel against the 7.9 GB `02_structures` and 9.8 GB
   `04_pockets` of the full stage outputs.
3. The GPU host runs `polyphar-dock` over that bundle — Docker on RunPod or vast.ai
   (where an instance *is* a container from your image), Apptainer on Colab and TRUBA.
4. Results come back (~400 MB per subtype) and merge into the run folder.

This works because docking is compute-bound, not data-bound, which the measurements
confirmed. Per-chunk `.done` markers already make the stage resumable, so a bundle can be
worked through across several sessions on a session-limited host.

### Provenance the dock manifest must record

With 5 replicates and a reimplemented scoring evaluator, per-run provenance matters more
than it did:

image digest; Uni-Dock version; `unidock-tools` version; CUDA runtime version; **GPU
model** (throughput and possibly scores vary by device); `search_mode`; `--batch_size`;
**replicate index and its seed**; the prepared-library hash; and the per-score mean and
standard error across replicates.

### Deliberately out of scope for v1

gnina, as a Vina-scoring cross-check and CNN rescorer, is a third image (its static
binary alone is ~2 GB). Add it only if the rescoring arm goes ahead.

## Tests to run before the production screen

All five fit in one Colab session for a handful of compute units, reusing
`unidock_test.tgz`. The first is the only one that can still overturn a decision.

1. **Pose-in-pocket validation (required).** Dock across the full box-size range and
   measure each best pose's centroid against the cavity's own `_gps_` grid points. The
   concern driving the 30 A cap is that a large box can swallow a neighbouring groove, so
   a good score need not involve the selected pocket — which would sever the link between
   pocket similarity and binding that the whole study rests on. This has never been
   tested; on CPU it was unaffordable. The box cap is being kept regardless (user
   decision), so this test now serves to **confirm the cap is sufficient**, not to decide
   whether to have one, and to quantify how often poses stray at the sizes we do keep.
2. **Replicate-count curve.** Run 8-10 replicates on a few targets and plot SE against N.
   N=5 was chosen by decision, not measurement; this confirms it is enough near decision
   thresholds, where it matters.
3. **Offset constancy.** n=580 shows smina and Uni-Dock agree (rho=0.951) but not that
   the +0.326 kcal/mol offset is flat across the score range. If it drifts, translating
   thresholds is not a simple shift. ~5,400 pairs are available from the existing 7,200
   smina reference scores in about 15 minutes. **Re-run after conformer preparation**,
   since the current figures used unprepared ligands.
4. **`balance` versus `detail`.** `balance` was chosen on the authors' recommendation.
   Confirm `detail` does not materially reorder ligands before committing ~600 GPU-hours.
5. **Conformer preparation effect.** Compare scores before and after preparation on a
   fixed target set, so the methodological change is quantified rather than assumed
   harmless.

Grid spacing 0.375 against 0.5 is no longer worth testing: 0.5 existed only to rescue the
6 proteins that fail grid-legality, and the user accepted losing proteins to the box
constraints.

## Audit findings, 2026-11-03 — THIS DESIGN IS NOT READY TO IMPLEMENT

An independent audit of this document, the 19 spike scripts and the pipeline code found
problems serious enough that **implementation must not start until items 1-4 below are
resolved**. Some are decisions only the user can make. Corrections already applied inline
are marked where they appear; what follows is what remains.

### Blocking — these change what gets built or what the result means

**1. Clusters are defined over pockets; docking uses one pocket per protein, and they
need not be the same pocket.** The panel algorithm picks each protein's *highest
drugScore* eligible pocket; cluster detection runs over *all* eligible pockets; `dock`
consumes exactly one pocket per row. So a protein can belong to cluster 4 via pocket X
while the pocket actually docked is pocket Y, chosen by an unrelated criterion. The
cumulative-binding statistic would then measure binding to pockets that are **not** the
ones whose similarity defines the cluster — severing the very mechanism this study rests
on, by a different route than the oversized-box problem this document already worries
about. This is not hypothetical: cluster 1 is 6 proteins over 9 pockets, cluster 2 is 10
over 11. The existing note about counting proteins rather than pockets addresses
double-counting, not this. **Three options, all scientific decisions: dock the cluster's
pocket; dock every eligible pocket of a cluster member; or define clusters over the
one-per-protein panel pockets.**

**2. The panel number is inconsistent and nothing downstream was recomputed.** The engine
section says 183 (after the 30 A cap and grid-legality); roughly a dozen other places
still say 209 / 645 pockets, including cluster coverage (31/209), the control donor pool
(~178), residualisation ("uses all 209"), and the experimental-only sensitivity subset
("about 50"). **153 of 645 eligible pockets (24%) are not dockable, and nobody checked how
many of the 51 clustered pockets survive** — the three 3-protein clusters sit exactly on
the `>= 3 proteins` rule and could vanish. The donor pool falls to ~152, worsening a
matching problem already called marginal. Related: the box section says *clamp and flag*,
the engine section says *drop 26 proteins*. Those are different operations. Neither the
cap nor grid-legality appears in **Pocket eligibility**, so as written `select` emits 209
rows and `dock` silently loses 26 — which violates this document's own anti-goal against
silently dropping proteins, and Goal 2 (panel size should be a property of the filters,
not of a GPU grid constraint).

**3. The 7 clusters are not 7 independent tests, and they are proposed as the
multiple-testing unit.** The member lists hold 46 protein-slots over 31 distinct
proteins: CENPI appears in four clusters, EME1 in three, ten more in two. With one docked
pocket per protein, **the same score is reused in up to four cluster statistics**.
Clusters 1 and 3 share four proteins; cluster 7 is nearly a subset of 5 ∪ 3. Across 7
subtypes × ~7 clusters × 2 disease positive sets this is ~100 non-independent tests with
**no correction scheme proposed at all** — the existing multiple-testing item addresses
only the legacy threshold × `min_targets` grid.

**4. The calibration that residualisation depends on mismatches production on four axes.**
Wrong engine (smina, not Uni-Dock — and offset constancy is still an open test), wrong
ligand geometry (unprepared), wrong sampling design (the document describes
`build_calib.py`'s maximin 24 AF / 24 experimental design, which `build_calib2.py`
explicitly supersedes with a proportional box-volume-stratified sample — and the
discarded space-filling property is exactly what a geometry regression needs), and wrong
ligand range. With 48 pockets and ~6 collinear descriptors the geometry term is weakly
identified before any of this; no power or identifiability analysis exists. **Either
re-run the calibration against production conditions, or abandon residualisation for
matching and say so.**

### Serious — would draw reviewer objections

**5. The `MITOTIC` gene set has no provenance and was used to tune the filters it then
validated.** It is a hardcoded 23-gene list in seven spikes, found nowhere else in the
repo — no source, no citation, no derivation, most plausibly read off the hub
neighbourhood, i.e. from the data. Worse, `panel_final.py` and `tm_variants.py` print a
`mito` column for every threshold combination in the sweeps that chose
`volume <= 2000, segments <= 3` and TM policy C. So the thresholds were picked partly by
how many mitotic genes survived, and the ~3x enrichment is **not independent evidence**.
Document its provenance or drop every claim that depends on it.

**6. Rejecting unfiltered clustering was quantitative; accepting filtered clustering is a
narrative.** Unfiltered was rejected on purity 0.19, 5% density, internal Pmax 0.760,
across four resolutions. Filtered is accepted on internal Pmax 0.708-0.751 — essentially
the same — plus reading seven gene lists, two readings of which were wrong (see the
cluster section). `families.py` computes no purity, no enrichment, no base rate, no null
model, and tests one resolution. **The document never states the quantitative criterion
that filtered clustering passes and unfiltered fails.**

**7. There is no positive control of any kind, and the document spins that as a strength.**
Calling the absence of known drug targets "cleaner for their absence" is rhetoric: it
removes the only internal scale for what a genuine binder scores against these receptors
under this engine. The study would report cumulative sub-threshold binding on receptors
that are ~75% AlphaFold models with a new engine carrying an uncharacterised +0.33
offset, with nothing anchoring score to real binding. The redocking control (co-crystal
ligands) is one HETATM pass over already-downloaded mmCIFs and should be **promoted from
"later if wanted" to a prerequisite.**

**8. The 209/645 figures were computed under TM policy A, not the stated default C.**
Every spike producing panel and cluster numbers uses `nseg >= 2` (drop multi-pass only).
Only `tm_variants.py` implements C, and on the superseded 241-protein panel. Under the
stated default both numbers are wrong by an unmeasured amount, including the per-pocket
overlap rule's effect on cluster membership.

**9. The engine validation is quoted as settled and un-validated in the same document.**
"Validated 2026-11-03" and rho=0.951 sit beside the admission that those figures used
unprepared ligands and need re-measuring. The same applies to **0.176 s/ligand**, from
which the entire 614 GPU-hour / 1,191-unit budget derives: measured on 4 targets × 145
ligands of unstated rotatable-bond composition. Vina-family throughput is dominated by
torsion count — this document has a whole section on flexible ligands stalling smina — so
if that subset is rigid-skewed the budget is wrong in the dangerous direction. Stale
figures also remain: Constraints still says "dock calibrated at 57 s/ligand", which is
the superseded benchmark and the basis of the `dock(d)` columns that informed the filter
choice, making those ~3.5x optimistic. And the CPU "~61 days for seven subtypes" implies
~3.9 net1-equivalents while every GPU figure uses 7.0 × net1 — the two cost models are
not on the same basis.

**10. Residualising on pocket geometry may remove the signal it protects.** Cluster
membership is a deterministic function of pocket similarity, itself a function of shape
and size — the quantities being residualised out. The choice between residualisation and
matching is currently to be settled by how well descriptors predict baseline score, which
is the wrong criterion. Report the association between each descriptor and cluster
membership first.

**11. The matched-control construction has never been feasibility-checked.** Four
simultaneous constraints against a ~152-protein donor pool for 46 control slots. The
pairwise-Pmax constraint is probably satisfiable; the 4-variable matching is the binding
one, which this document half-concedes. Running the construction on the actual eligible
table takes minutes, and without it "report achieved balance" may mean reporting that the
match failed after the docking is paid for.

**12. Pmax is size-normalised, so clusters are partly size-matched by construction.**
`Pmax <= min(n1,n2)/max(n1,n2)`, so Pmax >= 0.7 forces the two pockets within ~30% of
each other in size — across a 300-2000 A^3 band, a 6.7x span. The claim that pocket-level
size filters "correct the same confound at source" was never tested on the filtered set;
the size-conditional threshold was only rejected on the *unfiltered* network. This also
cuts against the control construction, which matches on volume while requiring pairwise
Pmax < 0.5 — the two constraints pull against each other.

### Implementation blockers

**13. Replacing smina with Uni-Dock breaks the dock stage's resumability and provenance
design, and this document does not address it.** A chunk's `.done` job key includes the
smina version and assumes determinism; with a non-deterministic engine and N replicates
it must gain a replicate index, and a `.done` marker no longer certifies *which* scores.
`RESULT_COLUMNS` has no replicate or SE dimension.
`tests/test_stage_dock.py::test_cpu_and_chunk_invariance` **cannot pass** with Uni-Dock
and is not on the retirement list. `run`'s preflight checks `smina --version` and will
exit 2 without it. And **there is no execution model for a remote GPU**: the DoGSiteScorer
problem was solved with an explicit host queue and worker; the "dock bundle" here has no
queue protocol, no cross-host manifest reconciliation, and no story for how the Colab
polling loop drives `polyphar dock`. This is a stage-sized piece of design presented as a
bullet.

**14. "Grid legality at 0.375 A spacing" is never defined** yet costs proteins. No limit
cited, no enforcement point named. A builder cannot implement it from this document.

**15. The rotatable-bond definition is non-standard and unpinned.** `ligsize.py` counts
single, non-ring, non-terminal heavy-atom bonds via bridge detection and does **not**
exclude amide C-N bonds, as the conventional RDKit/Lipinski definition does. So its count
runs systematically high, and the offered `obabel` alternative *does* exclude amides —
swapping implementations silently changes which 123 ligands are excluded. Pin the exact
definition in Methods. The "~106 legitimate flexible drugs vs ~17 peptides" split is an
eyeball reading of the top 25, stated as a measurement.

**16. "PocketMatch is not re-run" contradicts "re-run from `similarity` onward".** The
former is a leftover from the spike context, where net1's pair file already existed; for
the six other subtypes `similarity` has not run at all, and running it *is* running
PocketMatch. It is cheap (~205k pairs against 67.2M), so there is no real problem — but
the document tells a builder two incompatible things. Also: `select` reads
`05_similarity/Pocketsim.csv`, not `PocketMatch_score.txt`, and `Pocketsim.csv` already
holds all surviving pairs with Pmin/Pmax and no threshold — so this document names the
wrong input file. And because the similarity-side filter deliberately excludes TM policy
while the panel filter includes it, the cluster graph contains pockets the panel rejects,
needing a second filtering pass that is not specified.

**17. The pipeline silently drops `Pmax == 1.0` pairs — the strongest possible cluster
edges.** `parse_scores` (`similarity.py:714`) excludes them as a "legacy filter". Since
same-accession pairs are excluded separately, this discards *cross-protein* pocket pairs
with perfect Pmax, exactly the tightest links a cluster could have. This document's entire
contribution is clusters at Pmax >= 0.7 and it inherits the filter unexamined. How many
such pairs exist among the 645 eligible pockets is a one-line check nobody ran.

**18. Colab feasibility at 614 GPU-hours is asserted, not demonstrated.** The operating
pattern was validated on a run costing 0.91 units. Production is 1,191 units over 614
GPU-hours on a client that hung three times out of three, under a ~24-hour session
ceiling — roughly 26 attended session-days. TRUBA is named as the fallback but never
sized: no partition, no queue wait, and no Apptainer build validated against any of the
four constraints this document calls hard requirements.

### Repair order

1. Reconcile the panel to one number end-to-end and decide whether the 30 A cap flags or
   drops. Re-derive clusters, coverage, donor pool and sensitivity subsets on the
   **dockable** set. Nothing else is meaningful until this is done.
2. Resolve the cluster/docked-pocket mismatch (finding 1). Scientific decision.
3. State the multiple-testing scheme for the cluster tests.
4. Re-run the calibration against production conditions, or drop residualisation.
5. Add the redocking positive control.
6. Document `MITOTIC`'s provenance or drop the enrichment claims.
7. Design the remote-GPU execution model properly.

## Open questions

All six questions put to the user on 2026-11-03 are answered and folded in above.
Two things remain outstanding, neither blocking a plan:

- **Residualisation or matching as the primary adjustment.** A technical choice, decided
  from the calibration results rather than by the user. Residualisation is preferred if
  the pocket descriptors predict baseline score well, because it uses all 209 panel
  proteins instead of discarding to ~31 matched controls; matching then becomes a
  robustness check. The calibration run finishes 2026-11-03 ~04:55 host time.
- **A curated HCC drug list**, if the user has one from their literature work. It would
  supersede the guideline-intersection approach above.

## What is left to the builder

- CLI surface: flag names and defaults for thresholds, TM policy, cluster threshold and
  the `dock` selector.
- Whether annotation is a subcommand or a cached lazy computation; the
  `pocket_annotation.csv` schema and cache-invalidation key.
- What to do with the now-meaningless `Similarity` and `Is_Hub` columns.
- `clusters.csv` schema.
- Test fixtures and cases. The superseded `test_reproduces_2023_targets` regression can
  go, since the 2023 comparison is a non-goal.
- Whether `similarity` filtering is on by default. It changes the pair count dramatically
  (645 eligible pockets on net1 gives ~205k pairs against 67.2M), so the ETA model and
  the `estimate` table need rechecking.
- Migration: re-run from `similarity` onward for every subtype. **`pockets` must not be
  re-run.**
