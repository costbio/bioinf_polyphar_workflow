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

These are structurally coherent in a way the unfiltered communities never were. CDC6 and
ORC1 are both AAA+ ATPases of the pre-replication complex with homologous
nucleotide-binding folds, and CDT1 sits in the adjacent cluster. The fourth cluster holds
NDC80, SPC24 and SPC25 — three components of the same Ndc80 kinetochore complex with
related RWD domains. The homeobox proteins cluster separately, as DNA-binding grooves
should.

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

**This is the largest implementation risk in the design.** pdbfixer renumbers 1..N *and
inserts modelled-in residues*, so there is no constant offset from pdbfixer numbering back
to author numbering, let alone to UniProt. `Protein_info.csv`'s `SIFTS segments` gives
ranges (e.g. `30:175`), not a residue-level map. The annotation pass must therefore
reconstruct the mapping by aligning the pdbfixer output back to its input chain and then
through SIFTS — which may require retaining the raw SIFTS residue mapping at the
`structures` stage, or re-fetching it. Note that re-fetching must not disturb
`uniprot.csv` or `Protein_info.csv`, for the reasons in Constraints. Plan this task
first; the rest of the annotation pass depends on it.

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
- **The selector must enter `_stage_fingerprint`'s `opts`** (`dock.py:625`). It currently
  hashes `final_targets_sha256`, which does not change between selections, so resume
  would otherwise treat a different selection as already done.

### Docking cost

| scope | proteins | per subtype |
|---|---|---|
| all clusters (Pmax >= 0.7) | 31 | 0.67 d |
| clusters + matched controls | ~62 | ~1.3 d |
| whole panel | 209 | 4.5 d |

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
