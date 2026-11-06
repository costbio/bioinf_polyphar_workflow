# Methods: what polyphar does, and how it differs from 2023

Written for the manuscript's Methods section. Every difference from the 2023 scripts that
could change a reported number is here, with the reason. `README.md` is how to run the
thing; `docs/RUNBOOK.md` is the order to run it in.

## Old scripts → subcommands

| 2023 script | subcommand |
|---|---|
| `get_uniprot_data.py` | `polyphar uniprot` (with `polyphar ids` building the list) |
| `get_pdb_af.py` | `polyphar structures` |
| `fix_pdb.py` | `polyphar fix` |
| `get_pockets.py` | `polyphar pockets` (+ host `bin/polyphar-dogsite-worker`) |
| `get_pocketsimilarity.py` (README: `get_PMsimilarities.py`) | `polyphar similarity` |
| `Gene_Disease_Association/target_selection_analysis.py` | `polyphar select` |
| (none: the 2023 runs docked the raw library) | `polyphar ligands` |
| per-folder smina scripts (`smina_screen_*.py`) | `polyphar dock` |
| (none: the 2023 analysis had no size-matched test) | `polyphar stats` |

## Methods-relevant changes vs 2023

- **Ligand preparation (`polyphar ligands`, new).** The library is filtered and minimised
  once, before docking. **Dockability filter: `rotatable bonds > 15` OR `MW > 900` are
  excluded** — 123 of e-Drug3D's 2,083 records (5.9%), keeping 1,960; the library's limit
  is MW <= 2000, so it holds peptides outside Vina's domain (bacitracin, 1,421 Da, 35
  rotatable bonds). **The rotatable-bond definition is non-standard and must be reported as
  such:** single, non-ring (bridge-detected), non-terminal heavy-atom bonds, **amide C-N
  bonds included**, and the MW is a plain element-mass sum over the molfile. RDKit's and
  Open Babel's conventional (Lipinski) definitions exclude amides and would select a
  different library. Conformers are **minimised in place with MMFF94**, keeping the
  deposited geometry and stereochemistry, because Vina-family docking samples torsions but
  not ring pucker or stereo centres ("in place" means same molecule and same starting
  point, not small movement: the largest per-atom displacement has median 1.61 A and
  maximum 12.63 A). 23 records (metronidazole, nifedipine, chloramphenicol, ...) carry a
  nitro group whose nitrogen is left neutral in the molfile; their valence model is
  repaired (`[N+](=O)[O-]`) without moving an atom. 6 records — the boron-containing drugs
  (bortezomib, tavaborole, ixazomib, ixazomib citrate, crisaborole, vaborbactam), not metal
  complexes — have no MMFF94 parameters and keep their input geometry. MMFF94 has no
  chirality restraint, so every tetrahedral stereocentre's signed chirality volume is
  compared before and after minimisation and the minimised conformer is discarded if any
  flips: 1 record does (trabectedin, ID 1808). A C/N/O atom with |formal charge| > 2 is a
  malformed `M  CHG` line rather than chemistry (isosorbide dinitrate, ID 271, reads as a
  carbon at +15); such records are passed through flagged and counted, never silently.
  **Ligands keep their deposited protonation** (430 of the 1,960 prepared records carry a
  non-zero net formal charge) while receptors are protonated at pH 7.0 by pdbfixer — a
  deliberate asymmetry that has to be stated.
  `dock` then docks this prepared library and fingerprints it, so re-preparing re-docks,
  and records `library_prepared` in its manifest.

- **Pocket-pair filtering.** The 2023 script dropped every pair with `Pmin == 1` or
  `Pmax == 1`. Since PocketMatch normalises Pmin by the smaller pocket and Pmax by the
  larger, `Pmin >= Pmax` always, so that was effectively `Pmin != 1` and it discarded every
  *containment* pair. It is not kept: it cost two real hits above the `Pmax >= 0.7`
  threshold (CDC25C inside PLK1, 210 of 276 elements, Pmax 0.76, in net1 and net2), and no
  cross-protein `Pmax == 1` pair exists in any 2023 subtype. Row counts rise about 1% per
  subtype; the `Pmax >= 0.7` network is unchanged. Pairs sharing a **chain stem** are now
  dropped alongside same-accession pairs, because a chimeric chain SIFTS maps to two
  accessions is extracted into both accession folders and would otherwise self-compare at
  `Pmin == Pmax == 1`. Judge pairs on Pmax: `Pmin == 1` alone means containment, not
  identity.

- **Chain selection by SIFTS.** For each (accession, PDB entry), the chain is the one SIFTS
  maps to that accession (`SP_PRIMARY == accession`, strict, no ortholog fallback). Among
  matching chains, the one with the largest UniProt coverage wins, and ties go to the smallest
  chain ID. Entries with no such chain are listed in `unmapped.tsv` and skipped. The 2023 code
  took the first chain of every entry, so 15–22 % of its chains were a different protein.
  Chain IDs are author IDs (`auth_asym_id`). ProDy matches `chain=` against `label_asym_id`,
  so extraction selects by segment (ProDy stores `auth_asym_id` there).
- **cif annotation cross-check.** Each extracted chain's depositor annotation (`_struct_ref`
  UNP) is compared with SIFTS. A mismatch is a warning only and is listed in
  `annotation_mismatch.tsv`.
- **One protein per workbook row (AD).** A workbook cell with several accessions (AD
  `UNIPROTID`) is resolved to the single reviewed human accession whose primary gene name equals
  the row's `GeneID` (`--gene-column`); the others are listed in `excluded.tsv` as `same
  workbook row as ...`. The 2023 AD lists had one protein per row, and this rule reproduces
  them (it only affects the CALM3 row in cluster1). HCC has one accession per row and is
  unchanged.
- **Reviewed human entries only** (`Reviewed == reviewed`, organism 9606). Secondary
  accessions are renamed to their primary (`renamed.tsv`). Excluded entries are listed in
  `excluded.tsv`.
- **Current database versions**, recorded in the manifests: UniProt release, SIFTS
  Last-Modified, AlphaFold model versions (the latest, currently v6, located through the AlphaFold DB
  API instead of a hardcoded `v2` URL).
- **pdbfixer is seeded** (Reference platform), so its outputs are byte-reproducible. The
  legacy call sequence is kept, at pH 7.0.
- **Docking receptors are protonated once, by pdbfixer at pH 7.0.** A PDB receptor is the fix
  output, and an AlphaFold model gets hydrogens only (`addMissingHydrogens(7.0)`, seeded). The
  PDBQT (Open Babel, Gasteiger, rigid) keeps those hydrogens. The legacy Open Babel pH 7.4
  step (`CorrectForPH` + `addh`) is gone: it re-placed hydrogens at random on every
  conversion. Agreement with the 2023 pose-1 scores on 2mjw (20 ligands, today's smina):
  Spearman ρ 0.957 when the 2023 receptor PDBQT is reused, 0.994 for the new receptor path
  (mean |Δ| 0.070 kcal/mol). Exact 2023 scores cannot be reproduced (random hydrogen
  placement, older smina build). Receptor PDBQT sha256 values are recorded in the dock
  manifest.
- **Ligands**: smina `--addH 1` (its default, pinned explicitly), `--seed 42`. Kept from
  2023: exhaustiveness 12, num_modes 12, unpadded box over the pocket atoms.
- **Pocket detection with DoGSiteScorer 2** (default; see
  [Pocket detection](#pocket-detection-dogsitescorer--dogsite3)). This is the 2023 method and
  rule: `-s -i -y -d -w 4`, pockets and subpockets with drugScore >= 0.5. DoGSiteScorer
  2.0.0 reproduces the 2023 drugScore values exactly on the net2 regression chains, so the
  rerun's pockets are comparable with 2023's. What changes is **coverage**: the 2023 net2
  pocket stage has output for only about 495 of its 1301 models (~38 %), and 2.0.0 finds
  pockets for all sampled missing ones, so the 2023 similarity networks were built from a
  minority of the structures. Suggested Methods text: "Binding pockets and subpockets were
  detected with DoGSiteScorer (version 2.0.0) [citation] on each pdbfixer-repaired PDB chain
  and each AlphaFold model, with subpocket detection enabled. Pockets (and subpockets) with
  a drug score of at least 0.5 were retained for pocket similarity." DoGSite3 with a
  calibrated volume cutoff is the sensitivity analysis (`--pocket-tool dogsite3`).
- **Fusion / chimeric chains.** A chain SIFTS maps to several accessions (a fusion construct,
  e.g. a target fused to MBP or T4 lysozyme) is kept whole for pdbfixer and pocket detection.
  `Protein_info.csv` records the accession's SIFTS segments on that chain (`SIFTS segments`,
  `PDB_BEG:PDB_END` in author numbering) and the partners (`SIFTS other accessions`). A
  pocket with more than 50 % of its residues outside the segment(s) is dropped before
  `DrugScore.csv` (`failed.tsv`: `pocket outside target segment: <pct>% outside; partner(s)
  <accs>`; `counts.outside_segment`). Pocket residues are mapped back to author numbering
  through their heavy-atom coordinates in the extracted chain (pdbfixer renumbers residues
  but keeps the coordinates). When SIFTS leaves a bound empty or `None` (an unobserved
  end), the segment's SEQRES range (`RES_BEG..RES_END`) is mapped through the cif's
  `_pdbx_poly_seq_scheme` to the first and last observed author residue in it
  (`SIFTS segments source` = `seqres`, else `pdb`; counts `segments_from_seqres`). Only
  when that fails too (source empty, `segments_unverified`) are the chain's pockets kept
  and counted as `segment_unverified`.
- **Target selection**: no longer uses the similarity network at all. The 2023 rule took the
  top-degree node of the Pmax ≥ 0.7 graph plus its neighbours, which made the target count a
  property of how many structures were processed rather than of biology. The panel is now
  defined by explicit pocket filters — DoGSiteScorer's `drugScore`, volume 300–1500 Å³, a
  subpocket-parent rule, an AlphaFold pocket-pLDDT floor and a transmembrane policy — with
  **up to 2 pockets per protein** (picking one was measured to be the known ligand site only
  59% of the time). `similarity` now runs **after** `select`, on the panel's eligible pockets,
  and its output is descriptive: pocket similarity is applied post hoc to docking results, not
  used to choose targets. See `docs/superpowers/specs/2026-11-03-panel-selection-redesign-design.md`.


## Statistics (`polyphar stats`)

Pre-registered in `polyphar/data/preregistration.md`. Read it before changing a parameter.

The last stage tests the study's claim -- *approved HCC drugs bind more of the HCC panel,
at a useful level, than comparable molecules do* -- in a comparison where **ligand size and
pocket generosity cannot produce the result**. It is pre-registered in
`polyphar/data/preregistration.md` and designed in
`docs/superpowers/specs/2026-11-03-statistics-design.md`; read the pre-registration before
changing a parameter.

- **Within-target percentile.** Every target's ligands are ranked by `score_mean`
  (mid-ranks; more negative = better = percentile near 1), so every target's percentiles
  have mean 0.5. Pocket generosity is *removed*, not adjusted for.
- **The statistic** is `M = mean over PROTEINS of (max over that protein's targets of the
  percentile)` -- per protein, best pocket, so a protein contributing two pockets does not
  count twice.
- **Size matching** is quintile bins of heavy-atom count x tertile bins of rotatable-bond
  count, with the rotatable-bond count taken from `ligands.csv` (the pinned non-standard
  definition, amides included -- never RDKit's). The edges are computed from the library
  and written to `08_stats/strata.json` **before any docking score is read**, and compared
  against the pre-registered ones.
- **The test** is a one-sided within-stratum label permutation (10,000 draws, seeded):
  mean `M` over the positive set minus mean `M` over the matched controls.
- **The pre-registered prediction -- HCC-on-HCC > HCC-on-AD -- is tested, not inferred.**
  The two primary arms are each tested against *their own* controls, so reading the
  prediction off "one is significant and the other is not" is the difference-of-significance
  error. The `contrast` arm is a **paired within-stratum permutation on
  `M_primary - M_control`** over the ligands carrying an `M` on both panels. It carries no
  Holm correction. Note also that Holm makes the negative control *harder* to reject, which
  is backwards for an arm one wants to fail: both its raw and adjusted p are reported and
  the raw one is named as the honest figure.
- **Coverage.** A ligand needs `--min-coverage` (0.98) of the panel's proteins **and** of
  its targets. Both halves matter: at 0.90 an informatively-missing ligand's `M` is biased
  by up to +0.051, and `E[max of k uniforms] = k/(k+1)` means losing one pocket of a
  two-pocket protein moves `M` while protein coverage stays 1.000. Positives dropped for
  coverage are named in `failed.tsv` and `summary.md`.
- **Stratification removes the BETWEEN-stratum size confound, not the within-stratum one.**
  A within-stratum heavy-atom imbalance of ~0.6 sd can produce p < 0.05 from size alone.
  The test is not changed for that; `balance.csv` reports the positives' size against their
  controls' per stratum, the arm rows carry `balance_sd` / `balance_flag`, and a
  pre-registered `sensitivity` arm repeats the primary test on heavy-atom **deciles**.
- **The controls are already other approved drugs of the same size**: e-Drug3D is a
  library of approved and withdrawn drugs, not a generic screening set. That sentence is
  written into `results.csv`, `surface.md` and `summary.md` because it is load-bearing.
- **Replicate spread**: `M` is recomputed from each replicate's raw poses and reported per
  replicate (`stat_r0..2`). No analytic SE is derived for `M`, and `mean |delta|/sqrt(N)`
  is not one (`E|X-Y| = 1.128 sigma`).
- **Arms**: HCC drugs on this panel and HCC drugs on `--control-panel` are the two
  **primary** tests and are Holm-corrected; `contrast` is the prediction's own test;
  the targeted / cytotoxic split is a pre-planned subgroup, Holm within its own family;
  `sensitivity` is the deciles rerun; AD-drug arms are exploratory with `n` on the row.
  Per-drug elevation (`per_drug.csv`) is exploratory with Benjamini-Hochberg FDR over the
  resolvable subset.
- **`surface.csv` / `surface.md`** are the `score_threshold x min_targets` grid,
  descriptively, with **no p-values**. This is an **allowlist**, not a denylist: only the
  fields of `stats.SURFACE_COLUMNS` may be written, in the column list or smuggled through
  a row dict, and any other name is refused whether or not it looks like a p-value (a
  denylist let `p-value`, `adj.P.Val`, `Pr(>|t|)`, `pFDR` and a Cyrillic `\u0440_value`
  through). `surface.md` is guarded by the same check.
- **The positive-control gate.** `inputs/positive_control.json`
  (`median_self_rank_percentile`, `redock_rmsd_success_fraction`,
  `library_top_decile_fraction`, all fractions in `[0, 1]`, each with its own `n` and with
  `panel_targets_scored` / `panel_targets_total`) must show a median self-rank in the top
  5%, >= 50% of self pairs redocked within 2.0 A, >= 50% of co-crystal ligands in the top
  decile of the production library at their own target, and >= 98% of the panel scored --
  each over at least 10 observations, below which the verdict is `unknown` rather than
  `yes` or `no`. Failing or absent,
  every primary, subgroup and contrast row is marked uninterpretable -- the pre-registered
  consequence, enforced in code -- and `summary.md`'s arms table prints `(UNINTERPRETABLE)`
  beside each p. An out-of-range or non-numeric value is `unknown`, never `yes`.
- **The pre-registration is pinned by sha256** (`stats.PREREGISTRATION_SHA256`, recorded in
  every manifest): editing the document fails the test suite until the constant and an
  `## Amendments` entry are updated together.

### The positive control (`polyphar control`)

`stats` enforces a pre-registered gate: a FAILED positive control makes the primary result
UNINTERPRETABLE rather than reinterpreted, and an ABSENT one is not a pass. `polyphar
control` is what produces it. The receptors are ~63% AlphaFold models and the engine carries
an uncharacterised +0.33 kcal/mol offset from smina, so without this nothing anchors a score
to real binding; the usual anchor (the approved HCC drugs' own targets) is unavailable
because the subtype gene lists do not contain them.

It labels every panel target that is an experimental chain by the drug-like ligand sitting
in its pocket -- a HETATM residue of 12-60 heavy atoms, not solvent/ion/buffer/glycan/lipid,
with at least 70% of its heavy atoms within 4.5 A of the pocket lining; the criterion of
`docs/superpowers/specs/2026-11-03-panel-selection-spikes/rg_calibrate.py` (whose own 687
is a lower bound -- it derives the chain file by splitting the pocket name on `_res_P`,
which misses chimeric chains and every DoGSite3 pocket; `control` uses `Receptor_Path`
instead). One row per PHYSICAL ligand: a ligand filling a pocket and its nested subpocket
is counted once. Those ligands are extracted from the deposited
coordinates (Open Babel perceives the bonds; bond order and protonation may differ from the
depositor's chemistry) and redocked against the WHOLE panel through the same `dock`
machinery, writing to `08_control/docking/` so `07_docking` is untouched.

Three gated readouts, all in `08_control/positive_control.json` and installed at
`inputs/positive_control.json`, which is where `stats` reads them:

| readout | what it says | pass | min n |
|---|---|---|---|
| `median_self_rank_percentile` | where a ligand's own pocket ranks among the whole panel (chance 0.5) | <= 0.05 | 10 |
| `redock_rmsd_success_fraction` | share of self pairs whose best replicate-0 pose is within 2.0 A of the deposited coordinates (symmetry-insensitive, in place) | >= 0.50 | 10 |
| `library_top_decile_fraction` | share of co-crystal ligands in the top decile of the PRODUCTION library at their own target -- the axis `M` is built on | >= 0.50 | 10 |
| `panel_targets_scored / panel_targets_total` | how much of the panel carried a score at all | >= 0.98 | - |

`pose_in_pocket_fraction` is still reported, **descriptively, and gates nothing**: the
docking box is the pocket's own unpadded lining bounding box and the engine keeps the
ligand inside it, so random non-clashing centroids in that box pass ~79% against ~85% for
true deposited poses -- the old `>= 0.70` gate sat below its own null. RMSD <= 2 A has a
null of 0.0% on the same random placements.

Expect about **24** labelled ligands on a 482-target panel, not 687: only panel pockets can
have a self-rank, and the measured label rate on panel-style pockets is 13.5%.

The JSON also carries `source`, each readout's `n` and the caveat, because the numbers get
quoted by people who did not read the code:

> **Sanity check, not an optimisation target.** The home pocket is the **holo** conformation,
> induced-fit around this very ligand (DoGSiteScorer runs on the pdbfixer structure, i.e. the
> holo backbone with heterogens stripped), so this is the easiest benchmark in the field and
> says nothing about cross-docking or about the ~63% of the panel that is AlphaFold. Absence
> of a co-crystal ligand does NOT mean a pocket is undruggable -- most such structures were
> solved apo. The label is biased toward LARGE pockets and toward WELL-STUDIED proteins, and
> "drug-like" overstates it: across 1,142 net2 chains the commonest labelled codes include
> PSU 49 (pseudouridine), 7PE 22 (a PEG oligomer), NAP 20, FAD 18, ATP 12. And the self-rank
> is invariant to any monotone transform applied uniformly across targets, **including the
> +0.33 kcal/mol smina offset, which cancels exactly** -- so it cannot detect the defect its
> own rationale names.

A panel with no labelled pocket writes no JSON at all, and a rerun that can no longer produce
one deletes the stale file: `stats` must say "unknown: not run" rather than read a pass that
is no longer true.

### The positive set is an input, not code

Building it is clinical-guideline reading (EASL / NCCN / ESMO), so the stage requires
`inputs/positive_sets.csv` and fails with instructions when it is absent (`polyphar run`
stops in the preflight, before anything starts):

    set,ID,name,subgroup,source
    hcc,1234,SORAFENIB,targeted,NCCN 2026 / EASL 2025
    hcc,,LENVATINIB,targeted,NCCN 2026
    ad,2001,DONEPEZIL,,guideline AD therapy

`set` is `hcc` or `ad`; `ID` is the e-Drug3D library ID, or leave it empty and the `name`
is matched case-insensitively against the library; `subgroup` is `targeted` or `cytotoxic`;
`source` is free text. A row that does not resolve to a **dockable** library record is
never dropped silently: it is named in `failed.tsv` and in `summary.md`, because the
Methods must say which guideline drug was lost to undockability. The same is true of a
drug that resolves but does not clear the coverage rule.

`polyphar run` stops at **`dock`** by default for the same reason: the positive set will be
absent for most of the project's life, and the everyday command must not refuse the nine
stages it can do. `--to stats` includes the analysis, `--from stats` runs it afterwards,
and a run that ends at `dock` without a positive set prints and notifies one line saying
so.
`tests/fixtures/positive_sets.csv` is a worked example.

