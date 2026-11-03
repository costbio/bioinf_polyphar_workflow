# Runbook: the v2 reruns, one stage at a time

For running the seven subtypes on POSEIDON by hand, one stage at a time, checking the output
of each before starting the next. `README.md` stays the reference for what the stages do and
why; this file is only the order, the exact commands and what to look at afterwards.

Decided 2026-11-02: **DoGSiteScorer 2** (`drugScore >= 0.5`, the 2023 rule) is the primary
pocket tool, so no `--pocket-tool` argument is needed anywhere below. DoGSite3 with
`--min-volume 120` is the sensitivity analysis and needs a **separate run folder**.

## 0. Once per shell

Run everything inside `tmux`: the wrapper starts and stops the host DoGSiteScorer worker, and
losing the ssh connection mid-stage would kill it.

    tmux new -s polyphar          # or: tmux attach -t polyphar

    export REPO=/media/ilayda/Elements/adep_csb/repos/bioinf_polyphar_workflow
    export DATA=/media/ilayda/Elements/adep_csb/disease_subtypes_v2
    export POLYPHAR_VERSION=1.0
    export DOGSITE_BIN=/media/ilayda/Elements/adep_csb/dogsitescorer2/dogsitescorer-2.0.0/dogsite
    export NPROCS=64
    export POLYPHAR=$REPO/bin/polyphar-docker

Check the tools answer before starting a subtype:

    "$DOGSITE_BIN" --version >/dev/null && echo "dogsite ok"
    docker images polyphar        # polyphar:1.0 must exist

**Every command below is run from inside the subtype's run folder** — the wrapper takes the
run folder from the current directory.

## 1. The run folders

Inputs are already in place (workbook + ligand library per folder).

| Subtype | Run folder | `ids` arguments |
|---|---|---|
| HCC net1 | `$DATA/Liver_cancer/net1` | `--xlsx inputs/clusters_features_ppi.xlsx --sheet net1 --id-column "uniprot Entry" --reviewed-column "Reviewed/unreviewed"` |
| HCC net2 | `$DATA/Liver_cancer/net2` | as net1, `--sheet net2` |
| HCC net3 | `$DATA/Liver_cancer/net3` | as net1, `--sheet net3` |
| AD Cluster1 | `$DATA/Alzhemier/Cluster1` | `--xlsx inputs/PPIN_ROSMAP_clusterFeatures_v7.xlsx --sheet cluster1 --id-column UNIPROTID --gene-column GeneID` |
| AD Cluster2 | `$DATA/Alzhemier/Cluster2` | as Cluster1, `--sheet cluster2` |
| AD Cluster3 | `$DATA/Alzhemier/Cluster3` | as Cluster1, `--sheet cluster3` |
| AD Cluster4 | `$DATA/Alzhemier/Cluster4` | as Cluster1, `--sheet cluster4` |

The AD sheets are lowercase (`cluster1`), the folders capitalised (`Cluster1`), as in 2023.
`Alzhemier` is misspelled on purpose.

AD uses `--gene-column GeneID` because many `UNIPROTID` cells hold several accessions; it
makes the `uniprot` stage keep one protein per workbook row, the one whose primary gene name
equals `GeneID`. Do **not** flatten the workbook into a txt list: that loses the row→gene
association, and e.g. Cluster1 would keep CALM1/CALM2/CALM3 instead of only CALM3.

## 2. The stages, in order

Start with **HCC net1**. Then the other subtypes, in whatever order suits you.

    cd $DATA/Liver_cancer/net1

    # see the whole plan and the ETA of every stage without doing anything
    $POLYPHAR run --dry-run

    $POLYPHAR ids --xlsx inputs/clusters_features_ppi.xlsx --sheet net1 \
                  --id-column "uniprot Entry" --reviewed-column "Reviewed/unreviewed"
    $POLYPHAR uniprot
    $POLYPHAR structures
    $POLYPHAR fix
    $POLYPHAR pockets          # starts the host DoGSiteScorer worker itself
    $POLYPHAR annotate
    $POLYPHAR select
    $POLYPHAR similarity       # AFTER select: it runs on the panel's eligible pockets
    $POLYPHAR ligands          # prepares inputs/e-Drug3D_2083.sdf once, for every target
    $POLYPHAR dock             # 1960 prepared ligands (2083 source records, 123 excluded)

Each stage can be rerun safely: it skips what is already done. Add `--dry-run` to any stage
to see its units and ETA first.

## 3. What to check after each stage

Every stage writes `manifest.json` **last** — if it exists, the stage finished. `failed.tsv`
lists what did not work (header only = nothing failed). The counts are in the manifest and
are also printed and pushed to Pushover.

    cat 0*/manifest.json | python3 -m json.tool | head -40
    wc -l 0*/failed.tsv

| Stage | Dir | Look at |
|---|---|---|
| `ids` | `00_ids/` | `ids.txt` line count; `rows.tsv` for AD multi-accession rows |
| `uniprot` | `01_uniprot/` | `uniprot.csv` rows; `not_found.txt`; `excluded.tsv` (unreviewed / non-human); `renamed.tsv` |
| `structures` | `02_structures/` | `Protein_info.csv` rows; `failed.tsv` for proteins with no structure; `versions` in the manifest |
| `fix` | `03_pdbfixer/` | one PDB per extracted chain; `failed.tsv` for pdbfixer failures |
| `pockets` | `04_pockets/` | `DrugScore.csv`: `counts.selected` (`drugScore >= 0.5`), `outside_segment` drops, `TOOL` = `dogsite2` |
| `annotate` | `04_pockets/annotation/` | per-pocket geometric annotation of the selected pockets |
| `select` | `06_targets/` | the target table; which protein(s) became the hub, and whether the top degree was tied |
| `similarity` | `05_similarity/` | `Pocketsim.csv` pair count; `PocketMatch_score.txt` |
| `ligands` | `07_docking/library/` | `counts.kept` / `excluded`; `stereo_inversions` and `implausible_formal_charge` must be explainable; `ligands.csv` |
| `dock` | `07_docking/` | per-target results; `incomplete_targets` must be 0; `library_prepared` must be `true` |

### Expected numbers at `ids` / `uniprot`

Check these before letting a subtype run on — a wrong count here invalidates everything after.

| Subtype | Workbook rows | Expected after `ids` |
|---|---|---|
| net1 | 2060 | **461 reviewed rows → 459 unique accessions** (two accessions appear twice) |
| net2 | 652 | 156 |
| net3 | 869 | 197 |
| Cluster1 | 78 | 78 rows |
| Cluster2 | 88 | 88 rows |
| Cluster3 | 132 | 132 rows |
| Cluster4 | 21 | 21 rows |

**net1 is the whole point of the rerun.** The 2023 `net1UniprotID.txt` was the *net3* tab
(197 IDs). If `00_ids/ids.txt` comes out at 197, the wrong sheet was read — stop and fix it.
See `Gene_Disease_Association/NOTES_HCC_subtype1_rerun.md`.

AD Cluster3 should now contain **P85299 (PRR5)**. The 2023 list had a truncated `P8529`, so
2023 ran Cluster3 without it:

    grep -c P85299 00_ids/ids.txt      # expect 1

An AD row whose accessions cannot be resolved to one protein by gene name **stops the stage**
on purpose. Bring the row here rather than editing the workbook.

## 4. Things that will bite

- **One stage at a time means the ETA is per stage.** `run --dry-run` prints the whole table.
- **`pockets` needs the DoGSiteScorer licence.** The current evaluation activation on this
  machine expires **2026-11-16**. Get the pocket stages of all seven subtypes done inside
  that window; `dock` needs no DoGSiteScorer and has no deadline, so it can wait.
- **Docking dominates the total time.** ~57 s per ligand per receptor on this machine, 2083
  ligands, so roughly 33 CPU-hours per receptor before parallelism. Run `select` for every
  subtype first to learn how many receptors there actually are.
- **Resume treats existing output as done, and the stage fingerprint does not include the
  code version.** After any change to polyphar, use a fresh run folder or delete the stage
  directory — otherwise half the output is from the old code.
- **One pocket tool per run folder.** A run asking for the other tool fails before queuing.
  DoGSite3 sensitivity analysis → its own run folder.
- **Pocket file names are not 2023's.** Both tools now use `<pdb>_<chain>` as the base name
  (`9r1w_AAA_res_P_1_1.pdb`), because the 2023 DoGSiteScorer 2 scheme exceeded PocketMatch's
  32-character limit on current PDB entries with multi-character chain IDs. Pocket IDs and
  drugScores are unchanged.
- **If `uniprot` reports "not reachable"**, check the endpoint before the network: UniProt's
  `/stream` endpoint was returning intermittent HTTP 500s on 2026-11-02, which is why the
  stage uses `/uniprotkb/search` with pagination.

## 5. When something fails

`failed.tsv` names the item and the reason; the stage log is in `logs/`. The stage can be
rerun after the cause is fixed and will only redo what is missing. For a stage that stopped
with a message instead of a `failed.tsv` row (an ambiguous AD row, a pocket name PocketMatch
cannot keep apart), the message names every offending item.
