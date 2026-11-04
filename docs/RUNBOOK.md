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
    export POLYPHAR_VERSION=1.1
    export DOGSITE_BIN=/media/ilayda/Elements/adep_csb/dogsitescorer2/dogsitescorer-2.0.0/dogsite
    export NPROCS=64
    export POLYPHAR=$REPO/bin/polyphar-docker

Check the tools answer before starting a subtype:

    "$DOGSITE_BIN" --version >/dev/null && echo "dogsite ok"
    docker images polyphar        # polyphar:1.1 must exist (1.0 predates rdkit)

**Every command below is run from inside the subtype's run folder** — the wrapper takes the
run folder from the current directory.

## 0b. After pulling a new polyphar

**Stage fingerprints do not include the code version** — `run` only *warns* on a mismatch, it
does not redo anything. So after a pull, any stage dir written by a module that changed is
silently kept and is wrong. Compare the recorded commit against the stages that moved:

    for d in 0*/ 0*/annotation; do
      [ -f "$d/manifest.json" ] && python3 -c "import json,sys;m=json.load(open('$d/manifest.json'));print('$d', (m.get('git_commit') or '?')[:7])"
    done

Then **move** (do not delete) the dirs of changed stages aside, so a mistake is recoverable:

    mkdir -p superseded_$(date +%Y%m%dT%H%M%S)
    mv 05_similarity 06_targets superseded_*/

`ids`, `uniprot`, `structures`, `fix` and **`pockets`** were unchanged by the 2026-11-04
update, so their output stays valid — which matters, because `pockets` is the licence-gated
one. `annotate`, `select`, `similarity` changed, and `ligands` / `dock` / `control` / `stats`
are new.

**The image must be rebuilt after a pull that touches `environment.yml`.** The 2026-11-04
update added `rdkit`, and `ligands` now fails the whole stage when rdkit is missing rather
than failing 2,083 records one at a time:

    docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) \
        -f docker/Dockerfile -t polyphar:1.1 -t polyphar:dev .

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
    $POLYPHAR control          # redocks the co-crystal ligands -> inputs/positive_control.json

`$POLYPHAR run` stops after `dock`: that is the default `--to`. The analysis stage is not
in the default range because it needs a positive set that is clinical-guideline reading
(below), and a run that would otherwise do nine stages must not refuse because the tenth
cannot start. A run that ends at `dock` without one prints and pushes a single line saying
so.

### Docking on a host with no GPU

POSEIDON has no GPU. With the default `--engine smina` everything runs in place and `dock`
simply takes a long time. With `--engine unidock`, `polyphar dock` **fails fast** rather than
starting, and `polyphar run` completes every stage up to `dock`, stops cleanly and **exits 4**
(`EXIT_INCOMPLETE`) — not 0, so a wrapper or `&&` chain cannot read it as a finished run.

To use the GPU engine, move the work to a GPU host and bring the results back:

    $POLYPHAR dock-export --bundle 07_docking/export/<name> --max-units 60
    #  copy the bundle to the GPU host, then there:
    #      polyphar dock-worker --bundle <bundle>
    #  copy it back, then here:
    $POLYPHAR dock-merge  --bundle 07_docking/export/<name>

Size each export to what one session can finish: a dropped session loses at most the units in
flight, and the next export recomputes what is outstanding. **Do not export a second bundle
before merging the first** — both would contain the same units, both would be worked, and the
second merge's draws are dropped (logged in `discarded.tsv`, but the GPU time is spent).

    # the analysis, once inputs/positive_sets.csv exists
    $POLYPHAR stats --control-panel $DATA/Alzhemier/Cluster3
    # or, in one go:
    $POLYPHAR run --to stats
    # or, after a finished dock:
    $POLYPHAR run --from stats

Each stage can be rerun safely: it skips what is already done. Add `--dry-run` to any stage
to see its units and ETA first.

### What `stats` needs that nothing else does

| File (in the run folder) | What it is |
|---|---|
| `inputs/positive_sets.csv` | **required.** The clinical-guideline (EASL / NCCN / ESMO) systemic-therapy list, intersected with the dockable e-Drug3D records. Columns `set,ID,name,subgroup,source`; `set` is `hcc` or `ad`; `subgroup` is `targeted` or `cytotoxic`. Building it is guideline reading, not code — `polyphar stats` will not invent it, and `polyphar run` with `stats` in range stops in preflight without it. `tests/fixtures/positive_sets.csv` is a worked example. |
| `inputs/positive_control.json` | **produced by `polyphar control`** — do not hand-write it. An absent one is **not** a pass: every primary, subgroup and contrast result is then reported `unknown` and marked UNINTERPRETABLE. Gated on redocking RMSD <= 2 A, the library top-decile fraction, and the median self-rank, each with a minimum n and a 98% panel-coverage floor. `pose_in_pocket_fraction` is still written but **gates nothing** — its null (~0.85 on real pockets) sat above the old 0.70 threshold, so it separated nothing. `control` needs `dock` to have run: without it the library-percentile readout is n=0 and the verdict is `unknown`, not a pass. |
| `--control-panel <other run folder>` | the pre-registered negative control: the same drug sets on the other disease's panel. Without it only one of the two primary tests runs, the contrast arm is not produced, and every output says so. Pointing it at this run folder is refused. |

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
| `control` | `08_control/` | `positive_control.json` — `redock_rmsd_success_fraction`, `library_top_decile_fraction`, `median_self_rank_percentile`, and `panel_scored_fraction` (>= 0.98). `n_redocked` is expected in the **tens**, not 687: only panel targets are labelled, ~13.5% of them carry a qualifying ligand. `summary.md` carries the caveat — read it before quoting any number |
| `stats` | `08_stats/` | `summary.md` first. `counts.positive_control` must be `yes` or the result is UNINTERPRETABLE; `counts.strata_edges_preregistered` must be `true`; `counts.positives_unresolved` and `positives_below_coverage` name every guideline drug dropped (in `failed.tsv`); `counts.balance_flagged` must be 0, or read the within-stratum size balance in `balance.csv`; `counts.contrast_tests` must be 1 or the pre-registered prediction was not tested |

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
- **Docking dominates the total time.** With smina, ~57 s per ligand per receptor on this
  machine and 1,960 prepared ligands, so roughly 32 CPU-hours per receptor before parallelism.
  With Uni-Dock on an L4 it is 0.176 s/ligand and flat across box volume — about 46 hours for
  one subtype at one replicate, against ~34 days on CPU. Run `select` for every subtype first
  to learn how many receptors there actually are.
- **`--nprocs` does not divide a Uni-Dock ETA.** The GPU runs one unit at a time; the estimate
  table knows this. If a GPU row ever shrinks with `--nprocs`, that is a bug, not a speedup.
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
