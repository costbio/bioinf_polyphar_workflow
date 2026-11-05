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
    export PATH=$HOME/.local/bin:$PATH      # the colab CLI lives here

Check the tools answer before starting a subtype:

    "$DOGSITE_BIN" --version >/dev/null && echo "dogsite ok"
    docker images polyphar        # polyphar:1.1 must exist (1.0 predates rdkit)
    colab usage                   # Colab auth + the compute-unit balance

If `colab usage` opens an authorization URL, the token has expired: visit it, approve, and
paste the code back. It expires often enough to catch you mid-week.

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
        -f docker/Dockerfile -t polyphar:1.2 -t polyphar:dev .

### 2026-11-06 (later): five ways the wrapper behaved differently from the tool

Running the seven subtypes through the REAL Colab/Uni-Dock path -- rather than a local
engine -- turned up five defects in `bin/polyphar-docker` and `bin/polyphar-colab`. All
five share one shape: **the host wrapper quietly doing something other than what
`polyphar` itself does**, and none was reachable by any test that existed before, because
those two scripts had no tests at all.

| what | why it mattered |
|---|---|
| `--library` was dropped at the round trip | `dock-export`/`dock-merge` fell back to the DEFAULT library. On a run folder that has one, the GPU would dock a **different library than the run planned**, with the wrong sha in every unit key. |
| `colab exec` exits 0 even when the remote script raised | A refused `launch.py` looked like a started worker: poll, `NO-PIDFILE`, merge **0 units**, report "done". On a real subtype that is 200 units you believe are running. |
| `import_bundle` copied non-atomically | An interrupted merge left a partial SDF that the next merge archived as a discarded draw. No draw lost, but the discard log fills with debris. |
| the drop-dedup was a substring scan | O(log size) per unit, and loose enough that two units could satisfy it between them without either having been logged. |
| the wrapper's `--to` defaulted to the LAST stage | `polyphar run` stops after `dock`; the wrapper ran `control` and `stats` **unasked** and died in stats' preflight. |

Both execs now require the script's own last line as proof -- `worker pid=... RUNNING`,
and a new `setup: ready` -- rather than the exit status. `tests/test_dock_defaults.py`
fails if the wrapper's default stage ever diverges from `run.DEFAULT_TO_STAGE` again.

**Running the seven subtypes.** They dock through real Colab Uni-Dock, which is the point:

    DOGSITE_BIN=$PWD/tests/fixtures/fake_dogsite2 POLYPHAR_VERSION=1.4 \
        pytest tests/test_subtypes_e2e.py

Three traps worth knowing, each of which cost a run:
- **`fake_dogsite2`, not `fake_dogsite`.** Same file; under its plain name it answers
  `--help` as DoGSite3, which conflicts with the default `--pocket-tool dogsite2`.
- **The commit must be pushed.** The VM pip-installs polyphar from it.
- **`ALLOW_SILENT=1`** unless `$POLYPHAR_HOME/pushover.env` exists -- `launch.py` refuses
  to start a worker that cannot notify you, which is right for a real run.

`POLYPHAR_COLAB_KEEP=1` holds one session across several run folders, so the engine
install is paid once instead of once per bundle. Useful for a night of subtypes too.

### 2026-11-06: results come home continuously, not once at the end

The Colab round trip used to pack and download **only after the worker finished**. A VM
reclaimed mid-session therefore cost every unit drawn since it started — up to 200 units
and ~20 GPU-hours on net1. Capping a session below 24 h was never a real defence.

`bin/polyphar-colab` now pulls on **every poll** (~2 min). `bin/colab/fetch.py` sends only
the finished units the host does not already hold, so a lost VM costs at most one poll
interval. `dock-merge` still runs once at the end: a pulled unit is already safe on disk,
and merging is far more expensive than downloading.

- Nothing is marked on the VM. The host sends what it has (`/content/have.txt`) and the
  selection is recomputed each time, so a dropped upload, exec, download or extraction is
  retried on the next poll and loses nothing.
- When the VM has nothing the host lacks, the poll costs no transfer at all.
- A unit still being written has no `.done` marker and does not travel.
- Interrupting the wrapper is now safe at any moment: re-running it re-attaches if the
  session is alive, and **merges what was already pulled** if the session is gone.
- `COLAB_FETCH_EVERY=N` fetches every Nth poll instead of every one.

### 2026-11-05 (later): image 1.3, unmodelled pockets and the PocketMatch hang

`export POLYPHAR_VERSION=1.3`. Two further defects, both found when AD Cluster1's
`similarity` appeared to hang:

- **PocketMatch does not fail on an input it cannot handle, it SPINS** — 100 % CPU, ~1 MB
  RSS, no output, no error. One pocket of `UNK` residues stalled 11 of Cluster1's 66 jobs
  for 2.5 h. A single UNK residue is enough. `similarity` now drops such a pocket before
  PocketMatch sees it, and `--pm-timeout` (default 3600 s) bounds every job.
- **Pockets lined by unmodelled backbone reached the panel.** A traced main chain with no
  side chains (`UNK`, or a poly-ALA stub) has no chemistry: docking it scores a backbone
  groove. `pockets` now drops a pocket whose lining is mostly unmodelled, so it reaches
  neither `select`, nor `similarity`, nor `dock`, and `structures` rejects a chain that is
  unmodelled throughout (> 95 %).

Measured over the 6,124 eligible pockets of Cluster1 and net1/net2/net3: 6,118 have no
unmodelled residue and 6 are 100 % unmodelled (all 6BQ1 chain A). Only Cluster1 is
affected — **net1, net2 and net3 need no pocket re-run.** Chains are a continuum, no chain
is entirely unmodelled, and the heavily traced ones still give fully modelled pockets, so
the chain guard is deliberately set where it discards nothing on this data.

To re-apply the pocket filter to a run folder whose DoGSiteScorer output already exists,
delete the stage manifest and re-run `pockets`: the queue re-uses every existing result,
so nothing is sent to DoGSiteScorer again.

    rm 04_pockets/manifest.json
    $POLYPHAR pockets

### 2026-11-05: image 1.2, and why `select` and `ligands` must be re-run

Three defects were found while investigating why no net1 target could complete. All three
are fixed in **image 1.2**, so `export POLYPHAR_VERSION=1.2`.

- `select` gave every target row of a protein the receptor of its FIRST row, while each
  row's box came from its own pocket. 83 of net1's 479 rows paired a pocket with a
  receptor from another structure, and 9 of the 37 targets docked had **no receptor atom
  inside the box at all** — every ligand scored in vacuum. Receptor and Chain are now
  resolved per row, and both `select` and `dock` refuse a box the receptor does not
  occupy.
- `ligands` now excludes records carrying elements the scoring function cannot type. The
  six boron drugs were silently given no poses by Uni-Dock, which alone made every target
  permanently incomplete. The library is 1,954 records, not 1,960.
- `uniprot` no longer crashes on a deleted accession (14 of AD Cluster1's 131).

**Any 06_targets/ or 07_docking/ built before 2026-11-05 is invalid.** Resume treats an
existing stage as done and the fingerprints do not carry the code version, so the stage
directory must be moved aside before re-running:

    mkdir -p superseded_mispaired_$(date +%Y%m%dT%H%M%S)
    mv 06_targets 07_docking superseded_mispaired_*/

Done already for HCC net1, net2 and net3: `select` and `ligands` were re-run under 1.2 and
every row now pairs a pocket with its own receptor. net1's 109 docked units were moved
aside with them; they were drawn against the old library and some against empty boxes.

## 1. The run folders

Inputs are in place (workbook + ligand library per folder) and all seven have `00_ids/`.
As of 2026-11-05: HCC net1/net2/net3 are complete through `select` and `ligands` under
image 1.2; AD Cluster1 is through `pockets` and needs `select` and `ligands` re-run under
1.2; AD Cluster2/3/4 have run `ids` only.

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

    $POLYPHAR dock             # GPU is elsewhere: exports, works it there, merges back
    $POLYPHAR dock --dry-run   # what is still outstanding
    $POLYPHAR control          # same round trip, then the readouts

`$POLYPHAR run` stops after `dock`: that is the default `--to`. It does **not** stop at
`dock` on POSEIDON — the wrapper splits the run around it, sending the docking to the GPU
host and carrying on afterwards. See below. The analysis stage is not
in the default range because it needs a positive set that is clinical-guideline reading
(below), and a run that would otherwise do nine stages must not refuse because the tenth
cannot start. A run that ends at `dock` without one prints and pushes a single line saying
so.

### Docking: POSEIDON has no GPU, so the work travels

`dock` defaults to **Uni-Dock on a GPU**. POSEIDON has none, so `$POLYPHAR dock` does not
run the engine locally: the wrapper exports the outstanding units, runs them on a Colab
Pro+ GPU, and merges the results back. You type the same command you would anywhere else.

    $POLYPHAR dock                      # one round trip
    $POLYPHAR dock --dry-run            # what is still outstanding
    $POLYPHAR run                       # splits itself around dock automatically

**Re-run until `--dry-run` reports nothing outstanding.** Every stage underneath resumes,
so each run exports only what is still missing, and a dropped session costs at most the
units in flight.

Sizing: one unit is one (target, replicate) — 1,960 ligands, about 6 GPU minutes. The
wrapper takes 200 units per session (~20 h) unless `POLYPHAR_MAX_UNITS` says otherwise. The
number is truncated **down to a whole number of targets**, because a target whose
replicates are split across sessions can be drawn on two different GPU models, and that
difference would land inside the replicate spread that exists to measure the engine's own
jitter.

`--engine smina` keeps everything local and on CPU, which is far slower (~264 h for net1)
but needs no GPU and no Colab.

What it refuses to do, and why:

- **A second bundle while an earlier one is unmerged.** Both would hold the same units,
  both would be worked, and the second merge's draws are discarded with the GPU time
  already spent. Instead, if the Colab session still holds that bundle (checked by the
  sha256 of its `run_key.json`), re-running **re-attaches**: no export, no upload, no new
  worker, just polling and then fetch, download and merge. That is the recovery when the
  host terminal dies mid-run: the detached worker keeps docking and notifying, and the
  same `$POLYPHAR dock` brings its results home. If the session is gone, it says so and
  stops.
- **Start with uncommitted changes under `polyphar/`.** The worker installs polyphar from
  the commit it pins, so the GPU would run different code than this host planned the unit
  keys with.
- **Start a worker that cannot notify you.** Without `~/.polyphar/pushover.env` the worker
  would run up to 20 hours in silence. `ALLOW_SILENT=1` overrides.

`dock-merge` exiting 3 is normal: this bundle merged and other units are still outstanding.

Underneath, `bin/polyphar-colab` drives the round trip and
[`bin/colab/README.md`](../bin/colab/README.md) documents the individual steps — useful for
another GPU host or for debugging a session, but you should not need either.

### The control's docking goes the same way

`control` redocks the co-crystal ligands, and it **must** use the same engine as `dock`:
its library-percentile readout scores the co-crystal ligand against the production library
at the same target, so a control run on smina beside a Uni-Dock `dock` compares two score
scales (Uni-Dock reads +0.326 kcal/mol less negative) and biases the percentile silently.

So `control` needs the GPU too. It is not a separate stage triple — the same commands take
`--control`:

    $POLYPHAR control        # labels targets, writes the library, docks it on the GPU
    $POLYPHAR control        # units now recorded -> computes the readouts

The first `control` stops with a message rather than failing: the library it writes is the
input the export needs. The last one finds every unit recorded and goes straight to the
readouts.

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

| Subtype | Workbook rows | `ids` → `ids.txt` lines | `uniprot` → `uniprot.csv` rows |
|---|---|---|---|
| net1 | 2060 | **459** (461 reviewed rows, two accessions appear twice) | 459 |
| net2 | 652 | 156 | 156 |
| net3 | 869 | 197 | 197 |
| Cluster1 | 78 | 131 | **78** |
| Cluster2 | 88 | 146 | **88** |
| Cluster3 | 132 | 217 | **132** |
| Cluster4 | 21 | 33 | **21** |

All seven were materialised on 2026-11-04 and match this table.

**The AD counts do not narrow at `ids`, and that is not a bug.** `ids` keeps *every*
accession in a multi-accession `UNIPROTID` cell and records those rows in `00_ids/rows.tsv`
(`row`, `gene`, `ids`, `solo`). The collapse to one protein per workbook row happens in
`uniprot`, via `collapse_rows`, which drops the losers into `01_uniprot/excluded.tsv`.
So Cluster1 legitimately reads 78 rows and writes 131 lines: 41 single-accession rows plus
37 multi-accession rows, which collapse back to 41 + 37 = 78. Seeing 131 after `ids` and
concluding the gene matching failed is a mistake — check `uniprot.csv`, not `ids.txt`.

**net1 is the whole point of the rerun.** The 2023 `net1UniprotID.txt` was the *net3* tab
(197 IDs). If `00_ids/ids.txt` comes out at 197, the wrong sheet was read — stop and fix it.
See `Gene_Disease_Association/NOTES_HCC_subtype1_rerun.md`.

AD Cluster3 should now contain **P85299 (PRR5)**. The 2023 list had a truncated `P8529`, so
2023 ran Cluster3 without it:

    grep -c P85299 00_ids/ids.txt      # expect 1

An AD row whose accessions cannot be resolved to one protein by gene name **stops the stage**
on purpose. Bring the row here rather than editing the workbook.

### Reference results: HCC net1, 2026-11-04

The first full v2 subtype, as a sanity check for the other six. Run on POSEIDON at
`polyphar:1.1` / commit `a0289cd`.

| Stage | Counts | Failed |
|---|---|---|
| `ids` | 459 accessions from 2060 rows | 0 |
| `uniprot` | 459 kept | 0 |
| `structures` | 3596 downloads | 32 |
| `fix` | 3132 chains | 0 |
| `pockets` | 3548 models, 88846 pockets, **11606 selected**, 711 outside_segment | 13 |
| `annotate` | 11606 annotated (3589 AlphaFold / 8017 experimental), 1127 grid_illegal | 0 |
| `select` | 2371 eligible pockets, **panel 269 proteins / 479 targets** (172 AF / 97 exp) | 0 |
| `similarity` | 2802528 pairs → 2754032 rows | 3 |
| `ligands` | 2083 records → **1954 kept**, 129 excluded | 0 |

Two things in that table are worth knowing before reading your own:

- **`grid_illegal=1127` at `annotate` is not a docking risk.** Those are pockets whose
  Uni-Dock grid at 0.375 Å exceeds the engine's limits (128/axis, 531441 total), and
  Uni-Dock *exits* on a violation. `panel.py` rejects them with reason `grid`, and the
  funnel confirms it: `tm 2590 → grid 2385`. Verified on this run — of 2371 eligible
  pockets and 479 final targets, **zero** are grid-illegal, by the stored flag and by
  recomputing the grid independently from the box dimensions. Only 205 are dropped at the
  `grid` step because the `volume` filter (11606 → 4682) already removed most oversized
  clefts first.
- **The 3 `similarity` failures are PocketMatch heap corruption**, not clean rejections:
  `cabbage-core` died with SIGABRT (`malloc(): corrupted top size`, and a glibc malloc
  assertion) on `P03952/6t7p_A` and two `P07288/2zch_P` pockets. Three of 2368 pockets does
  not move the network, but memory corruption in the comparison binary is worth
  understanding rather than filing as "3 failed".

## 4. Things that will bite

- **Check `docker ps` before rebuilding the image.** A running container keeps the image it
  started with, so a rebuild cannot disturb a live run — but a job started after it gets
  different code than one started before, and nothing in the run folder records which. On
  2026-11-04 a net2 docking run was found still on a pre-rebuild image, running smina after
  the default had changed to Uni-Dock.
- **Run the test suite inside the image before trusting a rebuild.** On the development
  machine 71 tests skip for missing openbabel and rdkit, and that is where two defects hid
  on 2026-11-04 (41 failures in `test_stage_dock.py` that the local run called green):

      docker run --rm --user 1000:1000 -v $PWD:/src -w /src --entrypoint python \
          polyphar:1.1 -m pytest tests -q -m 'not repo_files and not docker'

  `repo_files` and `docker` are deselected because the image has neither `bin/` nor a
  docker CLI; mounting the checkout makes those tests run and then fail on missing `pgrep`.
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
- **If a stage reports UniProt "not reachable", check the endpoint before the network.**
  UniProt's REST service has twice returned HTTP 500 for a request shape that used to work,
  while the service as a whole stayed up:
  - 2026-11-02: intermittent 500s from `/stream`, which is why `uniprot` uses
    `/uniprotkb/search` with pagination.
  - 2026-11-04: `/uniprotkb/search` began rejecting `accession:A OR accession:B ...` queries
    of **20 terms or more** (10 still returned 200, 20/40/80 did not), which killed `select`
    mid-run on net1. `tmfetch` now uses `/uniprotkb/accessions`, the endpoint meant for
    fetching a known list, which takes them comma-separated and is documented to 100 per call.

  The tell in both cases is that a free-text query (`?query=insulin`) returns 200 while the
  stage's own request 500s. That means the request shape, not the service — retrying will
  never succeed. `tests/test_tmfetch.py` pins the accessions endpoint so it cannot regress.

## 5. When something fails

`failed.tsv` names the item and the reason; the stage log is in `logs/`. The stage can be
rerun after the cause is fixed and will only redo what is missing. For a stage that stopped
with a message instead of a `failed.tsv` row (an ambiguous AD row, a pocket name PocketMatch
cannot keep apart), the message names every offending item.
