# Runbook: the v2 reruns, one stage at a time

Running the seven subtypes on POSEIDON, one stage at a time. `README.md` is the reference
for what the stages do and why; this file is the order, the commands, and what to check.

**DoGSiteScorer 2** (`drugScore >= 0.5`) is the pocket tool, so no
`--pocket-tool` argument is needed below. DoGSite3 with `--min-volume 120` is the
sensitivity analysis and needs its own run folder.

## 0. Once per shell

Work inside `tmux`: the wrapper runs the host DoGSiteScorer worker, and losing the ssh
connection mid-stage would kill it.

    tmux new -s polyphar          # or: tmux attach -t polyphar

    export REPO=/media/ilayda/Elements/adep_csb/repos/bioinf_polyphar_workflow
    export DATA=/media/ilayda/Elements/adep_csb/disease_subtypes_v2
    export POLYPHAR_VERSION=1.5
    export DOGSITE_BIN=/media/ilayda/Elements/adep_csb/dogsitescorer2/dogsitescorer-2.0.0/dogsite
    export NPROCS=64
    export POLYPHAR=$REPO/bin/polyphar-docker
    export PATH=$HOME/.local/bin:$PATH      # the colab CLI lives here

    "$DOGSITE_BIN" --version >/dev/null && echo "dogsite ok"
    docker images polyphar
    colab usage                   # Colab auth + compute-unit balance

`colab usage` opening an authorization URL means the token expired: visit it, approve,
paste the code back. It expires often enough to catch you mid-week.

**Every command below runs from inside the subtype's run folder** — the wrapper takes the
run folder from the current directory.

### Images

| Tag | Commit | Notes |
|---|---|---|
| `1.1` | `76dd25b` | the HCC net1 reference results below were produced here |
| `1.2` | `a41a3f0` | receptor/pocket pairing, undockable elements, deleted accessions |
| `1.3` | `a77dc99` | unmodelled pockets, the PocketMatch hang |
| `1.4` | — | **do not use:** rebuilt repeatedly with different content, so the tag identifies nothing |
| `1.5` | `7897ece` | current; the Colab round trip hardened |

Build with a **new** version number, plus the commit as its own tag so content stays
identifiable even if a version tag is ever moved:

    sha=$(git rev-parse --short HEAD)
    docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) -f docker/Dockerfile \
        -t polyphar:<next> -t polyphar:latest -t polyphar:dev -t "polyphar:$sha" .

    docker run --rm --entrypoint sh "polyphar:$POLYPHAR_VERSION" -c 'echo $POLYPHAR_GIT_COMMIT'
    git -C "$REPO" rev-parse HEAD        # these two must agree

## 0b. After pulling

**Rebuild the image after any pull touching `polyphar/`** — the image carries the code, so
a stale one runs old stages under a new commit's name. The host wrappers
(`bin/polyphar-docker`, `bin/polyphar-colab`) are read from disk instead, so never pull
while a round trip is in flight or the run mixes two versions of them.

**Stage fingerprints do not include the code version** — `run` only *warns* on a mismatch.
Any stage dir written by a module that changed is silently kept and is wrong:

    for d in 0*/ 0*/annotation; do
      [ -f "$d/manifest.json" ] && python3 -c "import json,sys;m=json.load(open('$d/manifest.json'));print('$d', (m.get('git_commit') or '?')[:7])"
    done

**Move** (do not delete) the dirs of changed stages aside:

    mkdir -p superseded_$(date +%Y%m%dT%H%M%S)
    mv 05_similarity 06_targets superseded_*/

## 1. The run folders

Inputs are in place (workbook + ligand library per folder) and all seven have `00_ids/`.

| Subtype | Through | Next |
|---|---|---|
| HCC net1, net2, net3 | `06_targets` | `dock` |
| AD Cluster1 | `06_targets` | `dock` |
| AD Cluster2, Cluster3, Cluster4 | `00_ids` | `uniprot` onwards |

**No subtype has any recorded docking draw** — 0 `.done` markers across all seven
`07_docking/`. net1's earlier 109 units were set aside as `superseded_mispaired_*` when
the receptor/pocket pairing defect was found.

| Subtype | Run folder | `ids` arguments |
|---|---|---|
| HCC net1 | `$DATA/Liver_cancer/net1` | `--xlsx inputs/clusters_features_ppi.xlsx --sheet net1 --id-column "uniprot Entry" --reviewed-column "Reviewed/unreviewed"` |
| HCC net2 | `$DATA/Liver_cancer/net2` | as net1, `--sheet net2` |
| HCC net3 | `$DATA/Liver_cancer/net3` | as net1, `--sheet net3` |
| AD Cluster1 | `$DATA/Alzhemier/Cluster1` | `--xlsx inputs/PPIN_ROSMAP_clusterFeatures_v7.xlsx --sheet cluster1 --id-column UNIPROTID --gene-column GeneID` |
| AD Cluster2 | `$DATA/Alzhemier/Cluster2` | as Cluster1, `--sheet cluster2` |
| AD Cluster3 | `$DATA/Alzhemier/Cluster3` | as Cluster1, `--sheet cluster3` |
| AD Cluster4 | `$DATA/Alzhemier/Cluster4` | as Cluster1, `--sheet cluster4` |

The AD sheets are lowercase (`cluster1`), the folders capitalised (`Cluster1`).
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

**Results come home at every poll, not once at the end.** Nothing is marked on the VM: the
host sends what it already holds and the selection is recomputed each time, so a dropped
upload, exec, download or extraction is retried on the next poll and loses nothing. A unit
still being written has no `.done` marker and does not travel. Interrupting the wrapper is
safe at any moment — re-running re-attaches if the session is alive and merges what was
already pulled if it is gone. A reclaimed VM costs one poll interval.

A lost websocket is not a refusal: `colab exec` exits 0 when the cell ran and the script
refused (definitive), non-zero when the CLI itself gave up. The wrapper retries `setup.py`
(idempotent) and never re-runs `launch.py` — a second launch is refused with *a dock-worker
is already running* and would leave a worker burning GPU units unwatched; it polls instead.

| Variable | Default | What it does |
|---|---|---|
| `POLYPHAR_MAX_UNITS` | 200 | units per session |
| `POLYPHAR_COLAB_KEEP` | unset | hold the session across run folders, so the engine install is paid once |
| `COLAB_POLL_SECONDS` | 120 | poll interval; the poll is also the keep-alive |
| `COLAB_SETUP_TIMEOUT` | 1800 | client wait for the install (a fresh VM takes ~9 min) |
| `COLAB_LAUNCH_TIMEOUT` | 900 | client wait for the worker to start |
| `COLAB_FETCH_EVERY` | 1 | fetch every Nth poll |
| `COLAB_MAX_POLL_FAILS` | 10 | consecutive poll failures before merging what is here |
| `POLYPHAR_ENGINE_VERSION` | `Uni-Dock v1.2.0` | the engine string the unit keys are planned under |

**`POLYPHAR_ENGINE_VERSION` must match what the VM's binary reports of itself**, or
`dock-merge` warns *drawn by X, planned under Y* on every merge and the one check guarding
against a silently changed engine becomes noise. It is pinned on the host because this
machine has no GPU to ask and a key cannot be assigned after the work it identifies.
**Changing it supersedes every draw made under the old string** — change it and the VM's
engine together, between runs, never during one. The gap this leaves, by construction: a
different engine on the VM is recorded in each unit's `prov.json` and warned about, but
does not supersede.

Underneath, `bin/polyphar-colab` drives the round trip and
[`bin/colab/README.md`](../bin/colab/README.md) documents the individual steps — useful for
another GPU host or for debugging a session, but you should not need either.

**The seven subtypes end-to-end**, through real Colab Uni-Dock, ~60 min:

    DOGSITE_BIN=$PWD/tests/fixtures/fake_dogsite2 POLYPHAR_VERSION=1.5 \
        pytest tests/test_subtypes_e2e.py

Three traps, each of which has cost a run: use **`fake_dogsite2`**, not `fake_dogsite`
(same file, but under its plain name it answers `--help` as DoGSite3 and conflicts with
`--pocket-tool dogsite2`); **the commit must be pushed**, since the VM pip-installs
polyphar from it; and **`ALLOW_SILENT=1`** unless `$POLYPHAR_HOME/pushover.env` exists.

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

**net1 and net3 are easy to confuse**: both sheets exist in the same workbook and net3 has
197 reviewed IDs. If net1's `00_ids/ids.txt` comes out at 197, the wrong sheet was read —
stop and fix it before anything downstream runs.

AD Cluster3 must contain **P85299 (PRR5)**; a truncated `P8529` in an input list is not a
UniProt entry and silently drops the protein:

    grep -c P85299 00_ids/ids.txt      # expect 1

An AD row whose accessions cannot be resolved to one protein by gene name **stops the stage**
on purpose. Bring the row here rather than editing the workbook.

### Reference results: HCC net1

The first full v2 subtype, as a sanity check for the other six. `polyphar:1.1` / `a0289cd`.

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

`grid_illegal=1127` at `annotate` is not a docking risk: those pockets' Uni-Dock grids
exceed the engine's limits and `select` rejects them, so zero of the 479 final targets are
grid-illegal. The 3 `similarity` failures are PocketMatch heap corruption (SIGABRT in
`cabbage-core`), not clean rejections — three of 2368 pockets does not move the network,
but it is memory corruption, not "3 failed".

## 4. Things that will bite

- **Check `docker ps` before rebuilding.** A running container keeps the image it started
  with, so a rebuild cannot disturb a live run — but a job started after it gets different
  code, and nothing in the run folder records which.
- **Run the suite inside the image before trusting a rebuild.** On a dev machine ~78 tests
  skip for missing openbabel and rdkit, which is where defects hide:

      docker run --rm --user 1000:1000 -v $PWD:/src -w /src --entrypoint python \
          polyphar:1.5 -m pytest tests -q -m 'not repo_files and not docker'

  `repo_files` and `docker` are deselected because the image has neither `bin/` nor a
  docker CLI.
- **`pockets` needs the DoGSiteScorer licence**, and this machine's evaluation activation
  expires **2026-11-16**. Get all seven subtypes through `pockets` inside that window;
  `dock` needs no DoGSiteScorer and has no deadline.
- **Docking dominates the total time.** Uni-Dock on an L4 is 0.176 s/ligand and flat across
  box volume — about 46 h per subtype per replicate, against ~34 days on CPU with smina.
  Run `select` for every subtype first to learn how many receptors there are.
- **`--nprocs` does not divide a Uni-Dock ETA.** The GPU runs one unit at a time. A GPU row
  shrinking with `--nprocs` is a bug, not a speedup.
- **Resume treats existing output as done** and the stage fingerprint does not include the
  code version. After any change, use a fresh run folder or move the stage dir aside.
- **One pocket tool per run folder.** A run asking for the other fails before queuing.
- **Pocket file names must stay under 32 characters**, PocketMatch's limit, and it keeps
  base names only — so a collision merges distinct pockets. Both tools use `<pdb>_<chain>`
  as the base name for that reason; `pockets` runs the name guard at its end and fails with
  the offending names listed.
- **UniProt "not reachable" usually means the request shape, not the network.** Its REST
  service has twice returned HTTP 500 for a shape that used to work while the service
  itself was up; the tell is that `?query=insulin` returns 200 while the stage's request
  500s. Retrying never succeeds. `tests/test_tmfetch.py` pins the endpoints in use.

## 5. When something fails

`failed.tsv` names the item and the reason; the stage log is in `logs/`. The stage can be
rerun after the cause is fixed and will only redo what is missing. For a stage that stopped
with a message instead of a `failed.tsv` row (an ambiguous AD row, a pocket name PocketMatch
cannot keep apart), the message names every offending item.
