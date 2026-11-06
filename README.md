# polyphar: bioinformatics-based polypharmacology workflow

`polyphar` takes a list of UniProt accessions (one disease subtype) and runs:

    ids -> uniprot -> structures -> fix -> pockets -> annotate -> select -> similarity -> ligands -> dock

1. reviewed human UniProt records,
2. PDB chains (SIFTS-mapped) and AlphaFold models,
3. pdbfixer repair,
4. DoGSiteScorer 2 pockets (DoGSite3 as a sensitivity analysis),
5. [PocketMatch 2.0](https://github.com/1337deepesh/PocketMatch_v2.0) pocket-pair similarities,
6. per-pocket geometric annotation, then target selection by explicit pocket filters,
7. Uni-Dock docking of a ligand library (e-Drug3D) into the selected targets, on a GPU —
   `--engine smina` keeps it local and on CPU.

Everything runs through one CLI, `polyphar <stage>`, in a Docker image. Two things run
outside it: pocket detection, because the binary is licensed to the host (see
[Pocket detection](#pocket-detection-dogsitescorer--dogsite3)), and docking when this
machine has no GPU (see [Engine on another machine](#engine-on-another-machine)). Each
stage prints and notifies its work units and ETA before starting, and every stage resumes.

- [docs/RUNBOOK.md](docs/RUNBOOK.md) — running the seven v2 subtypes on POSEIDON.
- [docs/METHODS.md](docs/METHODS.md) — every difference from 2023 that could change a
  reported number, and the pre-registered statistics.

## First-time server setup

Once per machine, in order — each step checks what the next depends on. Needs `git`,
`docker` (your user in the `docker` group), `conda`/`mamba`, `tmux`, host `python3` >= 3.8.

    # 1. repo + host env (the DoGSite and regression tests run on the host)
    git clone https://github.com/costbio/bioinf_polyphar_workflow.git ~/repos/bioinf_polyphar_workflow
    REPO=~/repos/bioinf_polyphar_workflow && cd "$REPO"
    conda env create -f environment.yml && conda activate polyphar && pip install -e .

    # 2. DoGSiteScorer 2 on the host, licensed as the user who runs polyphar.
    #    --benchmark runs the preflight and stores the `pockets` calibration.
    #    Exit 3 = licence or installation problem; the binary's error is printed.
    export DOGSITE_BIN=/path/to/dogsitescorer-2.0.0/dogsite
    python3 "$REPO"/bin/polyphar-dogsite-worker --benchmark

    # 3. image. The tag is what POLYPHAR_VERSION selects; `dev` is also needed by
    #    tests/test_docker.py, which skips silently without it.
    docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) \
        -f docker/Dockerfile -t polyphar:1.0 -t polyphar:dev .
    export POLYPHAR_VERSION=1.0

    # 4. check the machine: PocketMatch's 32-bit binaries run in the image
    #    (`Exec format error` = no IA32 emulation, see Server checks), the wrapper runs as
    #    you with the real pocket binary, and RCSB/AlphaFold downloads work.
    pytest -W error tests/test_docker.py
    pytest -W error -m "dogsite and not regression" tests/test_stage_pockets.py
    pytest -W error tests/test_stage_structures.py -k end_to_end

    # 5. Pushover (format below), then a test notification
    mkdir -p ~/.polyphar && $EDITOR ~/.polyphar/pushover.env && chmod 600 ~/.polyphar/pushover.env
    set -a && . ~/.polyphar/pushover.env && set +a && pytest -m pushover tests/test_notify.py

    # 6. calibrate the container stages (pockets was done in step 2)
    mkdir -p ~/runs && cd ~/runs && NPROCS=64 "$REPO"/bin/polyphar-docker benchmark

    # 7. optional, ~minutes: end to end with the real binary, two proteins, 3 ligands.
    #    "no protein pair has Pmax >= 0.7" is a property of these two test proteins, not a
    #    fault -- rerun with POLYPHAR_E2E_PMAX=0.3.
    cd "$REPO" && pytest tests/test_run_e2e.py

**Reproducing the 2023 pockets** (DoGSiteScorer 2 only) tells you whether a newly licensed
binary behaves like the one used then. It reruns three net2 PLK1 chains and requires the
same pocket names with drugScore within 1e-4. If it fails, do not start a rerun until you
know why: pockets, and everything after them, would not be comparable.

    mkdir -p ~/adep_csb/disease_subtypes/Liver_cancer/net2/{pockets,pdbfixer/P53350}
    rsync -av workstation:/mnt/c/experiments/adep_csb/disease_subtypes/Liver_cancer/net2/pockets/DrugScore.csv \
        ~/adep_csb/disease_subtypes/Liver_cancer/net2/pockets/
    rsync -av 'workstation:/mnt/c/experiments/adep_csb/disease_subtypes/Liver_cancer/net2/pdbfixer/P53350/{2ogq_A,4o56_A,4h5x_B}.pdb' \
        ~/adep_csb/disease_subtypes/Liver_cancer/net2/pdbfixer/P53350/
    POLYPHAR_POCKET_TOOL=dogsite2 ADEP_DATA=~/adep_csb \
        pytest -W error -m "dogsite and regression" tests/test_stage_pockets.py

**For the DoGSite3 sensitivity analysis instead**: install it with the ZBH installer
(`/opt/ZBH/dogsite3/dogsite3`, `VERSION_INFO` beside it), activate once
(`dogsite3 --license <key>`), point `DOGSITE_BIN` at it, and add `--tool dogsite3` in step 2
and `POLYPHAR_POCKET_TOOL=dogsite3` for the runs. It needs its own run folder.

To run the seven v2 subtypes on POSEIDON stage by stage — exact commands, expected ID
counts, what to check — see **[docs/RUNBOOK.md](docs/RUNBOOK.md)**.

## Quick start

    docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) \
        -f docker/Dockerfile -t polyphar:dev .

One run folder per subtype, inputs in `inputs/`:

    REPO=~/repos/bioinf_polyphar_workflow
    mkdir -p ~/runs/net1/inputs
    cp clusters_features_ppi.xlsx e-Drug3D_2083.sdf ~/runs/net1/inputs/
    cd ~/runs/net1
    export POLYPHAR_VERSION=dev DOGSITE_BIN=/opt/ZBH/dogsite3/dogsite3 NPROCS=32
    "$REPO"/bin/polyphar-docker run --dry-run        # ETA table, does nothing
    "$REPO"/bin/polyphar-docker run --xlsx inputs/clusters_features_ppi.xlsx --sheet net1 \
        --id-column "uniprot Entry" --reviewed-column "Reviewed/unreviewed"

The HCC workbook has tabs `net1`/`net2`/`net3`. The Alzheimer workbook has no reviewed
column and needs `--gene-column`:

    "$REPO"/bin/polyphar-docker run --xlsx inputs/PPIN_ROSMAP_clusterFeatures_v7.xlsx \
        --sheet clusterN --id-column UNIPROTID --gene-column GeneID

**One workbook row is one protein.** Some AD cells list several accessions (the CALM3 row
holds CALM1, CALM2 and CALM3); `--gene-column` makes `uniprot` keep only the reviewed
accession whose primary gene name equals the row's gene, and stop with the candidates
listed if none or several match. HCC has one accession per row and keeps every reviewed one.

The wrapper takes the run folder from the current directory (or `$POLYPHAR_RUN_DIR`),
mounts it at `/work` and `~/.polyphar` (or `$POLYPHAR_HOME`) at `/polyphar_home`, runs as
your user and passes `NPROCS` through. Once `00_ids/ids.txt` exists, `run` needs no ids
options.

**Long runs** take hours to days. Keep the wrapper in the foreground under `tmux`. A tmux
session does not inherit variables exported before the server started, so pass them in the
command string (or `tmux new -e VAR=value`):

    tmux new -s net1 "POLYPHAR_VERSION=dev DOGSITE_BIN=/opt/ZBH/dogsite3/dogsite3 NPROCS=32 \"$REPO\"/bin/polyphar-docker run"

**Never `docker run -d`.** The wrapper must stay alive: it supervises the host pocket
worker and stops it when the container exits. To stop a run, stop the wrapper (Ctrl-C, or
`kill <wrapper pid>`). A raw `docker stop` allows only 10 s before SIGKILL, which can cut
the FAILED notification and run summary short — use `docker stop -t 60` if you must. An ssh
drop (SIGHUP) is handled like SIGTERM: the run ends with its notification, summary and pool
cleanup (exit 129).

### `polyphar run`

`run [--from STAGE] [--to STAGE]` runs the stages in order, or the contiguous range you give.
It accepts every stage's own options (`--xlsx/--txt/--sheet/--id-column/--reviewed-column/--gene-column`,
`--ph`, `--worker-timeout`, `--worker-dead-after`, `--pm-bin-dir`, `--volume-min`,
`--volume-max`, `--plddt-min`, `--tm-policy`, `--per-protein`, `--exclude`,
`--library`, `--seed`, `--replicates`, `--threads-per-job`,
`--unit-timeout`), plus `--max-failed-fraction` and `--keep-going`.

- Each stage goes through the same gate as its own subcommand. The gate computes the exact
  units, then prints and notifies the ETA before the first unit, then reports progress
  milestones and done or FAILED.
- `ids` runs when `--txt` or `--xlsx` is given. Without them, an existing `00_ids/ids.txt`
  counts as done, and a missing one is an error before anything starts.
- **Preflight.** Before the first stage, `run` checks what later stages need: with `dock` in
  range, the ligand library exists and `smina --version` works; with `similarity` in range,
  the PocketMatch binaries are there. A problem exits 2 before anything starts.
- **Resume.** Each stage's remaining units are computed first. A stage with 0 units left and
  a `manifest.json` whose recorded inputs (`inputs[].sha256`) all still match the current
  files is `skipped (done)`: it is not benchmarked, not run, and sends no notifications. Its
  manifest counts still feed the failure flags. If an input changed (for example a corrected,
  smaller `ids.txt` or `uniprot.csv`), the stage runs again and rebuilds its tables
  (`Protein_info.csv`, `DrugScore.csv`, ...) from the current inputs only, so removed
  proteins cannot reach similarity, selection or docking. A stage with 0 units but no
  manifest still runs, to finish its bookkeeping. A restarted run therefore shows ETAs only
  for the work that is left.
- **Code version.** Skipping a stage whose manifest was written by another code version
  (`git_commit`), or an unknown one, prints and notifies a warning (it does not block). An
  image built without `POLYPHAR_GIT_COMMIT` is warned about once at the start.
- Options of stages outside `--from`..`--to` get a one-line warning on stderr and are ignored.
- A stage that still has units left and no calibration on this machine is benchmarked first
  (`benchmark --stages <stage>`). `pockets` is the exception: the wrapper calibrates it on the
  host before `run` starts. If it is still uncalibrated, its ETA arrives at the 5 % milestone.
- The first failed stage stops the run with a priority-1 notification. `docker stop` (SIGTERM)
  and an ssh drop (SIGHUP) count as a failure at any point, including during an automatic
  benchmark.
- **Mass failures stop the run.** When `structures`, `fix` or `pockets` fails more than
  `--max-failed-fraction` (default 0.2) of its work items (for example an RCSB outage or a
  broken DoGSiteScorer license), the run stops with `run STOPPED at <stage>` (priority 1,
  exit 1) instead of docking a near-empty network for hours. `--keep-going` continues anyway.
  Expected exclusions are not failures: unmapped (accession, PDB) pairs and pockets dropped
  as `outside_segment`.
- A stage that finishes with failed units (`counts["failed"] > 0`), or `dock` with
  `incomplete_targets > 0`, does not stop the run. It is flagged in the final summary, which
  is then sent at priority 1.
- The summary (per stage: status, units, estimate, actual time, failed, incomplete targets)
  is printed, notified, and appended to `logs/run.log`. The notified text starts with the
  `failures:` line, because Pushover messages are cut at 1024 characters.
- `run --dry-run` prints the same table as `polyphar estimate` and does nothing else.

Exit codes of `polyphar run`:

| code | meaning |
|---|---|
| 0 | every stage finished, nothing flagged |
| 1 | a stage failed, or failed more than `--max-failed-fraction`, and the run stopped (143 when stopped by SIGTERM / `docker stop`, 129 by SIGHUP) |
| 2 | usage error (bad options, `--from` after `--to`, no ids source, failed preflight) |
| 3 | the run finished, but failed units or incomplete targets were flagged |
| 4 | the run stopped cleanly with work still outstanding: `--engine unidock` on a host with no GPU. Nothing failed; everything up to `dock` is done, and `dock-export` / `dock-worker` / `dock-merge` finish it |

## Run folder

One run folder = one subtype. Every stage reads and writes only inside it. Inputs are never
modified or deleted.

    <run>/inputs/                  ids source (xlsx or txt), ligand library sdf
    <run>/00_ids/ids.txt
    <run>/01_uniprot/uniprot.csv  not_found.txt  excluded.tsv  renamed.tsv
    <run>/02_structures/<UniProt>/{pdb}.cif {pdb}_{chain}.pdb AF-<acc>-F1.pdb
                       Protein_info.csv  sifts/  unmapped.tsv  annotation_mismatch.tsv
                       alphafold.tsv  rejected.tsv
    <run>/03_pdbfixer/<UniProt>/{pdb}_{chain}.pdb
    <run>/04_pockets/<UniProt>/...  DrugScore.csv  queue/
    <run>/05_similarity/Pocketsim.csv  PocketMatch_score.txt  cabbage/
    <run>/06_targets/final_targets.csv  centrality.csv  unique_pairs.csv
    <run>/07_docking/library/library_prepared.sdf  ligands.csv
    <run>/07_docking/<label>/replicates/r<i>.sdf  r<i>.log  r<i>.prov.json  r<i>.done
    <run>/07_docking/<label>/poses.sdf  (poses_r<i>.sdf with --poses all)
    <run>/07_docking/df_results_<label>.csv  discarded.tsv  discarded/  export/  obsolete/
    <run>/08_control/labels.csv  self_rank.csv  pose_in_pocket.csv  library_percentile.csv  summary.md
                     positive_control.json  docking/ (dock's tree, co-crystal ligands)
    <run>/08_stats/strata.json  m_statistic_<panel>.csv  results.csv  balance.csv
                   per_drug.csv  surface.csv  surface.md  summary.md
    <run>/logs/<stage>.log  run.log  notifications.log

Every stage directory gets `manifest.json` (written last; its presence marks the stage as
finished) and `failed.tsv` (`item, stage, reason`). Paths inside CSVs (`Pocket_Path`,
`Receptor_Path`, `File`) are relative to the run folder.

## One command per stage

All stages take `--run-dir` (default `/work` in the container, else the current directory),
`--nprocs` (default `$NPROCS`, else 4; results are identical for any value), `--dry-run`
(print units and ETA, then exit) and `--no-notify`.

    polyphar ids --xlsx inputs/clusters_features_ppi.xlsx --sheet net1 \
                 --id-column "uniprot Entry" --reviewed-column "Reviewed/unreviewed"
                 # AD: --xlsx inputs/PPIN_ROSMAP_clusterFeatures_v7.xlsx --sheet clusterN
                 #     --id-column UNIPROTID --gene-column GeneID   or: --txt inputs/ids.txt
    polyphar uniprot                                  # reviewed human entries only
    polyphar structures                               # SIFTS chains + AlphaFold F1 models
    polyphar fix [--ph 7.0]
    polyphar pockets [--pocket-tool dogsite2] [--drugscore-cutoff 0.5] [--min-volume 120] \
                     [--worker-timeout 600] [--worker-dead-after 21600]
                                                      # host worker runs DoGSiteScorer 2 / DoGSite3
    polyphar annotate                                 # per-pocket geometry, no DoGSiteScorer
    polyphar select [--volume-min 300] [--volume-max 1500] [--plddt-min 70] \
                    [--tm-policy C] [--per-protein 2] [--exclude O15350]
    polyphar similarity [--pm-bin-dir /opt/pocketmatch]   # AFTER select, on the panel
    polyphar ligands --library inputs/e-Drug3D_2083.sdf   # filter + MMFF94 conformers
    polyphar dock --library inputs/e-Drug3D_2083.sdf [--seed 42] [--replicates 3] \
                  [--threads-per-job 1] [--unit-timeout 43200]
                                                      # single host: plan, run, merge
    polyphar dock-export [--bundle DIR] [--max-units N]   # engine on another machine:
    polyphar dock-worker --bundle DIR                     #   export here, run there,
    polyphar dock-merge  --bundle DIR                     #   merge back here
    polyphar control [--control-replicates 1] [--engine smina]   # the positive control:
                                                      #   redock the co-crystal ligands
    polyphar stats [--positive-set inputs/positive_sets.csv] [--control-panel ../ad_run] \
                   [--permutations 10000] [--permutation-seed 42] [--min-coverage 0.98] \
                   [--positive-control inputs/positive_control.json]
    polyphar run [--from STAGE] [--to STAGE] [stage options...]
    polyphar estimate                                 # ETA table, next stage exact
    polyphar benchmark [--stages uniprot fix dock]    # per-machine calibration

Through the image: `bin/polyphar-docker <same arguments>`.

A single-stage subcommand that finishes with failed units (or `dock` with incomplete targets)
exits 3, like `run`, and sends its `done` notification at priority 1.

- **`polyphar estimate`** prints one row per remaining stage. The next stage is exact (from
  its inputs). Later stages are labelled `projected` and use values observed in this run
  where they exist, else priors from the 2023 HCC net2 run.
- **`polyphar benchmark`** times a small fixture workload per stage on one process and
  stores seconds per unit in `~/.polyphar/calibration.json`. Every finished stage also
  updates its entry. `pockets` is calibrated on the host by
  `bin/polyphar-dogsite-worker --benchmark --tool <tool>`, which the wrapper runs when there is
  no entry for the run's pocket tool (a calibration of the other tool is never used).

## Methods

[docs/METHODS.md](docs/METHODS.md) maps the 2023 scripts onto the subcommands and records
every difference that could change a reported number — ligand preparation and its
non-standard rotatable-bond definition, pocket-pair filtering, SIFTS chain selection, hub
ties, AlphaFold versions, seeds and protonation. Read it before writing Methods or quoting
a count.

## Docking units, replicates and resume

**One work unit is one (target, replicate): one engine invocation over the whole prepared
library.** `--replicates` (default 3) is how many independent draws each target gets;
replicate `i` runs with seed `--seed + i` and records it. The recorded seed is
**provenance, not a reproducibility claim** -- the GPU engine this stage is being moved to
scores 95-98% of ligands differently across identical same-seed repeats.

Because of that, a unit's identity includes its replicate index --
`(receptor sha, pocket sha, library sha, replicate, seed, search parameters, engine
version)` -- and **resume skips, it never redoes**. A `.done` marker certifies that the
draw was made and recorded, not what it scored; the result file
(`07_docking/<label>/replicates/r<i>.sdf`) is the only record of the scores, and its
marker `r<i>.done` sits beside it in the same directory, so a copy or a tar of
`replicates/` cannot separate a draw from its certificate. Raising `--replicates` adds
units and leaves the recorded ones alone; lowering it keeps the surplus records and
counts them as `units_present_outside_plan`.

The one legitimate repeat is a unit whose record is missing or truncated: it is discarded
and re-drawn under the **same** replicate index, and appended to
`07_docking/discarded.tsv` (append-only, across runs), because the discarded draw was a
real sample and dropping it is a small selection effect. Replacing a record whose key has
changed, or one left uncertified by an interrupted run, is logged the same way.

Two rules keep a discard from costing more than it says. **It is deferred**: what would be
thrown away is decided when the run is planned, but nothing is thrown away until the
replacement is actually about to be written, so a pass that is filtered, aborted, or does
not execute at all (`dock-merge`) costs no draws. And **the record is archived, not
deleted**: it moves to `07_docking/discarded/<timestamp>/<label>/r<i>.sdf` (with its log),
named in discarded.tsv's `archived` column, the same "moved, never deleted" policy orphan
labels follow.

A target counts as complete only when **every** replicate holds poses for every library
ligand. A record truncated on a clean `$$$$` boundary looks intact, so each replicate's
ligand coverage is checked, not only replicate 0's.

The receptor PDBQT the engine reads is deliberately **not** in the unit key (an Open Babel
upgrade must not invalidate every recorded draw). Instead, Open Babel's version is in
`07_docking/<label>/receptor.key`, so an upgrade re-prepares the receptor, and a PDBQT
whose sha256 differs from the previous manifest's is reported as a warning.

With `--threads-per-job 1` (the default), `NPROCS` units run at once and each does its 12
searches serially; `--threads-per-job 12` gives one unit all cores. `NPROCS` also sets how
many DoGSite3 / DoGSiteScorer jobs the host worker runs at once.

### Engine on another machine

`dock` needs the engine and the run folder on one host. When they are apart, the same work
splits in three around a bundle directory: `dock-export` plans here and writes the
outstanding units (replicate-major, capped by `--max-units` to one session's worth),
`dock-worker --bundle DIR` runs them there with no run folder, and `dock-merge --bundle
DIR` folds the results back in without needing an engine. A result is merged only when the
current plan has a unit with the same label, replicate **and key**, and a draw already on
record is never replaced.

`dock-merge` records the engine that actually drew the results (from the bundle's
`worker.json`) in the manifest as `engine_worker`, beside the `engine` the keys were
planned under, and warns up front when the bundle's options are not this merge's.

**`--max-units` does not coordinate between sessions.** Export only excludes units already
on record here, so two bundles exported before either is merged contain the same units. If
both are worked, the first merge wins and the second bundle's draw for that unit is
dropped -- logged in `discarded.tsv` with the bundle it came from, since dropping a real
sample is a selection effect. Merge a bundle before exporting the next one to avoid
spending engine time on draws that will be dropped.

## Pushover notifications

Create `~/.polyphar/pushover.env` (mode 600) with bare `KEY=value` lines and no quotes
(`docker --env-file` keeps quotes literally):

    PUSHOVER_APP_TOKEN=axxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    PUSHOVER_USER_KEY=uxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    chmod 600 ~/.polyphar/pushover.env

The wrapper passes the file to the container. Credentials are never in the image or the repo.
Notifications are optional and never fail a run. Titles start with `[<run folder name>]`.
These are sent:

- `<stage> starting`: units and ETA, before the first unit
- `<stage> 5% / 25% / 50% / 75%`: the 5 % one carries the revised ETA from the live rate
- `<stage> done`: time vs estimate, counts, failures
- `<stage> FAILED` (priority 1): the last log lines
- `pockets: no worker` (priority 1): no host worker heartbeat within `--worker-timeout`,
  repeated hourly; after `--worker-dead-after` (default 6 h) without any heartbeat the
  `pockets` stage fails (its jobs stay queued)
- `<stage> done` goes out at priority 1 when the stage has failed units or incomplete targets
- `run STOPPED at <stage>` (priority 1): too many failed work units (`--max-failed-fraction`)
- warnings: `run: code version unknown`, `<stage>: skipped stage from another code version`
- `run done` (priority 0), `run done with failures (N flagged)` (priority 1) and
  `run FAILED at <stage>` (priority 1): the `failures:` line, then the summary table
- a missing-credentials warning is logged once per process, not once per notification

Every notification, sent or not (for example without credentials), is also recorded in
`logs/notifications.log`.

## Pocket detection (DoGSiteScorer / DoGSite3)

DoGSiteScorer 2 (default, the primary tool) and DoGSite3 (the sensitivity analysis) are both
licensed for the host. A binary is never copied or mounted into the image and never committed
here.

- **Tool choice.** `--pocket-tool dogsite2|dogsite3` on `pockets` / `run` (default
  `$POLYPHAR_POCKET_TOOL`, else `dogsite2`), with `DOGSITE_BIN` pointing to that tool's
  binary. The wrapper reads the same option / variable and starts the host worker with
  `--tool <tool>`; the stage also writes the tool into every job file, and a worker refuses a
  job of the other tool (`failed`: `job is for X, this worker runs Y`).
- **The binary must be the chosen tool.** Before any job, the worker checks `DOGSITE_BIN`
  against `--tool` from what the binary says about itself: the `Tool:` line of `VERSION_INFO`
  next to it, else its `--help` text (DoGSite3 lists `--writeSiteResiduesPDB`, DoGSiteScorer
  2.0.0 `--druggability`). DoGSite3 must
  be confirmed; DoGSiteScorer 2 is refused only when the binary is DoGSite3. The preflight
  then also requires residue files named as the tool's are (`<out>_P_*_res.pdb` or
  `<out>_res_P_*.pdb`), not just exit 0 and `<out>_desc.txt`. Either failure exits 3 before
  a job is claimed.
- **One tool per run folder.** The tool is recorded in `04_pockets/TOOL` when the stage
  starts, in each `.done` and in the manifest (`versions.pocket_tool`,
  `versions.dogsite_version`). A run asking for the other tool than the folder's existing
  pocket outputs fails before anything is queued (`was made with --pocket-tool X, this run
  asks for Y`); queue records and manifests without a tool are DoGSiteScorer 2 outputs of
  earlier polyphar versions. Use a new run folder to switch tools.
- **Names (both tools).** `<out>` = `04_pockets/<UniProt>/<pdb>_<chain>` for a PDB chain and
  `04_pockets/<UniProt>/AF-<acc>-F1` for an AlphaFold model, so pocket files are e.g.
  `2mjw_B_P_1_res.pdb` (DoGSite3) or `2mjw_B_res_P_1.pdb` (DoGSiteScorer 2), and
  `AF-P25815-F1_P_1_2_res.pdb`. The UniProt folder already names the accession; only a chain
  extracted for several accessions (a chimeric chain) keeps the `<UniProt>_` prefix, because
  PocketMatch keeps base names only and they must stay unique. The longest realistic names
  (`AF-A0A0B4J2F0-F1_P_12_1_res.pdb`, 31 characters) fit PocketMatch's 32-character limit.
  The exception is a chimeric chain of a 10-character accession, which keeps the long
  prefix: `A0A0B4J2F0_7a8b_AB_P_12_1_res.pdb` has 33. `pockets` runs `similarity`'s name
  guard (≤ 32 characters, no collisions) on the selected pockets at its end and fails with
  the names listed in `failed.tsv` (no manifest), before `similarity` would.
- **DoGSiteScorer 2** (BioSolveIT; tested with 2.0.0), the **primary tool**: flags
  `-s -i -y -d -w 4`, outputs `<out>_desc.txt` / `<out>_res_P_<n>[_<m>].pdb`, selection
  `drugScore >= 0.5` (`--drugscore-cutoff`), which is the 2023 rule. Version 2.0.0
  reproduces the 2023 drugScore values exactly on the net2 regression chains (setup step 5),
  so a rerun is comparable with 2023.
  - **The base name is not 2023's.** 2023 wrote `<UniProt>_<pdb>_<chain>.pdb_res_P_<n>.pdb`.
    A current PDB entry with a multi-character author chain ID breaks PocketMatch's
    32-character limit that way (`P53350_9r1w_AAA.pdb_res_P_1_1.pdb` is 33), and PocketMatch
    keeps only base names, so truncation would merge distinct pockets. Decided 2026-11-02:
    both tools use the same base name, `<pdb>_<chain>` (see **Names** below), so the
    DoGSiteScorer 2 pocket file of that chain is `9r1w_AAA_res_P_1_1.pdb` (22). The
    accession is in the enclosing folder and in `DrugScore.csv`'s `UniProt ID` column, so
    nothing is lost; only file names differ from 2023.
- **DoGSite3** (ZBH, Universität Hamburg; tested with 1.3.0), the **sensitivity analysis**
  and not the primary tool: install with the ZBH installer (`/opt/ZBH/dogsite3/dogsite3`)
  and activate the license once with `dogsite3 --license <key>`. The worker runs
  `dogsite3 -p <model> -o <out> -s -d --writeSiteResiduesPDB -v 1`: subpockets on (`-s`, as
  with DoGSiteScorer 2), descriptors to `<out>_desc.txt` (`-d`), each (sub)pocket's complete
  residues to `<out>_P_<n>[_<m>]_res.pdb` (the PocketMatch input). The version is read from
  `VERSION_INFO` next to the binary (else `unknown`).
  - **Selection: `volume >= 120 Å³`** (`--min-volume`), the cutoff calibrated against the
    2023 drugScore selection (see below). DoGSite3 has no drugScore. Pockets and subpockets
    are judged alike, as drugScore was in 2023.
  - **Volumes are much smaller than DoGSiteScorer 2's.** On the same 20 pdbfixer chains of
    the 2023 net2 run, the largest DoGSite3 pocket had a median 0.24× the volume of the
    largest DoGSiteScorer 2 pocket (2OGQ_A, PLK1 kinase domain: 211 vs 546 Å³). Of 523
    DoGSite3 pockets in 50 random net2 models (40 chains, 10 AlphaFold models), only 11 (in
    10 models) have volume >= 300 Å³; DoGSiteScorer 2 kept 40 pockets by drugScore >= 0.5 on
    the 20 chains above. That is why the cutoff is calibrated rather than taken over. Check
    the selected-pocket counts (`counts.selected`) before similarity.
  - **Calibration on DoGSite3's own scale:** `scripts/calibrate_dogsite3_volume.py` runs a
    tool on the 2023 net2 inputs and finds the DoGSite3 volume cutoff that keeps the same
    fraction of pockets as drugScore >= 0.5 did in 2023 (pockets and subpockets on both
    sides), and compares the two selections by residue overlap. On 488 net2 models the kept
    fraction is f = 0.129 and the cutoff is **120 Å³** (the default). The two selections
    agree only loosely: precision 0.57 and recall 0.26 against the DoGSiteScorer 2
    selection, and only 54 % of the DoGSiteScorer 2 pockets have any DoGSite3 match. Treat
    DoGSite3 as a different method, not a drop-in substitute.
- **`DrugScore.csv`** (name kept from 2023) has `UniProt ID, Gene Name, Protein, Pocket`, the
  tool's descriptor columns, `Pocket Filename` (the residue file on disk), `Tool` and
  `Selected` (the tool's rule). The manifest records `selection_rule`. `similarity` uses the
  rows with `Selected` true and records that rule (`selection_rule` from the pockets
  manifest) in its manifest and fingerprint; a legacy `DrugScore.csv` without the column
  uses `drugScore >= 0.5` (recorded as `drugscore_cutoff`).
  Changing `--min-volume` / `--drugscore-cutoff` reruns `pockets` (selection only) on the
  next `run`.
- **Calibration** is per tool: `calibration.json`'s `pockets` entry records its `tool`; an
  entry of the other tool (no `tool` = DoGSiteScorer 2) leaves the stage uncalibrated, and the
  wrapper benchmarks the run's tool before `run`.
- **Host worker.** `bin/polyphar-dogsite-worker` needs only Python ≥ 3.8 and its standard
  library. For `pockets` and `run`, the wrapper starts it with `--watch` and waits until its
  preflight (one call of the tool on a bench chain) passes. It allows
  `POLYPHAR_PREFLIGHT_TIMEOUT` seconds (default 600) plus 60 s, and stops the worker on exit.
  A failing preflight points to a license or installation problem.
- **Timeouts and a dead worker.** A job running longer than `--job-timeout`
  (wrapper: `POLYPHAR_DOGSITE_JOB_TIMEOUT`, default 21600 s = 6 h) is killed with its process
  group and recorded as failed (`timeout after N s`). If the worker exits while the container
  runs and the queue is not `CLOSED`, the wrapper says so and stops the container, so the run
  ends with its FAILED notification instead of waiting for a dead worker. The worker logs to
  stdout and to `logs/dogsite-worker.log`; a broken stdout (closed terminal) never stops it.
- **Queue protocol** in `04_pockets/queue/`, with paths relative to the run folder:
  - The container writes `<id>.job` (`model`, `out`, `size`, `tool`).
  - A worker claims it by renaming it to `<id>.running.<host>.<pid>`.
  - The result is `<id>.done` (with `tool`, `tool_version`, `dogsite_sha256`) or
    `<id>.failed`.
  - The worker touches `WORKER.<host>.<pid>` as a heartbeat every 30 s (mtime set explicitly
    to the host's clock, so it also works on NFS).
  - The container writes `CLOSED` when every job has finished.
- **By hand** (for example, with polyphar in a host conda env instead of the image):

      bin/polyphar-dogsite-worker --run-dir ~/runs/net1 --dogsite "$DOGSITE_BIN" --nprocs 8 --watch
      bin/polyphar-dogsite-worker --benchmark --dogsite "$DOGSITE_BIN"   # calibrate pockets

  Add `--tool dogsite2` for DoGSiteScorer 2. Without `--watch` it runs the jobs queued now
  and exits.
- For tests and local end-to-end runs, `tests/fixtures/fake_dogsite` stands in for either
  binary (dogsite3 mode when `--writeSiteResiduesPDB` is given). It copies real outputs
  (three DoGSite3 1.3.0 output sets, or single net2 DoGSiteScorer 2 pockets), so its pocket
  geometry is not a prediction for the model it is given.

## Statistics (`polyphar stats`)

The last stage tests the study's claim — *approved HCC drugs bind more of the HCC panel, at
a useful level, than comparable molecules do* — in a comparison where ligand size and
pocket generosity cannot produce the result. The design is pre-registered in
`polyphar/data/preregistration.md` and written up in
[docs/METHODS.md](docs/METHODS.md#statistics-polyphar-stats); read it before changing a
parameter.

Three inputs it needs that nothing else does:

| Input | What it is |
|---|---|
| `inputs/positive_sets.csv` | **required, and hand-written.** The clinical-guideline (EASL / NCCN / ESMO) systemic-therapy list intersected with the dockable e-Drug3D records. Columns `set,ID,name,subgroup,source`. `polyphar stats` will not invent it and `run` stops in preflight without it. `tests/fixtures/positive_sets.csv` is a worked example. |
| `inputs/positive_control.json` | **produced by `polyphar control`** — never hand-written. Absent is **not** a pass: every result is then `unknown` and marked UNINTERPRETABLE. |
| `--control-panel <run folder>` | the pre-registered negative control, the same drug sets on the other disease's panel. Without it only one primary test runs and no contrast arm is produced. Pointing it at this run folder is refused. |

`control` redocks the co-crystal ligands and **must use the same engine as `dock`**: its
library-percentile readout scores the co-crystal ligand against the production library at
the same target, so mixing engines compares two score scales (Uni-Dock reads
+0.326 kcal/mol less negative than smina) and biases the percentile silently.

    polyphar control        # labels targets, writes the library, docks it
    polyphar control        # units recorded -> computes the readouts

The first call stops with a message rather than failing: the library it writes is the input
the export needs.

## Outputs and failures

- `failed.tsv` in each stage directory lists `item, stage, reason`. The manifest's `counts`
  has `failed`, the number of failed work units.
  - `structures`: unmapped (accession, PDB) pairs are listed in `failed.tsv` for visibility,
    but they are counted under `counts.unmapped`, not `failed`.
  - `pockets`: pockets dropped as outside the target's SIFTS segment are listed in
    `failed.tsv` and counted under `counts.outside_segment`, not `failed`.
  - `dock`: `failed` counts failed **units** ((target, replicate) pairs), `discarded` the
    draws thrown away and re-drawn in that run (see `discarded.tsv`).
    `incomplete_targets` counts targets with no
    results (a replicate not recorded, missing ligands, receptor preparation failure, or a
    missing receptor). An
    incomplete target gets no `df_results_<label>.csv`, and `failed.tsv` lists its
    missing ligand IDs. Treat `incomplete_targets > 0` as a failure.
- When `select` drops a target on a rerun, its docking outputs are moved to
  `07_docking/obsolete/<timestamp>/`, not deleted.
- `df_results_<label>.csv` has `Ligand, ID, Pocket, n, score_mean, score_sd, score_se,
  score_min, score_max`: one row per ligand, summarising that ligand's BEST score in each
  of the `n` replicates. `score_sd` is the sample sd (ddof=1) and `score_se = sd/sqrt(n)`;
  both are empty at `n == 1`. The raw per-replicate pose records stay in
  `<label>/replicates/r<i>.sdf`, so any other summary can be recomputed without re-docking.
  The retired `df_results_<label>_fixed.csv` (`Pose, SMINA_Score, Ligand, ID, Pocket`) was
  replicate 0 only: an engine-specific column name, and one number where there is now a
  distribution. A leftover `_fixed` file is swept into `obsolete/`.
  **Downstream code must key ligands on `ID`, not `Ligand`**: 35 different e-Drug3D drugs
  share the title `Molecule-1`.
- **Engines.** `--engine smina` (default) is the CPU single-host path; `--engine unidock`
  is Uni-Dock >= 1.2.0 on a CUDA GPU through `unidock_tools`. The engine and its search
  settings (`--search-mode`, default `balance` = exhaustiveness 384; `--batch-size`) are
  part of every unit key, so changing them supersedes the recorded draws instead of mixing
  two populations. Older Uni-Dock is refused outright: before 1.2.0 it silently corrupted
  ~20 % of results for boxes >= 30 A. The gate checks the binary that will ACTUALLY dock --
  `unidock_tools` spawns a hardcoded `unidock` from `$PATH` and never reads `$UNIDOCK_BIN`,
  so a `$UNIDOCK_BIN` that is not that same file is refused -- and it requires the
  `--version` banner to say `Uni-Dock` (a bare `1.2.5` could be AutoDock Vina) and to be a
  released build, not an `-rc`. Scores beyond +-100 kcal/mol are dropped as the known
  upstream failure (issue #205, reported against 1.2.0 itself) and listed in
  `<label>/replicates/r<i>.filtered.tsv`; 1.2.0 also fixed out-of-bounds poses carrying
  `FLT_MAX` affinities, which the same filter catches. A library whose records already
  carry a `<minimizedAffinity>` field is refused before any GPU time: `unidock_tools`
  passes input SD fields through, so that stale score would come back beside the real
  `<docking_score>`. `--engine unidock` on a host with no GPU fails immediately with the
  `dock-export` / `dock-worker` / `dock-merge` instructions, and `polyphar run` stops
  cleanly before `dock` there (exit 4) rather than failing.
- **ETA on a GPU.** Uni-Dock docks one unit at a time -- it batches the whole library onto
  the one device -- so `--nprocs` does not shorten its ETA and does not scale the
  seconds-per-ligand written back to `~/.polyphar/calibration.json`. The dock stage reports
  its own slot count for both.
- `--poses none|first|all` (default `first`) decides which replicates' poses are written
  out. Replicate 0's poses are what the pose-in-pocket validation reads.

## Reproducibility caveats

- **Resume treats existing outputs as done.** Outputs written by code before a fix (for
  example, unseeded pdbfixer) still count as done on resume. To get outputs of the current
  code, use a fresh run folder.
- **Stage skip checks do not include the code version.** `run` warns when it skips a stage
  made by another code version, but does not redo it. After upgrading polyphar, use a fresh
  run folder, or delete a stage's directory to recompute it.

## Server checks

- **PocketMatch** ships as upstream prebuilt 32-bit binaries and needs IA32 emulation in the
  host kernel. `Exec format error` from `similarity` means the kernel lacks it. Check on the
  server with `pytest tests/test_docker.py -k pocketmatch` (runs PocketMatch in the image).
- **AlphaFold / RCSB downloads**: `pytest tests/test_stage_structures.py -k end_to_end`
  downloads the structures of two proteins, including their AlphaFold models;
  `bin/polyphar-docker benchmark --stages structures` checks RCSB from inside the container.

## Tests

    conda env create -f environment.yml && conda activate polyphar && pip install -e .
    pytest -W error                    # unit tests; network tests run by default
    pytest tests/test_docker.py        # needs the polyphar:dev image
    POLYPHAR_E2E_PMAX=0.3 DOGSITE_BIN=$PWD/tests/fixtures/fake_dogsite \
        pytest tests/test_run_e2e.py   # full run through the image, ~3 min

The end-to-end test runs `polyphar run` on P25815 + P53350 with a 3-ligand library through
`bin/polyphar-docker`. It needs docker, the `polyphar:dev` image (or `$POLYPHAR_VERSION`) and
network access. With `fake_dogsite` the pockets are copies of a few real pockets, so the best
pair can stay below Pmax 0.7 (0.39 with the DoGSiteScorer 2 copies): set
`POLYPHAR_E2E_PMAX=0.3`, or `select` correctly stops the run with "no protein pair has
Pmax >= 0.7". With the real binary
(`DOGSITE_BIN=/opt/ZBH/dogsite3/dogsite3`), run it without `POLYPHAR_E2E_PMAX`.

Markers `dogsite`, `pocketmatch`, `docking`, `docker`, `regression` (`ADEP_DATA`) and
`pushover` skip unless their tool or environment is available.
