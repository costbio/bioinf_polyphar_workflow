# polyphar: bioinformatics-based polypharmacology workflow

`polyphar` takes a list of UniProt accessions (one disease subtype) and runs:

    ids -> uniprot -> structures -> fix -> pockets -> similarity -> select -> dock

1. reviewed human UniProt records,
2. PDB chains (SIFTS-mapped) and AlphaFold models,
3. pdbfixer repair,
4. DoGSite3 pockets (DoGSiteScorer 2 as an alternative),
5. [PocketMatch 2.0](https://github.com/1337deepesh/PocketMatch_v2.0) pocket-pair similarities,
6. target selection on the similarity network (hub + neighbours),
7. smina docking of a ligand library (e-Drug3D) into the selected targets.

Everything runs through one CLI, `polyphar <stage>`, in a Docker image. Pocket detection is the
exception: DoGSiteScorer 2 (or DoGSite3) runs on the host (see
[Pocket detection](#pocket-detection-dogsitescorer--dogsite3)). Each stage prints and
notifies its exact work units and ETA before it starts, and every stage can be resumed.

## First-time server setup

Do these once per machine, in this order. Each step checks something the next one depends on.
You need `git`, `docker` (your user in the `docker` group), `conda` or `mamba`, `tmux` and a
system `python3` ≥ 3.8 on the host.

1. **Clone the repo and create the host env.** The DoGSite3 tests and the regression
   checks run on the host.

       git clone https://github.com/costbio/bioinf_polyphar_workflow.git ~/repos/bioinf_polyphar_workflow
       REPO=~/repos/bioinf_polyphar_workflow && cd "$REPO"
       conda env create -f environment.yml && conda activate polyphar && pip install -e .

2. **Install and license DoGSiteScorer 2 on the host** (the primary tool), then check that
   it runs. Install its license as BioSolveIT documents, as the user who runs polyphar.
   `--benchmark` first runs the worker's preflight, which is one DoGSiteScorer 2 call on a
   bench chain. It then stores the `pockets` calibration (for this tool) in
   `~/.polyphar/calibration.json`. Exit code 3 means a license or installation problem, and
   the binary's error is printed.

       export DOGSITE_BIN=/path/to/dogsitescorer-2.0.0/dogsite
       python3 "$REPO"/bin/polyphar-dogsite-worker --benchmark     # --tool dogsite2 (default)

   For the DoGSite3 sensitivity analysis instead, install it with the ZBH installer
   (`/opt/ZBH/dogsite3/dogsite3`, with `VERSION_INFO` next to it), activate its license once
   (`dogsite3 --license <key>`), then set `DOGSITE_BIN` to it and add `--tool dogsite3` here
   and `POLYPHAR_POCKET_TOOL=dogsite3` (or `--pocket-tool dogsite3`) for the runs.

3. **Build the image.** The tag you choose is what `POLYPHAR_VERSION` selects later. Also tag
   it `dev`: `tests/test_docker.py` looks for `polyphar:dev` and skips silently without it.

       docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) \
           -f docker/Dockerfile -t polyphar:1.0 -t polyphar:dev .
       export POLYPHAR_VERSION=1.0

4. **Check the machine.**
   - The PocketMatch binaries run in the image. `Exec format error` means the kernel lacks
     IA32 emulation (see [Server checks](#server-checks)).
   - The wrapper runs as your user, with the real pocket-detection binary.
   - RCSB and AlphaFold downloads work.

       pytest -W error tests/test_docker.py
       pytest -W error -m "dogsite and not regression" tests/test_stage_pockets.py
       pytest -W error tests/test_stage_structures.py -k end_to_end

5. **DoGSiteScorer 2 only: compare it with the 2023 pockets.** This tells you whether a newly
   licensed DoGSiteScorer 2 binary behaves like the one used in 2023 (DoGSite3 is a different
   method, so it skips; set `POLYPHAR_POCKET_TOOL=dogsite2` and `DOGSITE_BIN` to the
   DoGSiteScorer 2 binary). It reruns three net2 PLK1 chains and
   requires the same pocket names, with drugScore within 1e-4. It needs only a small part of
   the 2023 data, so copy that from the workstation (adjust the host name):

       mkdir -p ~/adep_csb/disease_subtypes/Liver_cancer/net2/{pockets,pdbfixer/P53350}
       rsync -av workstation:/mnt/c/experiments/adep_csb/disease_subtypes/Liver_cancer/net2/pockets/DrugScore.csv \
           ~/adep_csb/disease_subtypes/Liver_cancer/net2/pockets/
       rsync -av 'workstation:/mnt/c/experiments/adep_csb/disease_subtypes/Liver_cancer/net2/pdbfixer/P53350/{2ogq_A,4o56_A,4h5x_B}.pdb' \
           ~/adep_csb/disease_subtypes/Liver_cancer/net2/pdbfixer/P53350/
       POLYPHAR_POCKET_TOOL=dogsite2 ADEP_DATA=~/adep_csb \
           pytest -W error -m "dogsite and regression" tests/test_stage_pockets.py

   If this fails, the DoGSiteScorer version differs from 2023. Don't start a rerun until
   that is understood: pockets, and everything after them, would not be comparable.

6. **Set up Pushover and send a test notification** (see
   [Pushover notifications](#pushover-notifications) for the file format):

       mkdir -p ~/.polyphar && $EDITOR ~/.polyphar/pushover.env && chmod 600 ~/.polyphar/pushover.env
       set -a && . ~/.polyphar/pushover.env && set +a && pytest -m pushover tests/test_notify.py

7. **Calibrate the container stages.** Run this from any directory. `pockets` was already
   calibrated on the host in step 2.

       mkdir -p ~/runs && cd ~/runs && NPROCS=64 "$REPO"/bin/polyphar-docker benchmark

8. **Optional: run end to end with the real binary** (~minutes, two proteins, 3 ligands).
   Run it without `POLYPHAR_E2E_PMAX` first. If `select` stops with "no protein pair has
   Pmax >= 0.7", that is a property of these two test proteins, not a fault: rerun it with
   `POLYPHAR_E2E_PMAX=0.3`.

       cd "$REPO" && pytest tests/test_run_e2e.py

Then follow the quick start: create a run folder, check `run --dry-run`, and start the run
in `tmux`.

To run the seven v2 subtypes on POSEIDON one stage at a time, with the exact commands, the
expected ID counts and what to check after each stage, see **[docs/RUNBOOK.md](docs/RUNBOOK.md)**.

## Quick start

Build the image once from the repo root. The git commit is recorded in every stage manifest.

    docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) \
        -f docker/Dockerfile -t polyphar:dev .

Make a run folder (one per subtype), put the inputs in `inputs/`, and run everything:

    REPO=~/repos/bioinf_polyphar_workflow           # this repository
    mkdir -p ~/runs/net1/inputs
    cp clusters_features_ppi.xlsx e-Drug3D_2083.sdf ~/runs/net1/inputs/
    cd ~/runs/net1
    export POLYPHAR_VERSION=dev DOGSITE_BIN=/opt/ZBH/dogsite3/dogsite3 NPROCS=32
    "$REPO"/bin/polyphar-docker run --dry-run        # ETA table, does nothing
    "$REPO"/bin/polyphar-docker run --xlsx inputs/clusters_features_ppi.xlsx --sheet net1 \
        --id-column "uniprot Entry" --reviewed-column "Reviewed/unreviewed"

The HCC workbook `clusters_features_ppi.xlsx` has tabs `net1`, `net2`, `net3`. The Alzheimer
workbook has no reviewed column and needs `--gene-column`:

    "$REPO"/bin/polyphar-docker run --xlsx inputs/PPIN_ROSMAP_clusterFeatures_v7.xlsx \
        --sheet clusterN --id-column UNIPROTID --gene-column GeneID

One workbook row is one protein. Some AD cells list several accessions (the CALM3 row holds
CALM1, CALM2 and CALM3), and `--gene-column` makes the `uniprot` stage keep only the reviewed
accession whose primary gene name equals the row's gene. The stage stops with the candidates
listed if none or several match. HCC needs no gene column: it has one accession per row, and
every reviewed row is kept.

The wrapper uses the current directory as the run folder (or `$POLYPHAR_RUN_DIR`). It mounts the run folder at `/work` and `~/.polyphar` (or
`$POLYPHAR_HOME`) at `/polyphar_home`, runs as your user, and passes `NPROCS` through.

**Long runs.** A full subtype takes hours to days. Run it in `tmux` (or `nohup`) and keep the
wrapper in the foreground. A tmux session does not inherit the variables you exported in
your shell when the tmux server is already running, so pass them inside the command string
(or with `tmux new -e VAR=value`):

    cd ~/runs/net1
    tmux new -s net1 "POLYPHAR_VERSION=dev DOGSITE_BIN=/opt/ZBH/dogsite3/dogsite3 NPROCS=32 \"$REPO\"/bin/polyphar-docker run"
    # or
    mkdir -p logs && nohup "$REPO"/bin/polyphar-docker run > logs/run.out 2>&1 &

To stop a run, stop the wrapper (Ctrl-C in its tmux window, or `kill <wrapper pid>`). A raw
`docker stop` gives the container only 10 s before SIGKILL, which can cut the FAILED
notification and the run summary short; if you must, use `docker stop -t 60 <container>`.
An ssh drop (SIGHUP) is handled like SIGTERM: the run ends with its FAILED notification,
summary and pool cleanup (exit 129).

(Once `00_ids/ids.txt` exists, `run` needs no ids options: see below.)

Do not use `docker run -d`. The wrapper must stay alive because it supervises the host
pocket worker, and stops it when the container exits.

### `polyphar run`

`run [--from STAGE] [--to STAGE]` runs the stages in order, or the contiguous range you give.
It accepts every stage's own options (`--xlsx/--txt/--sheet/--id-column/--reviewed-column/--gene-column`,
`--ph`, `--worker-timeout`, `--worker-dead-after`, `--pm-bin-dir`, `--pmax`, `--exclude`,
`--hub-tie-break`, `--library`, `--seed`, `--chunk-size`, `--threads-per-job`,
`--chunk-timeout`), plus `--max-failed-fraction` and `--keep-going`.

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
    <run>/07_docking/<label>/poses.sdf  df_results_<label>_fixed.csv  obsolete/
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
    polyphar similarity [--pm-bin-dir /opt/pocketmatch]
    polyphar select [--pmax 0.7] [--exclude O15350] [--hub-tie-break union]
    polyphar dock --library inputs/e-Drug3D_2083.sdf [--seed 42] [--chunk-size 50] \
                  [--threads-per-job 1] [--chunk-timeout 43200]
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

## Old scripts -> subcommands (for Methods)

| 2023 script | subcommand |
|---|---|
| `get_uniprot_data.py` | `polyphar uniprot` (with `polyphar ids` building the list) |
| `get_pdb_af.py` | `polyphar structures` |
| `fix_pdb.py` | `polyphar fix` |
| `get_pockets.py` | `polyphar pockets` (+ host `bin/polyphar-dogsite-worker`) |
| `get_pocketsimilarity.py` (README: `get_PMsimilarities.py`) | `polyphar similarity` |
| `Gene_Disease_Association/target_selection_analysis.py` | `polyphar select` |
| per-folder smina scripts (`smina_screen_*.py`) | `polyphar dock` |

### Methods-relevant changes vs 2023

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
- **Target selection**: Pmax ≥ 0.7, top-degree hub + neighbours, O15350 excluded, as in 2023.
  Two changes: protein pairs are de-duplicated as **unordered** pairs, and a tie at the top
  degree takes the **union of the tied hubs** (`--hub-tie-break union`, the default), with one
  row per protein. In 2023 the hub was tied in 5 of 6 subtypes, so row order decided it.

## Docking parallelism

smina runs with exhaustiveness 12. With `--threads-per-job 1` (the default), `NPROCS` smina
processes run at once and each does its 12 searches serially. That is the same total work as
`--cpu 12`, and the results are byte-identical (`test_cpu_and_chunk_invariance`). For a tail
of peptide-sized ligands that dominate the run time, use `--threads-per-job 12` so the last
chunks use all cores. `NPROCS` also sets how many DoGSite3 / DoGSiteScorer jobs the host worker
runs at once. Results do not depend on `--chunk-size`, but changing it re-docks everything (the chunk
is part of the job key).

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

## Outputs and failures

- `failed.tsv` in each stage directory lists `item, stage, reason`. The manifest's `counts`
  has `failed`, the number of failed work units.
  - `structures`: unmapped (accession, PDB) pairs are listed in `failed.tsv` for visibility,
    but they are counted under `counts.unmapped`, not `failed`.
  - `pockets`: pockets dropped as outside the target's SIFTS segment are listed in
    `failed.tsv` and counted under `counts.outside_segment`, not `failed`.
  - `dock`: `failed` counts failed **chunks**. `incomplete_targets` counts targets with no
    results (missing ligands, receptor preparation failure, or a missing receptor). An
    incomplete target gets no `df_results_<label>_fixed.csv`, and `failed.tsv` lists its
    missing ligand IDs. Treat `incomplete_targets > 0` as a failure.
- When `select` drops a target on a rerun, its docking outputs are moved to
  `07_docking/obsolete/<timestamp>/`, not deleted.
- `df_results_<label>_fixed.csv` has `Pose, SMINA_Score, Ligand, ID, Pocket`. **Downstream
  code must key ligands on `ID`, not `Ligand`**: 35 different e-Drug3D drugs share the title
  `Molecule-1`.

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
