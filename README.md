# polyphar: bioinformatics-based polypharmacology workflow

`polyphar` takes a list of UniProt accessions (one disease subtype) and runs:

    ids -> uniprot -> structures -> fix -> pockets -> annotate -> select -> similarity -> ligands -> dock

1. reviewed human UniProt records,
2. PDB chains (SIFTS-mapped) and AlphaFold models,
3. pdbfixer repair,
4. DoGSite3 pockets (DoGSiteScorer 2 as an alternative),
5. [PocketMatch 2.0](https://github.com/1337deepesh/PocketMatch_v2.0) pocket-pair similarities,
6. per-pocket geometric annotation, then target selection by explicit pocket filters,
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

## Old scripts -> subcommands (for Methods)

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

### Methods-relevant changes vs 2023

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

The last stage tests the study's claim -- *approved HCC drugs bind more of the HCC panel,
at a useful level, than comparable molecules do* -- in a comparison where **ligand size and
pocket generosity cannot produce the result**. It is pre-registered in
`docs/superpowers/specs/2026-11-04-preregistration.md` and designed in
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
