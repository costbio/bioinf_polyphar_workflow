# Polypharmacology Workflow Modernization Implementation Plan (v2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the pocket-similarity → target-selection → eDrug3D-docking workflow into one tested command-line program, `polyphar`, that runs in a Docker image against a local run folder. It must:
- select chains by SIFTS mapping
- record the version of every data source
- run every stage in parallel
- print and notify an ETA before each stage starts
- send progress through Pushover

Then rerun HCC subtype 1 (net1) on its real gene list, and after it every other subtype, on the user's 64-CPU server.

**Architecture:**
- A Python package `polyphar` with one module per stage under `polyphar/stages/` and one argparse entry point, `polyphar <stage> [options]`.
- Every stage reads and writes only inside one **run folder**, mounted at `/work`, using a fixed numbered layout. Every stage first computes its exact work units, prints and notifies an ETA, and then runs with a progress reporter that sends milestone notifications.
- The Docker image holds the Python env, smina, and the upstream 32-bit PocketMatch binaries.
- DoGSiteScorer stays on the host. The `pockets` stage in the container writes jobs into the run folder, and a standard-library host worker runs them and writes back results and markers. The container keeps ownership of the ETA, progress, notifications and parsing.
- A small host wrapper, `bin/polyphar-docker`, makes the image behave like a local command, and starts and stops the host worker automatically.

**Tech Stack:** Python 3.12, pandas, networkx, requests, biopython, prody, pdbfixer/OpenMM 8, openbabel, smina 2020.12.10, pytest. DoGSiteScorer 2 (host-installed, licensed). PocketMatch v2.0 at commit `4648bff3c4d233104c6d1f4ced351b5f7944972a`. Docker on an amd64 base with i386 runtime libraries. Pushover HTTP API.

**Spec:** There is no separate spec. The requirements are the "Decisions" and "Verified facts" sections below, plus `/home/onur/repos/Gene_Disease_Association/NOTES_HCC_subtype1_rerun.md`.

## Decisions (user-confirmed 2026-10-31)

- Scope: fix the pipeline, run HCC net1 first, then rerun net2, net3 and AD Cluster1/3/4 with the same pipeline.
- Data versions: use the current UniProt, PDB, SIFTS and AlphaFold DB, and record every version in a per-stage `manifest.json`.
- Apply SIFTS chain mapping. Input lists are UniProt accessions, reviewed (Swiss-Prot) human entries only. Maybridge is out of scope.
- eDrug3D docking and target selection become stages of this workflow, and are included in the image.
- Execution is subagent-driven (user, 2026-10-31).
- Delivery is a Docker image that acts as a CLI on a local folder. The user controls the server: Docker is available, and DoGSiteScorer will be installed and licensed **on the server host** before any run.
- Every stage reports an ETA **before** it starts. Progress goes to Pushover. akusei's `pushover-bash` is a curl wrapper, not a Python package, so the program calls the Pushover HTTP API directly through `requests` and needs no new dependency.

## Verified facts this plan relies on (checked 2026-10-31)

**Structures and pockets**
- **The old chain rule was wrong.** Wrong-chain rates against SIFTS: net1 195/1314, net2 257/1145, net3 195/1314 chains. Final targets that sat on a different protein include:
  - HCC: 2MJW_A (RAGE, should be S100P chain B), 7KZR_A (FANCA, should be UBE2T chain X), 5CVD_A (NTMT1, should be CENPA chains D/E), 2ZIX_A (MUS81, should be EME1 chain B).
  - AD: 2HQW_A (rat CALM1), 4DIN_A (mouse PKA Cα), 6UM2_A (bovine IGF2R), 6KDR_A (NTMT2, should be CENPB chains D/E), 8HA0 (no chain A; PTH1R is chain R).
  - In 2HQW, SIFTS maps the GRIN1 chain to *rat* P35439, so no chain maps to human Q05586.
- **net1 input:** the `net1` tab of `clusters_features_ppi.xlsx` (columns `gene name, uniprot Entry, Reviewed/unreviewed, …`) has 461 reviewed rows, which are **459 unique accessions**. The AD workbook `PPIN_ROSMAP_clusterFeatures_v7.xlsx` has tabs `cluster1..4`, a `UNIPROTID` column, and no reviewed column.
- **UniProt REST:** `https://rest.uniprot.org/uniprotkb/stream?query=(accession:X) OR …&fields=accession,id,gene_primary,reviewed,organism_id,xref_pdb&format=tsv`. It returns the headers `Entry, Entry Name, Gene Names (primary), Reviewed, Organism (ID), PDB`, with `;`-separated PDB IDs and a trailing `;`. The response header `x-uniprot-release` gives the release (`2026_03`).
- **AlphaFold DB:** `https://alphafold.ebi.ac.uk/api/prediction/{acc}` returns `latestVersion` 6 and a `pdbUrl`. Plain `curl` of the file URLs returned **HTTP 403** from the workstation, which is unresolved (Task 5). The 2023 runs used AF **v2**, and 4–5 "AF" files per 2023 run are saved HTTP error pages.
- **SIFTS:** `https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz` returns 200. Its first line is a comment, so read it with `header=1`. Columns: `PDB CHAIN SP_PRIMARY RES_BEG RES_END PDB_BEG PDB_END SP_BEG SP_END`.
- **PocketMatch:**
  - The upstream prebuilt **32-bit** binaries (`alpha-file_maker/Step0-cabbage-core`, `alpha-file_maker/Step0-END-FILE`, `Step3-PM_serial`) reproduce the upstream `sample_outputs` exactly in `i386/debian:bookworm-slim`: the cabbage file is byte-identical and 120/120 Pmin/Pmax match.
  - The native 64-bit build aborts with `*** buffer overflow detected ***`.
  - Names are truncated to 32 characters; the longest 2023 names are exactly 32, with no collisions.
  - The 2023 net2 scores have 5 numeric columns while upstream writes 6. Only columns 0–3 are used.
- **Case bug:** pocket files on disk keep the chain case (`P53350_2ogq_A.pdb_res_P_1.pdb`), but the legacy `DrugScore.csv` lowercases it. This only resolves on case-insensitive NTFS.
- **Pocketsim contract:** the columns are `,UniProt ID,Gene Name,Protein,Pocket,Pocket_Path,UniProt ID_2,Gene Name_2,Protein_2,Pocket_2,Pocket_Path_2,Pmin,Pmax`. `Protein` is `4O56_A` for PDB entries and `AF-P06493-F1` for AlphaFold. `Pocket` is `P_1` or `P_12_1`. `Gene Name` is the UniProt Entry Name (`PLK1_HUMAN`).

**Target selection**
- The canonical code is `Gene_Disease_Association/target_selection_analysis.py`, which produced `target_selection_results/*_final_targets.csv`. The older notebooks in `disease_subtypes/Network analysis/` are superseded and disagree with it: `analysis .ipynb` filters on Pmin, and `subtype_network.ipynb` reports the betweenness hub. For each subtype, the script:
  1. Sorts by Pmax descending and runs `drop_duplicates(['UniProt ID','UniProt ID_2'])`. This works on the **ordered** pair, so (A,B) and (B,A) can both survive, depending on file order.
  2. Keeps rows with `Pmax >= 0.7`.
  3. Builds an undirected networkx graph whose nodes are pockets, named `f"{UniProt}_{EntryName}_{Protein}_{Pocket}"`.
  4. Ranks nodes by `Degree`, with an unspecified tie order.
  5. Takes the top node as the hub, skipping `O15350` (p73, which is "problematic for docking"); the targets are the hub plus all its neighbours, with neighbours that are O15350 also skipped.
  6. Writes `final_targets.csv` with the columns `Node, UniProt_ID, Gene_Name, Organism, Structure_ID, Similarity`, where the hub gets `Similarity` 1.0. Chain and pocket survive **only inside `Node`**.
- `parse_node_string` finds the structure as the first `_`-part, after the UniProt ID, that starts with `AF-` or is 4 characters starting with a digit. That misreads Entry Names such as `2AAA_HUMAN` (Task 11 tests this).
- **The top hub was tied in 5 of 6 subtypes in 2023**, so the hub was decided by row order:
  - net1: PI15 (`O43692_-_AF-O43692-F1_P_1`) and TTK, both degree 12 → PI15
  - net2: p73 and S100P, degree 8 → S100P, after the p73 skip
  - net3: p73, S100P and TTK, degree 8 → S100P
  - AD Cluster3: TGFB1 and CENPB, degree 8 → TGFB1
  - AD Cluster4: every node has degree 1 (a 2-node network)

  The legacy AlphaFold-only proteins have the Entry Name `-` (e.g. `O43692_-_…`), because the old `get_pdb_af.py` took gene names only from PDB rows.
- **Decided (user, 2026-10-31) — hub ties: take the union of all tied hubs.** Every non-excluded node with the top degree is a hub, and the targets are all those hubs plus all their neighbours, de-duplicated. Implemented as `--hub-tie-break {union,weighted,betweenness}`, **default `union`**. Every tie is always reported in the manifest and the notification.

**Docking**
- Legacy eDrug3D docking lives in `docking_screen_data/eDrug3D_DOCK/eDrug3D_DOCK/<Disease>/<target>/`, with 34 folders of hand-edited, non-identical script copies.
- Settings:
  - receptor converted with pybel (`CorrectForPH(7.4)`, `addh`, partial charges)
  - box = the pocket's atom bounding box, unpadded
  - `smina --num_modes 12 --exhaustiveness 12 --cpu 12`, with an unseeded random seed that is logged and reused per ligand
  - library IDs assigned **by position**
- The receptor is the structure the pocket was computed on: the pdbfixer output for PDB chains, the raw model for AlphaFold.
- **Defects:**
  - **2xnn**: only 93 of 2,083 ligands were docked.
  - **TTK (AF-P33981-F1)** was docked against a different library version (223 unknown names), so its IDs are wrong.
  - 35 distinct molecules are titled `Molecule-1`, and `docking_screen_analysis.ipynb` groups by name.
  - 2xnn and TTK are targets of net1 and net3.
- **e-Drug3D library:** md5 prefix `65b16127`, 2,083 records, with fields `status, cas, name, ID`. `ID` is unique; `name` is not.
- **smina:** conda-forge `smina 2020.12.10`. The 2023 notebook read `df_results_<target>_fixed.csv` (`Pose, SMINA_Score, Ligand, ID, Pocket`) and never used the RMSD columns. **That file is retired** (see the output contract below): the engine is selectable, so an engine-specific column name cannot survive, and one score column cannot express a distribution over replicates. A leftover `_fixed` file is now an orphan and the dock stage sweeps it into `07_docking/obsolete/<ts>/`.

**Environment**
- Workstation: Docker 29.8 works and IA32 emulation works in WSL2.
- Server: user-controlled, 64 CPUs, Docker present.

## Global Constraints

- Python `3.12`. `environment.yml` is the only dependency source, for both the image and any host env. The package is installed with `pip install -e .`, which provides the console script `polyphar = polyphar.cli:main`.
- **Run-folder contract.** One run folder = one subtype. Every stage reads and writes only inside the run folder, which is `/work` in the image; `--run-dir` defaults to `/work`, else the current directory. Inputs a run needs are copied into `<run>/inputs/` first. The layout:
  ```
  <run>/inputs/            ids source (xlsx or txt), e-Drug3D library sdf
  <run>/00_ids/ids.txt
  <run>/01_uniprot/uniprot.csv  not_found.txt  renamed.tsv
  <run>/02_structures/<UniProt>/{pdb}.cif {pdb}_{chain}.pdb AF-<acc>-F1.pdb ; Protein_info.csv ; sifts/
  <run>/03_pdbfixer/<UniProt>/{pdb}_{chain}.pdb
  <run>/04_pockets/<UniProt>/… ; DrugScore.csv
  <run>/05_similarity/Pocketsim.csv  PocketMatch_score.txt
  <run>/06_targets/final_targets.csv  centrality.csv  unique_pairs.csv
  <run>/07_docking/<label>/poses.sdf ; <label>/replicates/r<i>.sdf ; df_results_<label>.csv
  <run>/logs/<stage>.log
  ```
  Every stage directory gets `manifest.json` and `failed.tsv` (`item, stage, reason`).
- **Never delete or modify inputs.** Never write outside the run folder, except `$POLYPHAR_HOME` (calibration) and temp dirs. Every stage is **resumable**: finished work units are skipped on rerun.
- **ETA before start.** Each stage computes its exact work units from its inputs, then prints and notifies `units, estimated duration, expected finish time` before doing any work. `--dry-run` prints that and exits. `polyphar estimate` prints the next stage exactly and the later stages as projections, labelled "projected".
- **Parallelism.** `--nprocs` defaults to `$NPROCS`, else 4. Results must be identical for any `--nprocs`.
- **Notifications** are optional and must never fail a run. Credentials come only from the env vars `PUSHOVER_APP_TOKEN` and `PUSHOVER_USER_KEY`, never from the image or the repo.
- **Kept legacy parameters:**
  - DoGSiteScorer flags `-s -i -y -d -w 4`; `drugScore >= 0.5`
  - AlphaFold `F1` models only, and **not** passed through pdbfixer
  - pdbfixer at pH 7.0, adding missing residues, atoms and hydrogens, and removing heterogens including water
  - PocketMatch filters: drop `NULL____`, drop same-UniProt pairs. **The 2023 `Pmin or
    Pmax == 1` drop was NOT kept** (commit `08d6e67`): it discarded cross-protein
    containment hits above the selection threshold. A chain-stem guard replaces the
    self-comparison coverage it incidentally provided.
  - target selection with `Pmax >= 0.7`, top-degree hub + neighbours, exclude `O15350`
  - docking with `--exhaustiveness 12 --num_modes 12`, receptor pH 7.4, unpadded pocket box, ligand SDF as supplied
- **New fixed choices:**
  - Strict SIFTS accession matching (`SP_PRIMARY == accession`, no ortholog fallback); one chain per (UniProt, PDB), chosen by maximum SP coverage, ties to the smallest chain ID.
  - Target selection de-duplicates on the **unordered** protein pair. Hub ties use `--hub-tie-break`, default `union` (user decision): every top-degree non-excluded node is a hub.
  - Docking uses a fixed recorded `--seed 42`, and ligand identity comes from the library `ID`.
- Output contracts kept for downstream notebooks:
  - `Pocketsim.csv` columns as above
  - `final_targets.csv` keeps the legacy columns and **adds** `Protein, Chain, Pocket, Pocket_Path, Receptor_Path`
  - `df_results_<label>.csv` (**not** `_fixed`) has the columns `Ligand, ID, Pocket, n, score_mean, score_sd, score_se, score_min, score_max`: one row per ligand, summarising that ligand's BEST score in each of the `n` replicates. `score_sd` is the sample sd (ddof=1), `score_se = sd/sqrt(n)`; both are empty at `n == 1`. The per-replicate pose records stay in `<label>/replicates/r<i>.sdf`, so any other summary can be recomputed without re-docking.
  - **Downstream code must key ligands on `ID`, not `Ligand`**: 35 different e-Drug3D drugs share the title `Molecule-1`.
- **DoGSiteScorer runs only on the host** (user decision 2026-10-31: the binary's license is host-bound). It is never copied into or mounted into the image, and never committed to the repo.
  - The container hands work to the host through a **file job queue** in the shared run folder, `<run>/04_pockets/queue/`.
  - `bin/polyphar-dogsite-worker` is a host-side script that needs only Python ≥ 3.8 with the standard library. It consumes the queue and runs the binary given by `--dogsite` (default `$DOGSITE_BIN`).
  - There are no sockets, no privileged containers, and no Docker-socket mounts.

## Review Focus

1. **A PDB entry listed in UniProt where no chain maps to that accession** (ortholog, e.g. 2HQW for Q05586). Skip and log it; never fall back to the first chain. Test: Task 5.
2. **A case-sensitive filesystem**, as in Docker and on the server. Pocket files must be found by their exact on-disk name. Test: Task 8.
3. **A docking job that dies part-way** (2xnn), or a ligand with no poses. Mark the target incomplete, list the missing IDs, and write no partial `df_results_*.csv`. Test: Task 12.
4. **Pocket names over 32 characters, or names that collide after truncation.** Fail before PocketMatch runs. Test: Task 10.
5. **An interrupted run that is restarted** (ssh drop, container killed). Finished units are skipped, half-written outputs are redone, and the new ETA counts only the remaining units. Tests: Tasks 5, 6, 8, 10 and 12 (the `resume` tests).

---

## File structure

| File | Responsibility |
|---|---|
| `pyproject.toml`, `environment.yml` (new) | Package metadata and the console script; pinned deps |
| `polyphar/cli.py` | argparse subcommands `ids uniprot structures fix pockets similarity select dock run estimate benchmark`, plus the common flags `--run-dir --nprocs --dry-run --no-notify` |
| `polyphar/common.py` | `RunLayout` (all paths above), `write_manifest`, `sha256_file`, `http_session`, `append_failure`, `PM_NAME_LIMIT` |
| `polyphar/eta.py` | `Calibration` load/save, `Estimate`, `estimate_stage`, `Progress` reporter |
| `polyphar/notify.py` | `notify(title, message, priority=0)` through the Pushover API |
| `polyphar/stages/{ids,uniprot,structures,fix,pockets,similarity,select,dock}.py` | One stage each; each exposes `units(layout, args) -> WorkUnits` and `run(layout, args, progress) -> dict` |
| `polyphar/benchmark.py`, `polyphar/data/bench/` | Fixture workloads that calibrate per-unit costs on the current machine |
| `docker/Dockerfile`, `bin/polyphar-docker` | Image; host wrapper (also starts and stops the worker) |
| `bin/polyphar-dogsite-worker` | Host-only, standard-library DoGSiteScorer queue worker |
| `tests/` | pytest; markers `network`, `regression` (`ADEP_DATA`), `dogsite` (`DOGSITE_BIN`, host only), `pocketmatch` (`PM_BIN_DIR`), `docking` (smina on PATH), `docker`, `pushover` (credentials set) |
| Removed | `get_uniprot_data.py`, `get_pdb_af.py`, `fix_pdb.py`, `get_pockets.py`, `get_pocketsimilarity.py`, `pdb_chain_uniprot.tsv`. The README maps each old script to its subcommand, for the manuscript's Methods. |

---

### Task 1: Package skeleton, CLI, run layout, manifests

**Files:**
- Create: `pyproject.toml`, `environment.yml`, `pytest.ini`, `polyphar/__init__.py`, `polyphar/cli.py`, `polyphar/common.py`, `tests/conftest.py`, `tests/test_common.py`

**Interfaces:**
- Produces:
  - `RunLayout(run_dir: Path)`, with properties for every path in the layout (`inputs`, `ids_txt`, `uniprot_csv`, `structures`, `protein_info`, `pdbfixer`, `pockets`, `drugscore_csv`, `similarity`, `pocketsim_csv`, `targets`, `final_targets_csv`, `docking`, `logs`) and `stage_dir(stage: str) -> Path`
  - `WorkUnits(count: int, cost_basis: str, sizes: list[float] | None)`. `sizes` gives per-unit weights, e.g. file bytes, used for ETA and the longest-job bound.
  - `sha256_file(path) -> str`
  - `http_session(retries=5, backoff=2.0) -> requests.Session`. Retries on 429/5xx, with the User-Agent `polyphar (+https://github.com/costbio/bioinf_polyphar_workflow)`.
  - `append_failure(stage_dir, item, stage, reason) -> None`
  - `write_manifest(stage_dir, stage, args: dict, inputs: list[Path], versions: dict, counts: dict, timing: dict) -> Path`
  - `PM_NAME_LIMIT = 32`
  - `cli.main(argv=None) -> int`, which registers the subcommands with their stage modules. The common flags are `--run-dir` (default `/work` if it exists, else cwd), `--nprocs` (default `int($NPROCS or 4)`), `--dry-run` and `--no-notify`.

- [ ] **Step 1:** Write `environment.yml` (env `polyphar`, from conda-forge): `python=3.12 pandas numpy networkx requests biopython prody pdbfixer openmm openbabel smina=2020.12.10 openpyxl pytest`. Then run `conda env create -f environment.yml && conda run -n polyphar pip install -e .`, and pin the versions that resolve.
- [ ] **Step 2:** Write failing tests:
  - `test_layout_paths`: `RunLayout(tmp).pocketsim_csv == tmp/"05_similarity/Pocketsim.csv"`.
  - `test_manifest_records_inputs_versions_timing`: the JSON has `stage`, the input sha256, `versions`, `counts`, `timing.estimated_s`, `timing.actual_s` and a non-empty `git_commit`.
  - `test_append_failure_header_once`.
  - `test_cli_lists_subcommands`: `polyphar --help` exits 0 and names all 11 subcommands.
- [ ] **Step 3:** Run `conda run -n polyphar pytest tests/test_common.py -v`. Expected: FAIL.
- [ ] **Step 4:** Implement. Subcommands without a stage yet print `not implemented` and exit 2.
- [ ] **Step 5:** Run `pytest tests/test_common.py -v`. Expected: PASS.
- [ ] **Step 6:** Commit: `build: polyphar package skeleton, run layout, manifests`

### Task 2: ETA engine, progress reporter, Pushover

**Files:**
- Create: `polyphar/eta.py`, `polyphar/notify.py`, `tests/test_eta.py`, `tests/test_notify.py`

**Interfaces:**
- Consumes: `WorkUnits` (Task 1).
- Produces:
  - `Calibration.load(home: Path) -> Calibration`, `.cost(stage) -> float | None` (seconds per unit-basis on one process), `.update(stage, seconds_per_unit)`, `.save()`.
    - The file is `$POLYPHAR_HOME/calibration.json`, default `~/.polyphar`, mounted at `/polyphar_home`.
    - Each entry stores `seconds_per_unit`, `n_observations`, the machine `nproc` and CPU model, and `updated`.
    - `.update` uses a running mean.
  - `estimate_stage(stage: str, units: WorkUnits, nprocs: int, calib: Calibration, parallel_efficiency: dict) -> Estimate(seconds: float | None, basis: str)`.
    - `seconds = max(sum(sizes) * c / slots, max(sizes) * c)` when `sizes` is set, else `count * c / slots`.
    - `slots` is set per stage: `nprocs` for most stages; `nprocs / 2` for similarity (the block-pair overhead); `nprocs // threads_per_job` for docking; `min(nprocs, 8)` for network stages.
    - When there is no calibration it returns `seconds=None` with `basis="uncalibrated — run polyphar benchmark"`.
  - `format_eta(est, now) -> str`, e.g. `"~3 h 20 min (finish ≈ Sat 14:05)"`.
  - `Progress(stage, total, est, notifier, milestones=(0.05, 0.25, 0.5, 0.75))`, with `.advance(n=1, size=None)`, `.done(counts)` and `.fail(exc)`.
    - It logs a line at most every 60 s.
    - At each milestone it sends a notification with the live-rate ETA. The 5% one is "revised ETA".
    - `.done` notifies the elapsed time against the estimate, the counts and the number of failures, and updates the calibration from the measured rate.
    - `.fail` notifies with priority 1 and the last 10 log lines.
  - `notify(title: str, message: str, priority: int = 0) -> bool`.
    - POSTs to `https://api.pushover.net/1/messages.json` with `token`, `user`, `title` (truncated to 250) and `message` (truncated to 1024).
    - Returns False, logging one warning, when the credentials are unset or the request fails. It never raises.
    - The title is prefixed with `[<run folder name>]`.

- [ ] **Step 1:** Write failing tests:
  - `test_estimate_uses_largest_job_bound`: sizes `[100, 1, 1]`, c=1, slots=8 → `100`.
  - `test_estimate_uncalibrated_none`.
  - `test_progress_milestones_fire_once`: with a stub notifier, advancing 100 units sends exactly 4 milestone calls plus none extra.
  - `test_progress_done_updates_calibration`.
  - `test_notify_without_credentials_returns_false`, with no network call made.
  - `test_notify_truncates`: a stubbed session sees `len(message) <= 1024`.
  - `@pushover test_notify_live`: sends `"polyphar test notification"` and expects `status == 1`.
- [ ] **Step 2:** Run `pytest tests/test_eta.py tests/test_notify.py -v`. Expected: FAIL.
- [ ] **Step 3:** Implement. Check the Pushover limits (1024/250) against the current API docs, and correct them here if they differ.
- [ ] **Step 4:** Run the tests. Expected: PASS. With credentials exported, the live test delivers one message to the phone.
- [ ] **Step 5:** Commit: `feat: ETA engine with calibration, milestone progress, pushover notifications`

### Task 3: `polyphar ids` — reproducible input lists

**Files:**
- Create: `polyphar/stages/ids.py`, `tests/test_stage_ids.py`

**Interfaces:**
- Produces: `build_id_list(xlsx, sheet, id_column, reviewed_column | None) -> list[str]`. It splits on `;`, `,` and whitespace, strips, de-duplicates in first-seen order, and keeps rows whose reviewed column equals `reviewed` (case-insensitive) when that column is given.
- CLI: `polyphar ids --xlsx inputs/<file> --sheet <tab> --id-column <col> [--reviewed-column <col>]` or `polyphar ids --txt inputs/<file>`. Writes `00_ids/ids.txt` and a manifest recording the sheet and columns. The units are trivial (no ETA line).

- [ ] **Step 1:** Write failing tests:
  - `test_dedup_reviewed_and_order`: rows `[("P1","reviewed"),("P2","unreviewed"),(" P1 ","Reviewed"),("P3;P4","reviewed")]` → `["P1","P3","P4"]`.
  - `@regression test_hcc_tabs_match_known_counts`: net1 → 459 IDs, with a set different from `Liver_cancer/net1/net1UniprotID.txt`; net2 equals `net2/net2.txt` (156); net3 equals `net3/net3.txt` (197).
- [ ] **Step 2:** Run the tests. Expected: FAIL. **Step 3:** Implement. **Step 4:** Run with `ADEP_DATA=/mnt/c/experiments/adep_csb`. Expected: PASS.
- [ ] **Step 5:** Commit: `feat(ids): reviewed accession lists from cluster workbooks with provenance`

### Task 4: `polyphar uniprot`

**Files:**
- Create: `polyphar/stages/uniprot.py`, `tests/test_stage_uniprot.py`, `tests/fixtures/uniprot_sample.tsv` (the real response for `P25815 Q9NPD8`, plus one unreviewed and one non-human row)

**Interfaces:**
- Produces:
  - `parse_uniprot_tsv(text) -> DataFrame`, with the trailing `;` of `PDB` stripped
  - `fetch_uniprot(accessions, session, batch=100) -> (DataFrame, release: str)`
  - `classify(requested, df) -> dict(kept, not_found, unreviewed, non_human, renamed)`. `renamed` maps a requested ID to its primary accession; detect these by re-querying the missing IDs one at a time.
  - `units() = WorkUnits(count=ceil(n/100), cost_basis="batch")`
  - Outputs in `01_uniprot/`: `uniprot.csv` (kept rows, same column names), `not_found.txt`, `renamed.tsv`, and `manifest.json` (versions: uniprot release).

- [ ] **Step 1:** Write failing tests: `test_parse_strips_trailing_semicolon` (P25815 PDB == `"1J55;1OZO;2MJW;7NMI"`), `test_classify_splits_categories`, `test_classify_reports_renamed` (a stubbed re-query maps `OLD1` → `P25815`), and `@network test_fetch_live` (the release matches `^\d{4}_\d{2}$` and the PDB list contains `2MJW`).
- [ ] **Step 2:** Run. Expected: FAIL. **Step 3:** Implement. **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: `feat(uniprot): UniProt REST stage with unreviewed/renamed reports`

### Task 5: `polyphar structures` — SIFTS chains, downloads, extraction

**Files:**
- Create: `polyphar/stages/structures.py`, `tests/test_stage_structures.py`
- Fixtures:
  - `tests/fixtures/sifts_subset.tsv`: the real rows for 2MJW, 7KZR, 5CVD, 2HQW and 1Q4K, plus a synthetic homodimer with equal A/B coverage, keeping the comment first line
  - `tests/fixtures/2mjw.cif`
  - `tests/fixtures/af_api_P25815.json`
- Delete: `get_uniprot_data.py`, `get_pdb_af.py`, `pdb_chain_uniprot.tsv`

**Interfaces:**
- Consumes: `01_uniprot/uniprot.csv` (`Entry`, `Entry Name`, `PDB`).
- Produces:
  - `download_sifts(dest_dir, session) -> (Path, last_modified: str)`
  - `load_sifts(path) -> DataFrame`
  - `select_chains(pairs[UniProt ID, PDB], sifts) -> (selected[UniProt ID, PDB, Chain ID, SP coverage], unmapped[UniProt ID, PDB, reason])`. The reason is `no SIFTS chain for accession` or `PDB absent from SIFTS`.
  - `download_mmcif(pdb_id, dest_dir, session) -> Path`. Gets `https://files.rcsb.org/download/{ID}.cif` and checks that the body starts with `data_`.
  - `extract_chain(cif, chain, out_pdb) -> dict(Resolution, Structure method)`. ProDy `parseMMCIF(chain=…)` → `writePDB`, keeping all models. A multi-character chain is written as `A`, while the filename keeps the original ID.
  - `fetch_alphafold(acc, dest_dir, session) -> (Path | None, version | None)`. Uses the API and the entry with `entryId == f"AF-{acc}-F1"`. On 404 or no F1 entry it writes nothing. It rejects bodies that don't start with a PDB record.
  - `units() = WorkUnits(count = n_pdb_entries + n_accessions, cost_basis="download")`
  - Outputs: `02_structures/` with `Protein_info.csv` (`UniProt ID, Gene Name, PDB, Chain ID, AlphaFold, Resolution, Structure method, SIFTS coverage, File`), `sifts/`, `failed.tsv` and `manifest.json` (SIFTS Last-Modified, the AF versions seen). Downloads run through a thread pool with `min(nprocs, 8)` workers, and existing valid files are skipped.

- [ ] **Step 1: Resolve the AlphaFold 403 first.** Run `conda run -n polyphar python -c "from polyphar.common import http_session; r=http_session().get('https://alphafold.ebi.ac.uk/files/AF-P25815-F1-model_v6.pdb'); print(r.status_code, r.text[:6])"`, on the workstation **and** on the server. Expected: `200 HEADER`. **If either gets 403, STOP and report**: do not change the User-Agent or work around a bot check. The candidate to bring to the user is the AFDB per-proteome archive on `ftp.ebi.ac.uk/pub/databases/alphafold/`, which is not yet verified.
- [ ] **Step 2:** Write failing tests:
  - `test_selects_target_chain_not_first`: (P25815, 2MJW) → B, (Q9NPD8, 7KZR) → X, (P49450, 5CVD) → D.
  - `test_ortholog_entry_unmapped` (Review Focus 1): (Q05586, 2HQW) is unmapped.
  - `test_absent_pdb_unmapped`.
  - `test_tie_breaks_alphabetically`.
  - `test_extract_chain_writes_only_requested_chain`: on `2mjw.cif` chain B.
  - `test_af_missing_model_writes_nothing`, `test_af_error_body_rejected`, `test_af_picks_F1`.
  - `test_resume_skips_existing` (Review Focus 5): a second run makes 0 downloads, and `units()` reports only the missing ones.
  - `@regression test_extraction_matches_2023_where_chain_was_right`: 1q4k chain B equals `Liver_cancer/net2/structures/P53350/1q4k_B.pdb`.
  - `@network test_end_to_end_two_proteins`: P25815 and Q9NPD8 → rows for 2MJW/B and 7KZR/X, and no row for 2MJW/A.
- [ ] **Step 3:** Run. Expected: FAIL. **Step 4:** Implement; delete the old scripts and the TSV. **Step 5:** Run. Expected: PASS.
- [ ] **Step 6:** Commit: `fix(structures): SIFTS-mapped chain extraction, AFDB API, resumable non-destructive downloads`

### Task 6: `polyphar fix`

**Files:**
- Create: `polyphar/stages/fix.py`, `tests/test_stage_fix.py`, `tests/fixtures/2mjw_B.pdb` (from Task 5)
- Delete: `fix_pdb.py`

**Interfaces:**
- Produces:
  - `select_inputs(structures_dir) -> list[Path]`: `*.pdb` files not starting with `AF-`
  - `fix_one(src, structures_dir, out_dir, ph=7.0) -> None`, using `from openmm.app import PDBFile` and the legacy call sequence
  - `units() = WorkUnits(count, cost_basis="bytes", sizes=[file sizes of the not-yet-fixed inputs])`
  - A process pool with `nprocs` workers. Output is written to a temp name and renamed when complete, so a killed run never leaves a truncated file that later looks done. Failures go to `failed.tsv`, and the input is kept.

- [ ] **Step 1:** Write failing tests: `test_select_inputs_excludes_af`, `test_fix_adds_hydrogens_and_drops_hets`, `test_failure_keeps_input`, and `test_resume_ignores_partial_tmp` (a leftover `*.tmp` is redone, and a completed file is skipped).
- [ ] **Step 2:** Run. Expected: FAIL. **Step 3:** Implement. **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: `fix(pdbfixer): parallel, atomic writes, never delete inputs`

### Task 7: Docker image and host wrapper

**Files:**
- Create: `docker/Dockerfile`, `bin/polyphar-docker`, `tests/test_docker.py`

**Interfaces:**
- Produces: image `polyphar:<version>`, with:
  - the env from `environment.yml` and the package installed
  - `ENTRYPOINT ["polyphar"]`
  - `/opt/pocketmatch/` holding the **upstream prebuilt 32-bit** binaries from the pinned commit, plus `Sample_pockets/` and `sample_outputs/`
  - `ENV PM_BIN_DIR=/opt/pocketmatch POLYPHAR_HOME=/polyphar_home`
  - `WORKDIR /work`
  - i386 runtime libraries (`libc6:i386 libstdc++6:i386`)
  - no DoGSiteScorer, and no credentials
- Base: a pinned `mambaorg/micromamba:<version>-debian12-slim` tag; confirm it exists with `docker pull` first.
- `bin/polyphar-docker [polyphar args…]` runs `docker run --rm --user "$(id -u):$(id -g)"` with:
  - `-v "$(realpath "${POLYPHAR_RUN_DIR:-$PWD}")":/work`
  - `-v "${POLYPHAR_HOME:-$HOME/.polyphar}":/polyphar_home`
  - `--env-file "${POLYPHAR_HOME:-$HOME/.polyphar}/pushover.env"`, only if that file exists
  - `-e NPROCS`
  - `${POLYPHAR_DOCKER_OPTS}`
  - and finally `polyphar:${POLYPHAR_VERSION:-latest} "$@"`

  Outputs are owned by the calling user. The container always runs in the foreground; for long runs, run the wrapper under `tmux` or `nohup`.

  **Worker lifecycle:** when the subcommand is `pockets`, or `run`/`benchmark` with a stage range that includes `pockets`, the wrapper:
  1. checks that `$DOGSITE_BIN` is executable, and exits with a clear message if not
  2. starts `bin/polyphar-dogsite-worker --run-dir <host run dir> --dogsite "$DOGSITE_BIN" --nprocs "${NPROCS:-4}" --watch` in the background
  3. stops the worker through a `trap` on EXIT/INT/TERM

  The worker's log goes to `<run>/logs/dogsite-worker.log`.

- [ ] **Step 1:** Write failing `@docker` tests:
  - `test_pocketmatch_reproduces_upstream_sample`: inside the image, the 16 sample pockets give 120/120 Pmin/Pmax identical to `sample_outputs`.
  - `test_wrapper_writes_into_local_folder`: `bin/polyphar-docker ids --txt inputs/ids.txt` in a tmp run folder creates `00_ids/ids.txt` owned by the current UID.
  - `test_image_has_no_dogsite_or_credentials`: no file named `dogsite*` exists in the image, and `env` has no `PUSHOVER_*`.
- [ ] **Step 2:** Run `pytest tests/test_docker.py -v`. Expected: FAIL. **Step 3:** Write the Dockerfile and wrapper; build with `docker build -f docker/Dockerfile -t polyphar:dev .`. **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: `build: docker image as CLI with upstream 32-bit PocketMatch, host wrapper`

### Task 8: `polyphar pockets` + host DoGSiteScorer worker (file queue)

**Files:**
- Create: `polyphar/stages/pockets.py`, `bin/polyphar-dogsite-worker`, `tests/test_stage_pockets.py`, `tests/test_dogsite_worker.py`
- Fixtures:
  - `tests/fixtures/fake_dogsite`: an executable stub that, given `-p <model> … -o <out>`, writes `<out>_desc.txt` (copied from a real net2 desc file) and one `<out>_res_P_1.pdb`. It exits 3 when the model's name contains `FAIL`, and sleeps when `FAKE_DOGSITE_SLEEP` is set.
  - `tests/fixtures/pockets_net2_P53350/`: real `*_desc.txt` files plus empty `*_res_*.pdb` placeholders with their exact on-disk names, and one AlphaFold example.
- Delete: `get_pockets.py`

**Queue protocol** (all paths in job files are **relative to the run folder**, so the container's `/work` and the host path both resolve them):
- **Directory:** `04_pockets/queue/`.
- **Job file:** `<job_id>.job` holds the JSON `{"model": "03_pdbfixer/P25815/2mjw_B.pdb", "out": "04_pockets/P25815/P25815_2mjw_B.pdb", "size": <bytes>}`. The container writes it through a temp name plus `os.replace`.
- **Claiming:** the worker claims a job by `os.rename(<id>.job → <id>.running.<hostname>.<pid>)`, which is atomic on one filesystem, so two workers can never take the same job.
- **Success:** if dogsite exits 0 and `<out>_desc.txt` exists, the worker writes `<id>.done` with the JSON `{"exit": 0, "seconds": s, "dogsite_sha256": …}`.
- **Failure:** otherwise it writes `<id>.failed` with `{"exit": n, "stderr_tail": …}`. Input files are never touched.
- **Stale claims:** at start, a `.running.<this host>.<pid>` whose pid is no longer alive is renamed back to `.job`.
- **Heartbeat:** the worker touches `queue/WORKER.<hostname>.<pid>` every 30 s.
- **Shutdown:** `--watch` keeps polling every 5 s until `queue/CLOSED` exists. The container writes that file once every job is done or failed.

**Interfaces:**
- Consumes: `03_pdbfixer/` (PDB chains) and `02_structures/` (the `AF-*.pdb` files and `Protein_info.csv`).
- Produces, container side (`polyphar/stages/pockets.py`):
  - Output naming: `04_pockets/<UniProt>/<UniProt>_<pdb>_<chain>.pdb` for PDB chains, `04_pockets/<UniProt>/AF-<acc>-F1.pdb` for AlphaFold.
  - `units() = WorkUnits(count, "bytes", sizes=[the model sizes without a _desc.txt])`.
  - `enqueue(layout) -> list[str]`: writes the job files for the not-yet-done models; this is idempotent.
  - `wait_for_queue(layout, progress, worker_timeout=600) -> dict(done, failed)`.
    - It reports progress from the `.done`/`.failed` files, weighted by `size`.
    - If no heartbeat is newer than `worker_timeout`, it sends **one** priority-1 notification, "No DoGSiteScorer worker is consuming the queue — start bin/polyphar-dogsite-worker on the host", and keeps waiting. Interrupting it is safe, because the stage is resumable.
    - Failed jobs go to `failed.tsv`.
  - `parse_pockets(pockets_dir, protein_info) -> DataFrame`: the legacy `DrugScore.csv` columns, with `Pocket Filename` built from the **actual output stem**, and every referenced file checked to exist.
  - The manifest records `dogsite_sha256`, taken from the `.done` files; if they disagree, the stage fails.
- Produces, host side: `bin/polyphar-dogsite-worker --run-dir <host path> --dogsite <binary> [--nprocs 4] [--watch] [--benchmark]`.
  - It uses only the Python ≥ 3.8 standard library (`subprocess`, `concurrent.futures`, `json`, `os`).
  - It runs `[dogsite, "-p", model, "-s", "-i", "-y", "-d", "-w", "4", "-o", out]` with `cwd` set to the output dir, up to `--nprocs` at a time.
  - **Preflight:** before claiming anything, it runs the binary on `polyphar/data/bench/2mjw_B.pdb` in a temp dir. If that fails, it exits non-zero with the binary's stderr and the hint `license/installation problem — see README §DoGSiteScorer`, without claiming any jobs.
  - `--benchmark` times the preflight on the three bench chains and writes `pockets.seconds_per_unit` (per byte) into `$POLYPHAR_HOME/calibration.json`. It writes through temp + replace, and keeps the other entries.

- [ ] **Step 1:** Write failing tests:
  - `test_pocket_filename_preserves_case` (Review Focus 2, in `tmp_path` on ext4).
  - `test_af_rows_named_without_prefix`.
  - `test_enqueue_idempotent_and_relative_paths`.
  - `test_resume_counts_only_missing_models`.
  - `test_worker_runs_queue_with_stub` (with `fake_dogsite`): every job is `.done` and the outputs are in place.
  - `test_worker_failure_recorded`: a model named `FAIL_A.pdb` gives a `.failed` with exit 3, and the input still exists.
  - `test_two_workers_never_double_claim`: 2 workers and 40 jobs give exactly 40 `.done` files, and the stub's invocation log shows 40 distinct models.
  - `test_stale_running_reclaimed`: a `.running.<host>.<dead pid>` becomes `.job` and is then done.
  - `test_preflight_failure_claims_nothing`: a stub that always fails leaves every job as `.job`, with a non-zero worker exit.
  - `test_wait_warns_once_without_worker`: a stub notifier with `worker_timeout=1` gets exactly one warning.
  - `test_wrapper_starts_and_stops_worker` (`@docker`): with `DOGSITE_BIN` set to the stub, `bin/polyphar-docker pockets` on a fixture run folder finishes, and no worker process remains afterwards.
  - `@dogsite test_real_dogsite_on_fixture` (needs `DOGSITE_BIN`; host only).
  - `@dogsite @regression test_matches_2023_pockets`: 3 files from `Liver_cancer/net2/pdbfixer/` give the same pocket names, with `drugScore` within 1e-4 of `net2/pockets/DrugScore.csv`.
- [ ] **Step 2:** Run `pytest tests/test_stage_pockets.py tests/test_dogsite_worker.py -v`. Expected: FAIL; the `@dogsite` tests SKIP.
- [ ] **Step 3:** Implement both sides, and give the worker its own `--help` text.
- [ ] **Step 4:** Run. Expected: PASS for everything except the `@dogsite` tests, which run on the server in Task 14. If the 2023 regression differs, STOP and report (it may be a DoGSiteScorer version difference).
- [ ] **Step 5:** Commit: `feat(pockets): host-side DoGSiteScorer via run-folder job queue; case-exact parsing`

### Task 9: `polyphar benchmark` — per-machine calibration

**Files:**
- Create: `polyphar/benchmark.py`, `polyphar/data/bench/` (`2mjw_B.pdb` plus two more extracted chains of different sizes, 16 PocketMatch sample pockets, `edrug3d_head20.sdf`, a 2mjw receptor and pocket), `tests/test_benchmark.py`

**Interfaces:**
- Consumes: the `fix` / `pockets` / `similarity` / `dock` internals (Tasks 6, 8, 10, 12) through their single-unit functions. It also uses `Calibration` (Task 2).
- Produces: `polyphar benchmark [--stages …]`, which runs each fixture workload on one process and stores `seconds_per_unit` for:
  - `fix`: per byte. `pockets` is calibrated on the host by `bin/polyphar-dogsite-worker --benchmark` (Task 8); the wrapper calls it when `benchmark` includes `pockets`.
  - `similarity`: per pocket pair, measured as all-vs-all of the 16 sample pockets = 120 pairs
  - `dock`: per ligand at exhaustiveness 12 on one thread
  - `uniprot` and `structures`: per request, timing 5 live requests

  It also stores the projection priors used by `polyphar estimate`: chains per accession, model bytes, kept pockets per model and targets per subtype. These come from the 2023 runs and are hard-coded with their source noted, then replaced by the observed values after each real stage. Stages that are unavailable, such as DoGSiteScorer without a mount, are reported as skipped. The whole benchmark takes under 5 minutes on 1 core.
- `polyphar estimate` prints a table for the stages after the last completed one. The next stage is exact (from `units()`); the later ones are projected from the priors.

- [ ] **Step 1:** Write failing tests: `test_benchmark_writes_calibration` (the `fix` and `similarity` entries exist and are > 0; inside the container `pockets` is reported as "host-calibrated" and left unchanged), and `test_estimate_table_marks_projection` (on a run folder with only `01_uniprot` done, the `structures` row is "exact" and `dock` is "projected").
- [ ] **Step 2:** Run. Expected: FAIL. **Step 3:** Implement. Do this task after Tasks 10 and 12 if they aren't done yet; the order in this list is for reading. **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: `feat: per-machine benchmark calibration and whole-run ETA projection`

### Task 10: `polyphar similarity` — PocketMatch, serial and block-pair parallel

**Files:**
- Create: `polyphar/stages/similarity.py`, `tests/test_stage_similarity.py`, `tests/fixtures/pm_scores_sample.txt` (8 real net2 lines: an AF–PDB pair, a same-UniProt pair, a `NULL____` row, a cross-protein Pmin==1 row and a same-UniProt Pmin==Pmax==1 row) with matching small `DrugScore.csv` / `Protein_info.csv`
- Delete: `get_pocketsimilarity.py`

**Interfaces:**
- Produces:
  - `select_pockets(pockets_dir, cutoff=0.5) -> list[Path]`
  - `check_pocket_names(pockets) -> None`. Raises, listing the names, when one is longer than 32 characters or two collide.
  - `run_pocketmatch(pockets, pm_bin_dir, work_dir) -> Path`. Runs `Step0-cabbage-core <abs path>` per pocket into `outfile.cabbage`, then `Step0-END-FILE`, then `Step3-PM_serial outfile.cabbage`, checking every exit code.
  - `plan_blocks(n, nprocs) -> list[(start, end)]`: the smallest `k` with `k(k-1)/2 >= nprocs`, and `k=1` for `nprocs=1`.
  - `run_pocketmatch_parallel(pockets, pm_bin_dir, work_dir, nprocs) -> Path`. Each block pair `i<j` runs all-vs-all on `B_i ++ B_j` and keeps the cross pairs. Within-block pairs come from exactly one job. The merged lines are ordered as the serial run would emit them.
  - `parse_scores(score_file, pockets, protein_info) -> DataFrame`: the Pocketsim contract. Filters NULL, same-UniProt and same-chain-stem pairs; the 2023 `== 1` drop is not kept (`08d6e67`).
  - `units() = WorkUnits(count = n*(n-1)//2, cost_basis="pair")`. A completed block-pair job leaves a `.done` marker and is skipped on rerun.
  - Outputs in `05_similarity/`: `Pocketsim.csv`, `PocketMatch_score.txt` and `manifest.json` (PocketMatch commit).

- [ ] **Step 1:** Write failing tests:
  - `test_name_guard_rejects_long_and_colliding` (Review Focus 4).
  - `test_parse_scores_schema_and_filters`.
  - `test_plan_blocks_cover_all`.
  - `@pocketmatch test_parallel_equals_serial_on_sample`: `nprocs=4` is byte-identical to `nprocs=1`.
  - `@pocketmatch @regression test_reproduces_net2_2023`: on `Liver_cancer/net2/pockets`, every pair's Pmin/Pmax matches `net2/PM_net2/PocketMatch_score.txt` (names compared case-insensitively, trailing `_` stripped). The rows keyed on `(UniProt ID, Protein, Pocket, UniProt ID_2, Protein_2, Pocket_2, Pmin, Pmax)` must be a **superset** of `Pocketsim_net2.csv` (every 2023 row still present; the extra rows all cross-protein with Pmin == 1 -- 5,678 of them on net2), for both `nprocs=1` and `nprocs=16`. Equality no longer holds since `08d6e67`. Report both wall-clock times.
- [ ] **Step 2:** Run. Expected: FAIL. **Step 3:** Implement. **Step 4:** Run `PM_BIN_DIR=… ADEP_DATA=… pytest tests/test_stage_similarity.py -v`, inside the image. Expected: PASS.
- [ ] **Step 5:** Commit: `fix(similarity): upstream PocketMatch in temp dirs, block-pair parallel, name guard`

### Task 11: `polyphar select` — target selection

**Files:**
- Create: `polyphar/stages/select.py`, `tests/test_stage_select.py`, `tests/fixtures/pocketsim_small.csv` (about 12 hand-written rows: a hub with 3 neighbours, an O15350 node with the top degree, a reversed duplicate (B,A) of an (A,B) pair, a degree tie, and an Entry Name `2AAA_HUMAN`)

**Interfaces:**
- Consumes: `05_similarity/Pocketsim.csv`.
- Produces:
  - `unique_protein_pairs(df, pmax=0.7) -> DataFrame`. Keys on the **unordered** pair `frozenset(UniProt ID, UniProt ID_2)` and keeps the max-Pmax row, with ties going to the smallest `(Node, Node_2)`. Then filters `Pmax >= pmax`.
  - `build_graph(df) -> nx.Graph`. Node IDs follow the legacy format `f"{UniProt}_{EntryName}_{Protein}_{Pocket}"`. Each node carries the attributes `uniprot, entry_name, protein, pocket, pocket_path`; **nothing is parsed back from the string**.
  - `rank_hubs(G, tie_break: str) -> DataFrame` with `Node, UniProt_ID, Gene_Name, Degree, Weighted_Degree, Degree_Centrality, Betweenness_Centrality, Closeness_Centrality`.
    - With `tie_break="weighted"` it sorts by `Degree desc, Weighted_Degree desc, Node asc`.
    - With `tie_break="betweenness"` it sorts by `Degree desc, Betweenness desc, Node asc`.
    - With `tie_break="union"` it ranks as `weighted`, and `select_targets` then uses every non-excluded top-degree node as a hub, with targets = the hubs plus all their neighbours, de-duplicated.
    - Also `tied_top(ranking) -> list[str]`.
  - `select_targets(G, ranking, exclude: set[str], tie_break: str, hub: str | None = None) -> DataFrame`. `hub` overrides the hub; it is used only by the regression test. The legacy columns `Node, UniProt_ID, Gene_Name, Organism, Structure_ID, Similarity`, plus `Protein, Chain, Pocket, Pocket_Path, Receptor_Path, Is_Hub`.
    - `Gene_Name`/`Organism` split the Entry Name at its last `_`.
    - `Structure_ID` is the PDB ID, or the `AF-…-F1` ID.
    - `Receptor_Path` is `03_pdbfixer/<UniProt>/<pdb>_<chain>.pdb` for PDB entries and `02_structures/<UniProt>/AF-<acc>-F1.pdb` for AlphaFold.
    - The hub is the first non-excluded node in the ranking and gets `Similarity` 1.0; excluded neighbours are dropped.
  - CLI: `polyphar select [--pmax 0.7] [--exclude O15350] [--hub-tie-break union]`. The default `--exclude` is `O15350`, recorded in the manifest with the reason "p73: problematic for docking (legacy decision)". Writes `06_targets/{final_targets,centrality,unique_pairs}.csv`. The manifest and the notification name any tie at the top. The units are trivial.

- [ ] **Step 1:** Write failing tests:
  - `test_reversed_pair_deduplicated`: (A,B) and (B,A) give one edge.
  - `test_hub_skips_excluded_and_neighbours_exclude_it`.
  - `test_union_of_tied_hubs`: with two non-excluded nodes tied at the top degree, `tie_break="union"` gives targets = both hubs (`Is_Hub` True, `Similarity` 1.0) + the union of their neighbours, with no duplicates; a neighbour of both hubs gets its **maximum** Pmax to any hub as `Similarity`. With `tie_break="weighted"` there is exactly one hub.
  - `test_tie_break_deterministic`: shuffling the input rows 10 times gives the same hub every time.
  - `test_entry_name_starting_with_digits`: for the `2AAA_HUMAN` node, `Gene_Name == "2AAA"` and `Structure_ID` is the real PDB ID.
  - `test_final_targets_columns_and_receptor_paths`.
  - `@regression test_reproduces_2023_targets`: run on the 2023 `Pocketsim_net2.csv`, `Pocketsim_net3.csv` and `Alzhemier/Cluster{1,3,4}/PM_*/Pocketsim.csv`, then compare the `Node` set and `Similarity` with `Gene_Disease_Association/target_selection_results/<disease>_<subtype>_final_targets.csv`. Pass `hub=<the legacy hub>` through a test-only override, because the legacy hubs were picked by row order out of ties. Expected: equal. If a subtype differs, the test prints the nodes that differ and whether the unordered-pair change or a hub tie explains it. **STOP and report** before changing anything.
- [ ] **Step 2:** Run. Expected: FAIL. **Step 3:** Implement. **Step 4:** Run with `ADEP_DATA=…`. Expected: PASS, or a reported, explained difference.
- [ ] **Step 5:** Commit: `feat(select): target selection stage with order-independent pairs, deterministic hub, explicit receptor paths`

### Task 12: `polyphar dock` — eDrug3D smina stage

**Files:**
- Create: `polyphar/stages/dock.py`, `tests/test_stage_dock.py`
- Fixtures: `tests/fixtures/dock_2mjw/` holding the legacy `2mjw_A.pdb`, `P25815_2mjw_A.pdb_res_P_1.pdb` and the first 20 pose-1 rows of `df_results_2mjw.csv`, plus `edrug3d_head20.sdf`. These are used only to check that the legacy run is reproduced.

**Interfaces:**
- Consumes: `06_targets/final_targets.csv` (`Node, Pocket_Path, Receptor_Path, Protein, Pocket`) and `inputs/<library>.sdf`.
- Produces:
  - `read_library(sdf) -> DataFrame[ID, name, cas, status, record_index]`. Text-based parsing; raises on a missing or duplicate `ID`.
  - `write_titled_library(sdf_in, sdf_out)`: records copied verbatim, with each title set to `EDRUG3D_<ID>`.
  - `receptor_to_pdbqt(pdb, pdbqt, ph=7.4)`: the legacy pybel steps.
  - `pocket_box(pocket_pdb) -> (center, size)`: unpadded.
  - `run_smina(receptor, ligands, center, size, out_sdf, log, seed, cpu, exhaustiveness=12, num_modes=12)`. Raises on a non-zero exit.
  - `parse_poses(out_sdf, library, label) -> DataFrame[Pose, Score, Ligand, ID, Pocket]` (`Score`, not `SMINA_Score`: the engine is selectable).
  - Jobs are (target × chunk of `--chunk-size 50` ligands), run with `--threads-per-job 1` (default) under a `ProcessPoolExecutor(nprocs // threads_per_job)`. A chunk is done when its `.done` marker exists.
  - Results are merged per target into `07_docking/<label>/poses.sdf` and `07_docking/df_results_<label>.csv`, where `label = f"{Protein}_{Pocket}"`.
  - `units() = WorkUnits(count = n_remaining_chunks × chunk_size, cost_basis="ligand")`.
  - The manifest records `smina --version`, the library md5 and the seed.
  - CLI: `polyphar dock --library inputs/e-Drug3D_2083.sdf [--seed 42] [--chunk-size 50] [--threads-per-job 1]`.

- [ ] **Step 1: Pin legacy equivalence and determinism first.** Write `@docking` tests:
  - `test_reproduces_legacy_2mjw_scores`: dock `edrug3d_head20.sdf` into the fixture receptor and box with `seed=-496630900` (from the legacy log), `cpu=12`, exhaustiveness 12 and 12 modes. The pose-1 score per ligand must equal the legacy value for all 20.
  - `test_cpu_and_chunk_invariance`: one job at `cpu=12` gives the same table as 2 chunks at `cpu=1`.
- [ ] **Step 2:** Run. Expected: FAIL. **Step 3:** Implement the single-job functions. **Step 4:** Run.
  - Expected: PASS.
  - If the legacy test fails, **STOP and report** the score differences; don't tune parameters to force agreement.
  - If the invariance test fails, set the default `--threads-per-job 12` and record why in the README.
- [ ] **Step 5:** Write failing tests:
  - `test_library_ids_unique_and_complete`.
  - `test_placeholder_titles_resolved_by_id`: two records titled `Molecule-1` give two distinct IDs with their real names.
  - `test_incomplete_target_flagged` (Review Focus 3): a chunk missing 3 IDs → `failed.tsv` lists them, and no `df_results` file is written for that target.
  - `test_resume_skips_done_chunks` (Review Focus 5): `units()` after a partial run counts only the remaining chunks.
- [ ] **Step 6:** Implement the rest. Run `pytest tests/test_stage_dock.py -v`. Expected: PASS.
- [ ] **Step 7:** Commit: `feat(dock): parametrised smina stage, ID-keyed ligands, chunked parallel and resumable`

### Task 13: `polyphar run`, README

**Files:**
- Modify: `polyphar/cli.py`, `README.md`, `/mnt/c/experiments/adep_csb/CLAUDE.md` (the stage 1 section)
- Test: `tests/test_run_e2e.py`

**Interfaces:**
- Produces: `polyphar run [--from <stage>] [--to <stage>]`, which runs `ids → uniprot → structures → fix → pockets → similarity → select → dock` in order.
  - Before each stage it prints and notifies that stage's ETA, which is exact at that point.
  - When the next stage is uncalibrated, it runs `benchmark --stages <that stage>` automatically first.
  - It stops at the first failed stage, with a priority-1 notification.
  - At the end it notifies a summary table: per stage, the estimate, the actual time and the failures.
  - `polyphar run --dry-run` equals `polyphar estimate`.
- README sections:
  - quick start with `bin/polyphar-docker`
  - the run-folder layout
  - one command per stage
  - old script → subcommand mapping (for Methods)
  - running long jobs on the server: `tmux new -s net1 'NPROCS=32 bin/polyphar-docker run'`, or `nohup … > logs/run.out &`. The wrapper must stay alive because it supervises the host DoGSiteScorer worker.
  - Pushover setup: a `~/.polyphar/pushover.env` file with mode 600, holding `PUSHOVER_APP_TOKEN` and `PUSHOVER_USER_KEY`
  - DoGSiteScorer: host install, `DOGSITE_BIN`, the host worker and queue protocol, how to run the worker by hand (`bin/polyphar-dogsite-worker --run-dir … --watch`), and the license note

- [ ] **Step 1:** Write a failing `@network @dogsite @pocketmatch @docking test_e2e_two_proteins`. A run folder with `inputs/ids.txt` = P25815, P53350 and a 20-ligand library goes through `polyphar run` and produces all 8 stage dirs with manifests, a non-empty `final_targets.csv` and docking results. A stub notifier receives one start/ETA message per stage before that stage's first unit, verified by timestamp order.
- [ ] **Step 2:** Run. Expected: FAIL, or SKIP off-server. **Step 3:** Implement and write the docs. **Step 4:** Run on the server. Expected: PASS.
- [ ] **Step 5:** Commit: `feat: polyphar run orchestration with per-stage ETA gate; docs`

### Task 14: Server setup (operational — on the 64-CPU server)

- [ ] **Step 1:** Install and license DoGSiteScorer on the host, and `export DOGSITE_BIN=/path/to/dogsite` in `~/.polyphar/env.sh`. Check: `python3 bin/polyphar-dogsite-worker --run-dir /tmp/dsw --dogsite "$DOGSITE_BIN"` passes its preflight. Then run `pytest -m dogsite` on the host with `DOGSITE_BIN` set. Expected: PASS, including the 2023 regression. If that differs, STOP and report.
- [ ] **Step 2:** Clone the repo and run `docker build -f docker/Dockerfile -t polyphar:1.0 .`. Run `pytest -m "docker"` on the host. Expected: PASS. This also proves IA32 emulation on the server kernel. An `Exec format error` means the kernel lacks IA32 emulation: STOP, because PocketMatch would then have to run on the workstation.
- [ ] **Step 3:** Run `NPROCS=32 bin/polyphar-docker benchmark`. This calibrates the container stages, and the host worker's `--benchmark` calibrates `pockets`. Expected: `~/.polyphar/calibration.json` has entries for every stage.
- [ ] **Step 4:** Create `~/.polyphar/pushover.env` (mode 600) and run `pytest -m pushover` inside the image. Expected: a phone notification.
- [ ] **Step 5:** Run the full suite in the image with `NPROCS=32`. Record the similarity and docking speedups from the Task 10 and 12 tests in the README.

### Task 15: Rerun HCC subtype 1 (operational)

- [ ] **Step 1:** Create the run folder `~/polyphar_runs/Liver_cancer/net1/inputs/` with `clusters_features_ppi.xlsx` and `e-Drug3D_2083.sdf`. Then run `bin/polyphar-docker ids --xlsx inputs/clusters_features_ppi.xlsx --sheet net1 --id-column "uniprot Entry" --reviewed-column "Reviewed/unreviewed"`. Expected: 459 IDs.
- [ ] **Step 2:** Run `bin/polyphar-docker estimate`, and report the whole-run projection to the user **before** starting.
- [ ] **Step 3:** Run `tmux new -s net1 'NPROCS=32 bin/polyphar-docker run'`. Expected: Pushover messages per stage. Afterwards, report from the manifests:
  - counts and `failed.tsv` sizes
  - proteins with at least one pocket at drugScore ≥ 0.5 (the corrected "95")
  - the hub, any tie, and the number of targets
  - docking completeness
- [ ] **Step 4:** Copy the run to `/mnt/c/experiments/adep_csb/disease_subtypes_v2/Liver_cancer/net1/` with `rsync -a`. Leave the old folders untouched.
- [ ] **Step 5:** Hand off to `Gene_Disease_Association`, which is out of scope here. `docking_screen_analysis.ipynb` must read `07_docking/df_results_<label>.csv` (columns `Ligand, ID, Pocket, n, score_mean, score_sd, score_se, score_min, score_max` -- **not** the retired `_fixed` file, which the dock stage now sweeps into `obsolete/`), parse `label` (`<Protein>_<Pocket>`), and group by `ID`, not `Ligand`. `target_selection_analysis.py` is superseded by `polyphar select`.

### Task 16: Rerun net2, net3, AD Cluster1/3/4 (operational — after Task 15)

- [ ] **Step 1:** For each subtype, create its own run folder and run `ids`:
  - HCC: `--sheet net2|net3` with the same columns as Task 15.
  - AD: `--xlsx PPIN_ROSMAP_clusterFeatures_v7.xlsx --sheet clusterN --id-column UNIPROTID`; the `uniprot` stage drops unreviewed entries. Compare each AD list with the existing `Alzhemier/ClusterN/ClusterN.txt`, and **STOP and report** if it differs beyond unreviewed entries.
- [ ] **Step 2:** For each subtype: `estimate`, then report, then `run`, then copy to `disease_subtypes_v2/<Disease>/<subtype>/`.
- [ ] **Step 3:** Report a table per subtype of proteins, chains, pockets ≥ 0.5, pairs, hub and targets, next to the 2023 numbers.
