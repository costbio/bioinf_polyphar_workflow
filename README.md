# Bioinformatics-based polypharmacology workflow

This repository contains Python scripts to start a workflow for computational prediction of binding pocket similarities using <a href="https://github.com/1337deepesh/PocketMatch_v2.0">PocketMatch 2.0</a>
between protein structures of a list of genes provided as input.

The scripts should be run in the order given below. 

The operations performed by each of them are described below:

1. **get_uniprot_data.py** reads gene names from a file, performs UniProt queries to retrieve gene data, saves the data to a CSV file.
2. **get_pdb_af.py** reads the CSV file produced by **get_uniprot_data.py** , downloads protein structures of the genes from the PDB and AlphaFold databases, extracts relevant information from the structures, and performs various data processing tasks. The final output includes CSV files containing the merged protein information, chain information, and filtered mapping information between PDB chains and UniProt IDs. SIFTS database mappings found in pdb_chain_uniprot.tsv are used to match PDB chains with UniProt IDs.
3. **fix_pdb.py** takes a directory path containing PDB files downloaded by **get_pdb_af.py**, fixes missing atoms in each PDB file using the pdbfixer library, and saves the fixed PDB files in the specified output directory.
4. **get_pockets.py** performs pocket detection on protein structures prepared by **fix_pdb.py**, utilizing DoGSiteScorer 2, and generates a consolidated CSV file with the detected pockets.
5. **get_PMsimilarities.py** receives the CSV file produced by **get_pockets.py**, computes similarity scores between each pair of pockets, and produces a CSV file, which is the final output of the workflow.

## Docker image and wrapper

Build the image from the repo root (the git commit is recorded in stage manifests):

    docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) -f docker/Dockerfile -t polyphar:dev .

Run any subcommand against a run folder (outputs are owned by you):

    POLYPHAR_RUN_DIR=/path/to/run POLYPHAR_VERSION=dev bin/polyphar-docker ids --txt inputs/ids.txt

The image contains no DoGSiteScorer and no credentials; Pushover settings are read from `$POLYPHAR_HOME/pushover.env` on the host. Tests: `pytest tests/test_docker.py` (needs the built `polyphar:dev` image).

## Docking (`polyphar dock`)

`polyphar dock --library inputs/e-Drug3D_2083.sdf [--seed 42] [--chunk-size 50] [--threads-per-job 1]` docks every ligand into every target in `06_targets/final_targets.csv` with smina. It keeps the legacy settings: `--exhaustiveness 12 --num_modes 12`, and an unpadded box over the pocket's atoms. Each ligand is identified by its library `ID`, never by its title or name.

- **Receptors are protonated once, by pdbfixer at pH 7.0.**
  - A PDB receptor is the fix-stage output.
  - An AlphaFold model gets hydrogens only, from pdbfixer `addMissingHydrogens(7.0)` with a seeded RNG.
  - The PDBQT is written with Open Babel. It keeps those hydrogens, uses Gasteiger charges, and is rigid. It is byte-deterministic, and its sha256 is recorded in the manifest.
  - The legacy `CorrectForPH(7.4)` + `addh` step is gone. It re-placed hydrogens at random on every conversion, and that alone moved scores by up to 0.9 kcal/mol.
- **`--threads-per-job` defaults to 1.** A ligand's result does not depend on smina's `--cpu`, on its chunk, or on its position in the chunk (`test_cpu_and_chunk_invariance`). Single-threaded jobs on `nprocs` processes are therefore safe and use the cores best.
- **Incomplete targets.** If any ligand ID of a target has no poses, the target is incomplete. `07_docking/failed.tsv` lists the missing IDs, and no `df_results_<label>_fixed.csv` is written for that target.
- **Ligand hydrogens.** smina's ligand hydrogen addition is pinned explicitly with `--addH 1`, smina's default. It is recorded in the manifest and in the job key.
- **`--chunk-size`.** Results do not depend on the chunk size, but changing it re-docks everything.
- **`--chunk-timeout`.** The default is 43200 s. A smina call that runs longer is killed, and its chunk is failed with `timeout after N s` plus the path to its log.
- **Dropped targets.** When `select` drops a target, its old `df_results_<label>_fixed.csv` and `<label>/` folder are moved, not deleted, to `07_docking/obsolete/<timestamp>/`. They are listed under `obsoleted` in the manifest.
