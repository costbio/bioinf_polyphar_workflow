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
