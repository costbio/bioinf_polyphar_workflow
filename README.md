# Bioinformatics-based polypharmacology workflow

This repository contains Python scripts to start a workflow for computational prediction of binding pocket similarities using <a href="https://github.com/1337deepesh/PocketMatch_v2.0">PocketMatch 2.0</a>
between protein structures of a list of genes provided as input.

The scripts should be run in the order given below. 

The operations performed by each of them are described below:

1. **get_uniprot_data.py** reads gene names from a file, performs UniProt queries to retrieve gene data, saves the data to a CSV file.
2. **get_pdb_af.py** reads the CSV file produced by **get_uniprot_data.py** , downloads protein structures of the genes from the PDB and AlphaFold databases, extracts relevant information from the structures, and performs various data processing tasks. The final output includes CSV files containing the merged protein information, chain information, and filtered mapping information between PDB chains and UniProt IDs. SIFTS database mappings found in pdb_chain_uniprot.tsv are used to match PDB chains with UniProt IDs.
3. **fix_pdb.py:** TO-DO
4. **get_pockets.py:** TO-DO
