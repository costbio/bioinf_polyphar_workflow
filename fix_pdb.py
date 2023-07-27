import os
import re
import argparse
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

def pdbfixer(out_path, in_path, ph=7.0):
   
    # Recursively traverse the directory tree and find all PDB files with the .ent extension

               
    pdb_list = []
    for root, dirs, files in os.walk(in_path):
        for file in files:
            if re.match(r'^(?<!AF)\w+\.pdb$', file):
                pdb_list.append(os.path.join(root, file))
    
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        
    for pdb_path in pdb_list:
        try:
            # Create a PDBFixer instance for the current file
            fixer = PDBFixer(filename=pdb_path)
            # Do something with the fixer instance, such as adding missing atoms
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(ph)
            fixer.removeHeterogens(keepWater=False)
            #fixer.removeChains(indices)
            # Get the protein directory relative to the protein_path
            protein_dir = os.path.dirname(os.path.relpath(pdb_path, in_path))
            # Construct the output file path
            output_path = os.path.join(out_path, protein_dir,  os.path.basename(pdb_path))
            # Make sure the directory exists
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            # Save the fixed PDB file
            with open(output_path, 'w') as outfile:
                PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
        except Exception as e:
            print(f"An error occurred while processing file {pdb_path}: {e}")
            # delete the error file
            os.remove(pdb_path)
            # save the name of the deleted file in a text file
            with open('fixer_deleted_files.txt', 'a') as f:
                f.write(os.path.basename(pdb_path) + '\n')
            continue

if __name__ == '__main__':
    # Define the command line arguments
    parser = argparse.ArgumentParser(description='Fix missing atoms in PDB files.')
    parser.add_argument('-out_path', type=str, help='output directory path')
    parser.add_argument('-in_path', type=str, help='input directory path')
    parser.add_argument('-p', '--ph', type=float, default=7.0, help='pH value for adding missing hydrogens')
    args = parser.parse_args()

    # Call the pdbfixer function with the command line arguments
    pdbfixer(args.out_path, args.in_path, args.ph)