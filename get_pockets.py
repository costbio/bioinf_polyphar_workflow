import io
import pandas as pd
import numpy as np
import subprocess
import os
import glob
import re
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import argparse

def getPocketsSingleCore(argument_chunk):

    pdb_list_chunk = argument_chunk[0]
    ds_path = argument_chunk[1]
    out_path = argument_chunk[2]

    #print(pdb_list_chunk)
    for model_path in pdb_list_chunk:
        #print(model_path)
        protein_name = os.path.splitext(os.path.basename(model_path))[0]
        protein_folder = os.path.basename(os.path.dirname(model_path))
        output_folder = os.path.join(out_path, protein_folder)
        output_file = os.path.join(output_folder, "{}_{}.pdb".format(protein_folder, protein_name))
        os.makedirs(output_folder, exist_ok=True)
        print(ds_path+"/dogsite", "-p", model_path, "-s", "-i", "-y", "-d", "-w", "4", "-o", output_file)
        p = subprocess.Popen([ds_path+"/dogsite", "-p", model_path, "-s", "-i", "-y", "-d", "-w", "4", "-o", output_file])
        p.wait()

def getPockets(out_path, pdb_path, alphafold_path, ds_path, nprocs=4):

    pdb_list = []
    for root, dirs, files in os.walk(pdb_path):
        for file in files:
            if re.match(r'.+\.(pdb)$', file):
                pdb_list.append(os.path.join(root, file))

    # Add alphafolds from a different path
    pdb_list2 = []
    for root, dirs, files in os.walk(alphafold_path):
        for file in files:
            if re.match(r'^AF-.+\.pdb$', file):
                pdb_list2.append(os.path.join(root, file))
    pdb_list.extend(pdb_list2)
    
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # Split pdb_list into nprocs chunks.
    pdb_list_chunks = np.array_split(pdb_list, nprocs)

    argument_chunks = zip(pdb_list_chunks,[ds_path]*nprocs,[out_path]*nprocs)
    with ThreadPoolExecutor(max_workers=nprocs) as executor:
        for chunk in argument_chunks:
            future = executor.submit(getPocketsSingleCore, chunk)

    for dirpath, dirnames, filenames in os.walk(out_path):
        for filename in filenames:
            if "AF" in filename:
                old_file_name = os.path.join(dirpath, filename)
                new_file_name = os.path.join(dirpath, "AF" + filename.split("AF")[1])
                os.rename(old_file_name, new_file_name)
                    
def parsePockets(folder_path, alphafold_path):
    # Get a list of all subdirectories in the folder
    subdirs = [os.path.join(folder_path, d) for d in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, d))]

    # Define an empty DataFrame
    df = pd.DataFrame()

    # Iterate over the subdirectories and read each text file into a DataFrame
    for subdir in subdirs:
        txt_files = glob.glob(os.path.join(subdir, '*_desc.txt'))
        for txt in txt_files:
            # Read the text file into a DataFrame
            df2 = pd.read_csv(txt, sep='\t')

            # Get the basename of the file
            txt_basename = os.path.basename(txt)
            

            # Add columns to the DataFrame with the basename of the file and the folder name
            folder_name = os.path.basename(subdir)
            df2['Protein'] = txt_basename
            df2['UniProt ID'] = folder_name
        
            # Concatenate the current DataFrame with the previous ones
            df = pd.concat([df, df2], ignore_index=True)

    # Rename the 'frame' and 'folder' columns and save the DataFrame to a CSV file
    col_c = df.pop('Protein')
    df.insert(0, 'Protein', col_c)
    col_c = df.pop('UniProt ID')
    df.insert(0, 'UniProt ID', col_c)
    df = df.rename(columns={'name': 'Pocket'})
    df['Protein'] = df['Protein'].apply(lambda x: x.split('_', 1)[1] if '_' in x and not x.startswith('AF') else x)
    #df['Protein'] = df['Protein'].str[3:]
    df['Protein'] = df['Protein'].str.upper()
    a1=df[["Protein","1"]]= df['Protein'].str.split(".", n=1, expand=True)
    df=df.drop("1", axis=1)
    df['Protein'] = df['Protein'].str.replace('-F1', '')

    # Read df_pdb_af from pdb_path
    df_pdb_af = pd.read_csv(os.path.join(alphafold_path,'Protein_info.csv'))

    gene_to_id = dict(zip(df_pdb_af['UniProt ID'], df_pdb_af['Gene Name']))
    # Add a new column to the DataFrame that contains the corresponding ID for each gene name
    df['Gene Name'] = df['UniProt ID'].map(gene_to_id)
    col = df.pop('Gene Name')
    df.insert(1, 'Gene Name', col)
    df['Pocket Filename'] = df.apply(lambda x: x['Protein'] + ('-F1' if 'AF-' in x['Protein'] else '') + '.pdb_res_' + x['Pocket'] + '.pdb' if 'AF-' in x['Protein'] else x['UniProt ID'] + '_' + x['Protein'].lower() + '.pdb_res_' + x['Pocket'] + '.pdb', axis=1)
    df.to_csv(os.path.join(folder_path, 'DrugScore.csv'))

    return df

if __name__ == '__main__':
    # Define the command line arguments
    parser = argparse.ArgumentParser(description='Detect putative ligand binding pockets on input protein structures using Dogsitescorer.')
    parser.add_argument('-out_path', type=str, help='output directory path')
    parser.add_argument('-pdb_path', type=str, help='input directory path containing PDB structures.')
    parser.add_argument('-alphafold_path', type=str, help='input directory path containing AF2 structures.')
    parser.add_argument('-ds_path', type=str, help='input directory path containing  dogsitescorer executable.')
    parser.add_argument('-nprocs', type=int, default=4, help='number of processors to be used in calculation.')
    args = parser.parse_args()

    out_path = os.path.abspath(args.out_path)
    pdb_path = os.path.abspath(args.pdb_path)
    alphafold_path = os.path.abspath(args.alphafold_path)
    ds_path = os.path.abspath(args.ds_path)
    nprocs = int(args.nprocs)
    
    getPockets(out_path=out_path, pdb_path=pdb_path, alphafold_path=alphafold_path, ds_path=ds_path, nprocs=nprocs)

    parsePockets(out_path, alphafold_path)
