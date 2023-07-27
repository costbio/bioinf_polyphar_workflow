from bioservices import UniProt
import argparse
import pandas as pd
import numpy as np
import os
import pandas as pd
from pathlib import Path
from Bio.PDB import *
from Bio.PDB import PDBList
from tqdm import tqdm
from Bio import PDB
from prody import *
import glob
import warnings
warnings.filterwarnings("ignore")

def download_pdb(df, out_path):
    
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    
    pdb_rows = []
    data_list = []  # Move this line outside the loop

    df = pd.read_csv(df)

    for i, row in df.iterrows():
        pdb_value = row['PDB']
        if isinstance(pdb_value, str):
            pdb_list = pdb_value.split(';')
            new_rows = [[row['Entry']]*len(pdb_list),
                        [row['Entry Name']]*len(pdb_list), pdb_list]
            new_df = pd.DataFrame(new_rows).transpose()
            new_df.columns = ["UniProt ID", "Gene Name", "PDB"]
            pdb_rows.append(new_df) 
            pdbl = PDBList()
            for pdb in pdb_list:
                gene = row['Entry'].split(';')[0]
                gene_dir = os.path.join(out_path, gene)
                if not os.path.exists(gene_dir):
                    os.mkdir(gene_dir)   
                pdbl.retrieve_pdb_file(pdb, pdir=gene_dir, file_format='mmCif')
                pdb_file = os.path.join(gene_dir,pdb + ".cif") 
                if not os.path.exists(pdb_file):
                    pdb_file = os.path.join(gene_dir, pdb + ".cif")
                parser = MMCIFParser()
                for filename in os.listdir(gene_dir):
                    if filename.endswith(".cif"):
                        pdb_id = filename[:-4]  # Remove the ".pdb" extension
                        pdb_file = os.path.join(gene_dir, filename)
                        data = parser.get_structure(pdb_id, pdb_file)
                        header_dict = {"PDB": data.header.get("idcode"),
                                       #"Structure name": data.header.get("name"),
                                       "Resolution": data.header.get("resolution"),
                                       "Has missing residues": data.header.get("has_missing_residues"),
                                       "Structure method": data.header.get("structure_method"),}
                        data_list.append(header_dict)  # Append dictionary to data_list
         
            for dirpath, dirnames, filenames in os.walk(gene_dir):
                for filename in filenames:
                    if filename.endswith(".cif") and filename.startswith("pdb"):
                        old_file_name = os.path.join(dirpath, filename)
                        new_file_name = os.path.join(dirpath, filename[3:])
                        os.rename(old_file_name, new_file_name)
    
    b = pd.concat(pdb_rows)
    res = pd.DataFrame(data_list)  # Use data_list to create res dataframe
    merged_df = pd.merge(b, res, on='PDB', how='outer')
    merged_df = merged_df.drop_duplicates()
    merged_df = merged_df.fillna('-')
    #merged_df = merged_df.dropna()
    merged_df = merged_df.dropna(subset=["UniProt ID", "Gene Name", "PDB", "Has missing residues", "Structure method"])
    merged_df.to_csv(os.path.join(out_path,'Protein_info.csv'))
    return merged_df

def getAlphaFold(merged_df, out_path):
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    alphafold_ID = 'AF-P28335-F1'
    database_version = 'v2'
    new_values = set(merged_df['UniProt ID'])
    new_IDs = []
    
    # Loop over the unique new values and generate a new ID and folder for each value
    for new_value in new_values:
        # Replace the middle part of the alphafold_ID with the new value
        prefix, middle, suffix = alphafold_ID.split('-')
        middle = new_value
        new_ID = '-'.join([prefix, middle, suffix])
        
        # Create a folder with the name of the middle value inside the output directory
        new_folder = os.path.join(out_path, middle)
        if not os.path.exists(new_folder):
            os.mkdir(new_folder)
        
        # Append the new ID to the list
        new_IDs.append(new_ID)
        
        
    for alphafold_ID in new_IDs:
        # Generate the model and error URLs for the current ID
        model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_ID}-model_{database_version}.pdb'
        error_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_ID}-predicted_aligned_error_{database_version}.json'

        # Create a folder with the name of the middle value inside the output directory
        middle = alphafold_ID.split('-')[1]
        middle_folder = os.path.join(out_path, middle)
        if not os.path.exists(middle_folder):
            os.mkdir(middle_folder)
            
        # Download the model and error files for the current ID and save them in the appropriate folder
        model_path = os.path.join(middle_folder, f'{alphafold_ID}.pdb')
        error_path = os.path.join(middle_folder, f'{alphafold_ID}.json')
        os.system(f'curl {model_url} -o {model_path}')
        os.system(f'curl {error_url} -o {error_path}')
        
            
        data_list = []
        for protein_name in new_values:
            protein_folder = os.path.join(out_path, protein_name)
            id_to_gene = dict(zip(merged_df['UniProt ID'], merged_df['Gene Name']))
            gene_name = id_to_gene.get(protein_name)
            for filename in os.listdir(protein_folder):
                if filename.endswith(".pdb"):
                    pdb_id = filename[:-4]
                    pdb_file = os.path.join(protein_folder, filename)
                    parser = PDBParser(PERMISSIVE=True, QUIET=True)
                    data = parser.get_structure(pdb_id, pdb_file)
                    header_dict = {"UniProt ID": protein_name,
                                   "Gene Name": gene_name,
                                   "PDB": data.header.get("idcode") or "-",
                                   "Resolution": data.header.get("resolution"),
                                   "Has missing residues": data.header.get("has_missing_residues"),
                                   "Structure method": data.header.get("structure_method")}
                    data_list.append(header_dict)


    a = pd.DataFrame(data_list)
    combined_df = pd.concat([a, merged_df], ignore_index=True)
    combined_df['AlphaFold'] = [False if val != '-' else True for val in combined_df['PDB']]
    combined_df= combined_df.reindex(columns=['UniProt ID', "Gene Name", 'PDB',"AlphaFold","Resolution","Has missing residues","Structure method"])
    return combined_df

def chains(combined_df, out_path):
    # Create a PDB parser object
    parser = MMCIFParser()

    # Initialize a list to store the data
    data = []

    # Iterate over each subfolder in the directory
    for subfolder in os.listdir(out_path):
        subfolder_path = os.path.join(out_path, subfolder)
        if os.path.isdir(subfolder_path):
            # Iterate over each PDB file in the subfolder
            for pdb_file in os.listdir(subfolder_path):
                if pdb_file.endswith(".cif"):
                    # Load the protein structure from the PDB file
                    pdb_id = os.path.splitext(pdb_file)[0]
                    structure = parser.get_structure(pdb_id, os.path.join(subfolder_path, pdb_file))
                
                # Iterate over each chain in the structure and append its ID and molecule name to the list
                for chain in structure.get_chains():
                    data.append({"PDB": pdb_id.upper(), "Chain ID": chain.id})

    # Create a pandas DataFrame from the list of data
    df = pd.DataFrame(data)
    final_df=pd.merge(combined_df, df, on='PDB', how='outer')
    final_df = final_df.fillna('-')
    final_df = final_df.drop_duplicates()
    final_df.to_csv(os.path.join(out_path,'Protein_info.csv'))
    return final_df

def chain_info(final_df):
    df_pdb_uniprot_chain = pd.read_csv('pdb_chain_uniprot.tsv',delimiter='\t',header=1)
    ids = final_df['Entry'].tolist()   # replace 'column_name' with the name of the column containing the IDs
    df_pdb_af = df_pdb_uniprot_chain[df_pdb_uniprot_chain['SP_PRIMARY'].isin(ids)]
    df_pdb_af.to_csv(os.path.join(out_path,'Protein_info.csv'))
    return df_pdb_af

def extract_required_chains(final_df, out_path):
    ent_files = glob.glob(out_path+'/**/*.cif', recursive=True)

    for ent_file in ent_files:
        try:
            pdb_code = os.path.splitext(os.path.basename(ent_file))[0]
            syst = parseMMCIF(ent_file)

            final_df_pdb = final_df[final_df['PDB'] == pdb_code.upper()]
            chain = final_df_pdb['Chain ID'].values[0]
            uniprot_id = final_df_pdb['UniProt ID'].values[0]

            syst_chain = syst.select('chain %s' % chain)
            writePDB(os.path.join(out_path,uniprot_id,pdb_code+'_'+chain+'.pdb'), syst_chain)
        except Exception as e:
            print(f"An error occurred while processing file {ent_file}: {e}")
            # delete the error file
            os.remove(ent_file)
            # save the name of the deleted file in a text file
            with open('deleted_files.txt', 'a') as f:
                f.write(os.path.basename(ent_file) + '\n')
            continue
                             
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Query PDB and AF2 structure database for protein structures of genes found in input dataframe')
    parser.add_argument('-df_uniprot', type=str, help='path to data frame containing gene names.')
    parser.add_argument('-out_path',type=str,help='path to output folder path.')

    args = parser.parse_args()

    df_uniprot = Path(args.df_uniprot).resolve()

    if not df_uniprot.is_file():
        raise ValueError(f"Invalid file: {df_uniprot}")

    out_path = args.out_path

    merged_df = download_pdb(df_uniprot, out_path)
    print('Downloaded structures from PDB.')
    combined_df = getAlphaFold(merged_df, out_path)
    print('Downloaded structures from AF2 database.')
    final_df = chains(combined_df, out_path)
    extract_required_chains(final_df, out_path)
    print('Extracted required chains from PDB files.')
    print('Done.')


    