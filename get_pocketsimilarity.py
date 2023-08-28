#get_pocketsimilarity.py1
import pandas as pd
import numpy as np
import subprocess
import os
import glob
import re
import sys
from pathlib import Path
import shutil
import argparse




def getPMsimilarities(df_path,final_path, out_path,pm_path,run_path, pockets_path,drugScore_cutoff=0.5):
    
    df_pockets = pd.read_csv(os.path.join(pockets_path,'DrugScore.csv'))
    df_pockets = df_pockets[df_pockets['drugScore'] >= drugScore_cutoff]
    
    cwd = os.getcwd()
    os.chdir(pm_path)

    pockets = []
    for root, dirs, files in os.walk(pockets_path):
        for index, row in df_pockets.iterrows():
            filename = row["Pocket Filename"]
            if os.path.exists(os.path.join(root, filename)):
                pockets.append(os.path.join(root, filename))

    for i in pockets:
        subprocess.call("cp %s $PWD" %i ,shell=True)
   
    if not os.path.exists(final_path):
        os.mkdir(final_path)


    pockets_bn = [os.path.basename(pocket) for pocket in pockets]
    
    df_sims_list = list()
    for i in range(len(pockets_bn)):
        for j in range(i+1, len(pockets_bn)):
            df_sim_list = list()
            
    os.chdir(run_path)
    cmd = ' '.join(["bash", "Step0-cabbage.sh", "Sample_pockets/"])
    pocket_match = subprocess.call(cmd,shell=True)
    
    os.chdir(out_path)
    subprocess.call("cp cabbage-file_maker/outfile.cabbage $PWD",shell=True)
    subprocess.run(["./Step3-PM_typeA",  "outfile.cabbage"])
    
    for file_path in Path(pm_path).glob('*'):
        try:
            if file_path.is_file():
                file_path.unlink()
        except Exception as e:
            print(f'Error deleting {file_path}: {e}')
    
    with open('PocketMatch_score.txt', 'r') as f:
        valid_lines = []
        bad_lines = []
        for i, line in enumerate(f):
            try:
                fields = line.strip().split(' ')
                # Process the fields here
                # ...
                valid_lines.append(fields)
            except:
                # If there's a parsing error, store the bad line in a separate list
                bad_lines.append(line)
        # Convert the lists to dataframes
        df= pd.DataFrame(valid_lines)

    df = df[[0,1,2,3]]
    df.columns = ['Pocket1','Pocket2','Pmin','Pmax']
    df = df.convert_dtypes()
    df= df[df['Pmax'] != 'NULL____']
    df= df[df['Pmin'] != 'NULL____']
    df["Pmin"] = pd.to_numeric(df["Pmin"])
    df["Pmax"] = pd.to_numeric(df["Pmax"])
    df= df[df['Pmax'] != 1]
    df= df[df['Pmin'] != 1]
    df_orig = df
    #df= df[df['Pocket1'].str[:8] != df['Pocket2'].str[:8]]
    #df[["Pocket1","1","2"]]=df['Pocket1'].str.split("___" , expand=True )
    #delete_col=["1","2"]
    #df.drop(delete_col, axis=1)
    #df[["Pocket2","1","2"]]=df['Pocket2'].str.split("___" , expand=True )
    #delete_col=["1","2"]
    #df=df.drop(delete_col, axis=1)
    
    # Remove trailing underscores from pocket filenames.
    df['Pocket1'] = df['Pocket1'].str.rstrip('_')
    df['Pocket2'] = df['Pocket2'].str.rstrip('_')
    # Define a lambda function that extracts the pocket path from the pockets list based on pocket name
    get_pocket_path = lambda pocket_name: next((pocket for pocket in pockets if pocket_name in pocket), None)
    df['Pocket_Path'] = df['Pocket1'].apply(get_pocket_path)
    df['Pocket_Path_2'] = df['Pocket2'].apply(get_pocket_path)
    col = df.pop('Pocket_Path')
    df.insert(1, 'Pocket_Path', col)
    col = df.pop('Pocket_Path_2')
    df.insert(3, 'Pocket_Path_2', col)

    # Read df_pdb_af from pdb_path
    df_pdb_af = pd.read_csv(os.path.join(df_path,'Protein_info.csv'))


    dfa=df.sort_values(by='Pmax',ascending=False)
    df1 = dfa[['Pocket1']].copy()
    df1['UniProt ID'] = df1['Pocket1'].apply(lambda x: x.split('_', 1)[0] if '_' in x and not x.startswith('AF') else '')
    df1 = dfa[['Pocket1']].copy()
    df1['UniProt ID'] = df1['Pocket1'].apply(lambda x: x.split('_', 1)[0] if '_' in x and not x.startswith('AF') else '')
    df1['Pocket'] = df1['Pocket1'].apply(lambda x: x.split('_', 1)[1] if '_' in x and not x.startswith('AF') else x)
    df1 = df1[['UniProt ID', 'Pocket']]
    dfa = dfa.join(df1)
    dfa = dfa.reindex(columns=['UniProt ID', "Pocket","Pocket_Path","Pocket1","Pocket2","Pocket_Path_2","Pmin","Pmax"])
    dfa=dfa.drop("Pocket1", axis=1)
    df2 = dfa[['Pocket2']].copy()
    df2['UniProt ID_2'] = df2['Pocket2'].apply(lambda x: x.split('_', 1)[0] if '_' in x and not x.startswith('AF') else '')
    df2['Pocket_2'] = df2['Pocket2'].apply(lambda x: x.split('_', 1)[1] if '_' in x and not x.startswith('AF') else x)
    df2 = df2[['UniProt ID_2', 'Pocket_2']]
    dfa = dfa.join(df2)
    col = dfa.pop('UniProt ID_2')
    dfa.insert(3, 'UniProt ID_2', col)
    col = dfa.pop('Pocket_2')
    dfa.insert(4, 'Pocket_2', col)
    dfa=dfa.drop("Pocket2", axis=1)
    a1=dfa[["Protein","Pocket"]]= dfa['Pocket'].str.split(".", n=1, expand=True)
    col = dfa.pop('Protein')
    dfa.insert(1, 'Protein', col)
    dfa['Pocket'] = dfa['Pocket'].str.replace('pdb_res_', '')
    a1=dfa[["Pocket","1"]]= dfa['Pocket'].str.split(".", n=1, expand=True)
    dfa=dfa.drop("1", axis=1)
    dfa['Protein'] = dfa['Protein'].str.upper()
    a1=dfa[["Protein_2","Pocket_2"]]= dfa['Pocket_2'].str.split(".", n=1, expand=True)
    col = dfa.pop('Protein_2')
    dfa.insert(5, 'Protein_2', col)
    dfa['Pocket_2'] = dfa['Pocket_2'].str.replace('pdb_res_', '')
    a1=dfa[["Pocket_2","1"]]= dfa['Pocket_2'].str.split(".", n=1, expand=True)
    dfa=dfa.drop("1", axis=1)
    dfa['Protein_2'] = dfa['Protein_2'].str.upper()
    dfa['P'] = dfa['Protein'].apply(lambda x: ''.join(re.findall(r'[A-Za-z]+\d+', x)[:-1]) if '-' in x else '')
    mask = dfa['UniProt ID'] == ''
    dfa.loc[mask, 'UniProt ID'] = dfa.loc[mask, 'P']
    dfa['P'] = dfa['Protein_2'].apply(lambda x: ''.join(re.findall(r'[A-Za-z]+\d+', x)[:-1]) if '-' in x else '')
    mask = dfa['UniProt ID_2'] == ''
    dfa.loc[mask, 'UniProt ID_2'] = dfa.loc[mask, 'P']
    dfa=dfa.drop("P", axis=1)
    #Create a dictionary that maps gene names to IDs
    gene_to_id = dict(zip(df_pdb_af['UniProt ID'], df_pdb_af['Gene Name']))
    # Add a new column to the DataFrame that contains the corresponding ID for each gene name
    dfa['Gene Name'] = dfa['UniProt ID'].map(gene_to_id)
    col = dfa.pop('Gene Name')
    dfa.insert(1, 'Gene Name', col)
    dfa['Gene Name_2'] = dfa['UniProt ID_2'].map(gene_to_id)
    col = dfa.pop('Gene Name_2')
    dfa.insert(6, 'Gene Name_2', col)
    #avoid getting same genes scores
    dfa= dfa[dfa['UniProt ID']!= dfa['UniProt ID_2']]
    df_sims=dfa
    
    df_sims.to_csv(os.path.join(final_path,'Pocketsim.csv'))
    result=glob.glob(out_path +'/*.txt')
    
    for text_path in result:
        shutil.copy(text_path, final_path)
        os.remove(text_path)
    
    os.chdir(cwd)
    return df_sims




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run getPMsimilarities function.')
    parser.add_argument('-df_path', type=str, help='Path to input DataFrame.')
    parser.add_argument('-final_path', type=str, help='Path to output directory.')
    parser.add_argument('-out_path', type=str, help='Path to PocketMatch directory.')
    parser.add_argument('-pm_path', type=str, help='Path to PocketMatch samplepockets directory.')
    parser.add_argument('-run_path', type=str, help='Path to PocketMatch run directory.')
    parser.add_argument('-pockets_path', type=str, help='Path to pocket files.')
    parser.add_argument('-drugScore_cutoff', type=float, default=0.5, help='DrugScore cutoff.')
    args = parser.parse_args()

    out_path = os.path.abspath(args.out_path)
    pm_path = os.path.abspath(args.pm_path)
    df_path = os.path.abspath(args.df_path)
    final_path = os.path.abspath(args.final_path)
    run_path = os.path.abspath(args.run_path)
    pockets_path = os.path.abspath(args.pockets_path)
    drugScore_cutoff = float(args.drugScore_cutoff)

getPMsimilarities(df_path,final_path, out_path,pm_path,run_path, pockets_path,drugScore_cutoff)