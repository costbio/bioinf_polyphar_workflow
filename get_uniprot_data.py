from bioservices import UniProt
import argparse
import pandas as pd
import numpy as np
import pandas as pd
from pathlib import Path
import re
import math
import os
import time

def checkgenes(gene_path, in_csv):
    # Open the file containing the words to search for
    with open(gene_path, 'r') as file:
        words = [line.strip() for line in file]

    # Open the file to search and read its contents
    with open(in_csv, 'r') as file:
        contents = file.read()

    # Search for the words using regular expressions
    matches = []
    not_found = []
    for word in words:
        pattern = r'\b{}\b'.format(word)
        if re.search(pattern, contents):
            matches.append(word)
        else:
            not_found.append(word)
     
     # Generate the output file path based on the input file name
    input_dir = os.path.dirname(in_csv)
    input_filename = os.path.basename(in_csv)
    output_filename = os.path.splitext(input_filename)[0] + '_not_found.txt'
    output_path = os.path.join(input_dir, output_filename)

     # Save the not found elements to a text file
    with open(output_path, 'w') as file:
        file.write('\n'.join(not_found))

    # Print the matched words and the words that were not found
    if matches:
        print('The following words were found in the file:', ', '.join(matches))
        print(len(matches))
    else:
        print('None of the words were found in the file.')
    if not_found:
        print('The following words were not found in the file:', ', '.join(not_found))
        print(len(not_found))
    else:
        print('All of the words were found in the file.')

def getUniprotData(gene_path, out_csv):
    
    #cwd = os.getcwd()
    #os.chdir(gene_path)
    
    # Read in the txt file containing gene names
    with open(gene_path,"r") as f:
        lines = f.readlines()
        gene_names = [line.rstrip('\n') for line in lines]

    #for filename in os.listdir(gene_path):
    #    genes = os.path.join(gene_path, filename)
    
    #if os.path.isfile(genes):
    #    with open(genes, "r") as file:
    #        gene_names.extend([line.strip() for line in file.readlines()])
        
    if not gene_names:
        raise ValueError(f"No gene names found in directory {gene_path}")    
    
    # Split gene names into blocks of 5
    n_blocks = math.ceil(len(gene_names) / 5)
    gene_blocks = [gene_names[i*5:(i+1)*5] for i in range(n_blocks)]
    
    service = UniProt()
    
    # Query UniProt for each gene block and concatenate the resulting dataframes
    dfs = []
    for block in gene_blocks:
        query = " OR ".join(block).split(" OR ")
        df = service.get_df(query, organism="Homo sapiens")
        df = df[df['Gene Names (primary)'].isin(query)]
        dfs.append(df)
        time.sleep(15)


    
    df= pd.concat(dfs)
    df.to_csv(out_csv)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Query UniProt for gene data')
    parser.add_argument('-gene_file', type=str, help='path to file containing gene names.')
    parser.add_argument('-out_csv',type=str,help='path to output csv file.')

    args = parser.parse_args()

    gene_path = Path(args.gene_file).resolve()

    if not gene_path.is_file():
        raise ValueError(f"Invalid file: {gene_path}")

    out_csv = args.out_csv
    df = getUniprotData(gene_path, out_csv)

    checkgenes(gene_path, out_csv)
    