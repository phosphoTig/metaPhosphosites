# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 13:53:13 2024

@author: Tigist Tamir, on machine = Luisa
"""
import pandas as pd
import os
from multiprocessing import Pool, cpu_count

# Define the path to your folder containing the CSV files
folder_path = './'

# Define the path to the processed folder
processed_folder = './processed'
os.makedirs(processed_folder, exist_ok=True)

# Read the phosphorylation site information
pSTY = pd.read_csv('PSP_annotated_STY_pdb.csv', sep=',')

# Ensure Full List column is string type
pSTY['Full List'] = pSTY['Full List'].astype(str)

# Create a mapping dictionary from Full List to uID
uid_mapping = dict(zip(pSTY['Full List'], pSTY['uID']))

# Function to process each file
def process_file(filename):
    file_path = os.path.join(folder_path, filename)
    df = pd.read_csv(file_path, sep=',', dtype={'PDB': str, 'res1': str, 'seqnum1': str, 'res2': str, 'seqnum2': str}, low_memory=False)
    
    # Check if required columns are present
    required_columns = {'PDB', 'res1', 'seqnum1', 'res2', 'seqnum2'}
    if not required_columns.issubset(df.columns):
        print(f"File {filename} does not contain all required columns. Skipping.")
        return filename, pd.DataFrame()  # Return filename and an empty DataFrame
    
    # Create the PDBseq and PDBseqTrue columns
    df['PDBseq_res1'] = df['PDB'] + '_' + df['res1'] + '_' + df['seqnum1']
    df['PDBseqTrue_res1'] = df['PDB'] + '_' + df['res1'] + '_' + df['seqnum1']
    df['PDBseq_res2'] = df['PDB'] + '_' + df['res2'] + '_' + df['seqnum2']
    df['PDBseqTrue_res2'] = df['PDB'] + '_' + df['res2'] + '_' + df['seqnum2']
    
    # Fill NaN values with a placeholder to avoid issues during filtering
    df['PDBseq_res1'] = df['PDBseq_res1'].fillna('NA')
    df['PDBseqTrue_res1'] = df['PDBseqTrue_res1'].fillna('NA')
    df['PDBseq_res2'] = df['PDBseq_res2'].fillna('NA')
    df['PDBseqTrue_res2'] = df['PDBseqTrue_res2'].fillna('NA')
    
    # Filter based on PDBseq and PDBseqTrue
    df['isPhosPDBseq_res1'] = df['PDBseq_res1'].apply(lambda x: 1 if x in pSTY['Full List'].values else 0)
    df['isPhosPDBseqTrue_res1'] = df['PDBseqTrue_res1'].apply(lambda x: 1 if x in pSTY['Full List'].values else 0)
    df['isPhosPDBseq_res2'] = df['PDBseq_res2'].apply(lambda x: 1 if x in pSTY['Full List'].values else 0)
    df['isPhosPDBseqTrue_res2'] = df['PDBseqTrue_res2'].apply(lambda x: 1 if x in pSTY['Full List'].values else 0)
    
    filtered_df = df[(df['isPhosPDBseq_res1'] == 1) | (df['isPhosPDBseqTrue_res1'] == 1) | (df['isPhosPDBseq_res2'] == 1) | (df['isPhosPDBseqTrue_res2'] == 1)]
    
    # Map uID from pSTY to the filtered dataframe
    filtered_df['uID'] = filtered_df['PDBseq_res1'].map(uid_mapping).fillna(filtered_df['PDBseqTrue_res1'].map(uid_mapping).fillna(filtered_df['PDBseq_res2'].map(uid_mapping).fillna(filtered_df['PDBseqTrue_res2'].map(uid_mapping))))
    
    return filename, filtered_df

# Get the list of CSV files in the folder
csv_files = [filename for filename in os.listdir(folder_path) if filename.endswith('.csv')]

# Use multiprocessing Pool to process files in parallel
if __name__ == "__main__":
    with Pool(cpu_count()) as pool:
        results = pool.map(process_file, csv_files)

    # Store the processed dataframes in a dictionary
    filtered_dfs_dict = {filename: df for filename, df in results}

    # Concatenate all non-empty filtered dataframes
    non_empty_dfs = [df for df in filtered_dfs_dict.values() if not df.empty]
    final_df = pd.concat(non_empty_dfs, ignore_index=True)

    # Save the final concatenated dataframe to a single CSV file
    final_df.to_csv('filtered_phos_residues.csv', index=False)

    # Save individual filtered dataframes to separate CSV files in the processed folder
    for filename, df in filtered_dfs_dict.items():
        if not df.empty:
            output_file = os.path.join(processed_folder, f'filtered_{filename}')
            df.to_csv(output_file, index=False)
