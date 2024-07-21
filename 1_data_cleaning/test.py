import pynmrstar
import pandas as pd
from helpers import *
from pprint import pprint
import pickle

# Configuration dictionary to control which parts of the script should run
config = {
    "run_multi_pdb_indices": False,
    "run_multi_peptide_indices": False,
    "run_fetch_pdbs": False,
    "run_save_backbone_shifts": False,
    "run_test_data": False,
    "run_save_angles": False,
    "run_put_data_together": False
}

# Load the CSV files that were downloaded from BMRB along with the NMR-STAR files
directory = './bmrb_entries'
csv_name = 'query_grid.csv'
df = pd.read_csv(csv_name)

# Remove the data that either have multiple PDBs or contain multimers
# In practice, these two result almost the same result but not identical. Both are needed.
if config["run_multi_pdb_indices"]:
    save_list_to_file(get_multi_pdb_indices(csv_name=csv_name), 'multi_pdb_indices.pkl')
    save_list_to_file(get_multi_peptide_indices(directory=directory, csv_file=csv_name), 'multi_peptide_indices.pkl')

# Modify the dataframe to exclude multi-PDB or multimer structures
multi_peptide_indices = load_list_from_file('multi_peptide_indices.pkl')
multi_pdb_indices = load_list_from_file('multi_pdb_indices.pkl')
indices_to_drop = merge_lists(multi_pdb_indices, multi_peptide_indices)
clean_df = df.drop(indices_to_drop).reset_index(drop=True)

# Fetch PDBs for the clean data and store it in the same folder as the NMR-STAR file
if config["run_fetch_pdbs"]:
    fetch_pdbs(clean_df=clean_df)

# Save the extracted backbone chemical shifts
if config["run_save_backbone_shifts"]:
    save_backbone_shifts(directory=directory, clean_df=clean_df)

# Test the data
if config["run_test_data"]:
    bmrb_id = 4984
    shifts = read_shifts(directory=directory, bmrb_id=bmrb_id)
    features = read_features(directory=directory, bmrb_id=bmrb_id)
    pprint(features)
    exit()

# Save features in each directory
if config["run_save_angles"]:
    save_angles(directory=directory, clean_df=clean_df)

# Combine the data resulting in global_features.csv which will be copied to the prep_data folder for further processing
if config["run_put_data_together"]:
    put_data_together(directory=directory, clean_df=clean_df)
