import pynmrstar
import pickle
import pandas as pd
import requests
from Bio import PDB
import glob
import os

amino_acids = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
]


def read_nmr_star(directory, bmrb_id):
    directory = directory if directory.endswith('/') else directory + '/' # add trailing slash to the directory if doesn't exist
    bmrb_id = str(bmrb_id) # make sure id is 
    file_path = directory + 'bmr' + bmrb_id + '/' + 'bmr' + bmrb_id + '_3.str'
    # Load the NMR-STAR file
    entry = pynmrstar.Entry.from_file(file_path)
    return entry


def save_list_to_file(data, filename):
    """
    Save a list to a file using pickle.

    Args:
        data (list): The list to save.
        filename (str): The name of the file to save the list to.
    """
    with open(filename, 'wb') as file:
        pickle.dump(data, file)


def load_list_from_file(filename):
    """
    Load a list from a file using pickle.

    Args:
        filename (str): The name of the file to load the list from.

    Returns:
        list: The loaded list.
    """
    with open(filename, 'rb') as file:
        data = pickle.load(file)
    return data


def get_multi_pdb_indices(csv_name):
    """
    Extract indices of rows with multiple PDB IDs from a CSV file.

    Args:
        csv_name (str): The name of the CSV file to read.

    Returns:
        list: A list of indices where the "pdb_ids" column contains multiple PDB IDs.
    """
    import pandas as pd

    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_name)

    multi_pdb_indices = []
    for index, row in df.iterrows():
        pdb_ids = df.loc[index, "pdb_ids"].split(",")
        if len(pdb_ids) > 1:
            multi_pdb_indices.append(index)

    return multi_pdb_indices


def get_multi_peptide_indices(directory, csv_file):
    """
    Get indices of entries with multiple peptide sequences from a CSV file.

    Args:
        directory (str): The directory where the NMR-STAR files are stored.
        csv_file (str): The CSV file containing the entry IDs.

    Returns:
        list: A list of indices with multiple peptide sequences.
    """
    df = pd.read_csv(csv_file)
    indices = []
    for index, row in df.iterrows():
        #print(index)
        bmrb_id = df.loc[index, 'Entry_ID']
        entry = read_nmr_star(directory, bmrb_id)
        sequence = get_fasta(entry)
        if len(sequence) != 1:
            print(sequence)
            indices.append(index)
    return indices


def merge_lists(list1, list2):
    """
    Merge two lists without duplicates using list comprehension to preserve order.

    Args:
        list1 (list): The first list.
        list2 (list): The second list.

    Returns:
        list: A merged list without duplicates.
    """
    merged_list = list1 + [item for item in list2 if item not in list1]
    return merged_list


def download_pdb(pdb_id, save_path):
    """
    Download a PDB file from the RCSB PDB website.

    Args:
        pdb_id (str): The PDB ID of the file to download.
        save_path (str): The local path where the file should be saved.
    """
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response = requests.get(url)

    if response.status_code == 200:
        with open(save_path, 'wb') as file:
            file.write(response.content)
        print(f'{pdb_id}.pdb has been downloaded and saved to {save_path}')
    else:
        print(f'Failed to download {pdb_id}.pdb')


def fetch_pdbs(clean_df):
    """
    Fetch PDB files for each entry in the DataFrame and save them to specific paths.

    Args:
        clean_df (pd.DataFrame): A DataFrame containing 'pdb_ids' and 'Entry_ID' columns.

    This function iterates over each row in the DataFrame, extracts the PDB ID and BMRB ID,
    constructs the save path, and calls the `download_pdb` function to download the PDB file.
    It prints the BMRB ID and PDB ID for each entry.

    Example:
        clean_df = pd.DataFrame({
            'pdb_ids': ['1abc', '2def'],
            'Entry_ID': [1234, 5678]
        })
        fetch_pdbs(clean_df)
    """
    for index, row in clean_df.iterrows():
        pdb_id = clean_df.loc[index, 'pdb_ids']
        bmr_id = clean_df.loc[index, 'Entry_ID']
        save_path = f'./bmrb_entries/bmr{bmr_id}/{pdb_id}.pdb'
        download_pdb(pdb_id, save_path)
        
        print(f'BMRB ID: {bmr_id}, PDB ID: {pdb_id}')


def get_fasta(entry):
    sequence = []
    for saveframe in entry:
        if saveframe.category == 'entity':
                if saveframe.get_tag('_Entity.Polymer_type')[0] ==  'polypeptide(L)':
                    sequence.append(saveframe.get_tag('_Entity.Polymer_seq_one_letter_code')[0].replace('\n', ''))
    return sequence


def get_backbone_shifts1(entry):
    """
    Extract backbone chemical shifts from an NMR-STAR entry.

    Args:
        entry (pynmrstar): The NMR-STAR entry read by pynmrstar.

    Returns:
        list: A list of backbone chemical shifts. Each element is a list containing:
              [residue_number, residue_type, H_shift, N_shift, C_shift, Ca_shift]
    
    Example:
        directory = './bmrb_entries'
        bmrb_id = '4023'
        entry = read_nmr_star(directory=directory,bmrb_id=bmrb_id)
        sequence = get_fasta(entry)[0]
        length = len(sequence)
        data = get_backbone_shifts(entry,sequence_length=length)
        pprint(data)
    """
    data = []
    fasta = get_fasta(entry)[0]
    for saveframe in entry:
        if saveframe.category == 'assigned_chemical_shifts':
            for residue in range(0,len(fasta)):
                print(residue)
                residue_number = residue + 1  # the residue number in the iteration
                H_shift = None
                N_shift = None
                C_shift = None
                Ca_shift = None  # Ca shift
                current_residue_type = None

                for i in range(len(saveframe.get_tag('_Atom_chem_shift.ID'))):  # iterate over the chemical shift table
                    current_residue_number = int(saveframe.get_tag('_Atom_chem_shift.Seq_ID')[i])

                    if current_residue_number == residue_number:
                        current_residue_type = saveframe.get_tag('_Atom_chem_shift.Comp_ID')[i]
                        atom_id = saveframe.get_tag('_Atom_chem_shift.Atom_ID')[i]
                        shift_val = saveframe.get_tag('_Atom_chem_shift.Val')[i]

                        if atom_id == 'H':
                            H_shift = shift_val
                        elif atom_id == 'N':
                            N_shift = shift_val
                        elif atom_id == 'C':
                            C_shift = shift_val
                        elif atom_id == 'CA':
                            Ca_shift = shift_val

                current = [residue_number, current_residue_type, H_shift, N_shift, C_shift, Ca_shift]
                data.append(current)
    return data


def get_backbone_shifts(entry):
    """
    Extract backbone chemical shifts from an NMR-STAR entry.

    Args:
        entry (pynmrstar.Entry): The NMR-STAR entry read by pynmrstar.

    Returns:
        list: A list of backbone chemical shifts. Each element is a list containing:
              [residue_number, residue_type, H_shift, N_shift, C_shift, Ca_shift]

    """
    data = []
    fasta = get_fasta(entry)[0]
    sequence_length = len(fasta)

    for saveframe in entry:
        if saveframe.category == 'assigned_chemical_shifts':
            # Extract the necessary tags only once
            seq_ids = saveframe.get_tag('_Atom_chem_shift.Seq_ID')
            comp_ids = saveframe.get_tag('_Atom_chem_shift.Comp_ID')
            atom_ids = saveframe.get_tag('_Atom_chem_shift.Atom_ID')
            shift_vals = saveframe.get_tag('_Atom_chem_shift.Val')
            
            # Create a dictionary to store shifts for each residue
            shifts_dict = {i: {'H': None, 'N': None, 'C': None, 'CA': None} for i in range(1, sequence_length + 1)}
            residue_types = {i: None for i in range(1, sequence_length + 1)}

            for i in range(len(seq_ids)):
                residue_number = int(seq_ids[i])
                atom_id = atom_ids[i]
                shift_val = shift_vals[i]
                residue_type = comp_ids[i]
                
                if residue_number in shifts_dict:
                    shifts_dict[residue_number][atom_id] = shift_val
                    residue_types[residue_number] = residue_type
            
            for residue_number in range(1, sequence_length + 1):
                shifts = shifts_dict[residue_number]
                current = [
                    residue_number, 
                    residue_types[residue_number], 
                    shifts['H'], 
                    shifts['N'], 
                    shifts['C'], 
                    shifts['CA']
                ]
                data.append(current)

    return data


def save_backbone_shifts(directory,clean_df):
    for index,row in clean_df.iterrows():

        bmrb_id = clean_df.loc[index, 'Entry_ID']
        print(bmrb_id)
        entry = read_nmr_star(directory=directory,bmrb_id=bmrb_id)
        data = get_backbone_shifts(entry)
        filename = f'{directory}/bmr{bmrb_id}/shifts.pkl'
        save_list_to_file(data=data,filename=filename)

def read_shifts(directory,bmrb_id):
    filename = f'{directory}/bmr{bmrb_id}/shifts.pkl'
    return load_list_from_file(filename=filename)

def read_features(directory,bmrb_id):
    filename = f'{directory}/bmr{bmrb_id}/features.pkl'
    return load_list_from_file(filename=filename)

def read_pdb_and_create_feature_vectors(directory, bmrb_id):
    pattern = os.path.join(directory, f'bmr{bmrb_id}', '*.pdb')
    pdb_files = glob.glob(pattern)
    
    if not pdb_files:
        print(f"No PDB files found for BMRB ID {bmrb_id}")
        return []  # Return an empty list to indicate no feature vectors

    pdb_file = pdb_files[0]
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    model = structure[0]  # Use the first model
    polypeptides = PDB.PPBuilder().build_peptides(model)
    feature_vectors = []

    for poly in polypeptides:
        phi_psi = poly.get_phi_psi_list()
        for i, residue in enumerate(poly):
            residue_number = residue.get_id()[1]
            current_residue_type = residue.get_resname()

            # Get psi and phi angles for the current residue
            current_phi, current_psi = phi_psi[i]

            # Get preceding and following residues if they exist
            preceding_residue_type = poly[i-1].get_resname() if i > 0 else None
            following_residue_type = poly[i+1].get_resname() if i < len(poly) - 1 else None

            # Get psi and phi angles for the preceding residue if it exists
            preceding_phi, preceding_psi = phi_psi[i-1] if i > 0 else (None, None)

            # Get psi and phi angles for the following residue if it exists
            following_phi, following_psi = phi_psi[i+1] if i < len(poly) - 1 else (None, None)

            # Create the feature vector
            feature_vector = [
                residue_number,
                current_residue_type,
                preceding_residue_type,
                following_residue_type,
                current_psi,
                current_phi,
                preceding_psi,
                preceding_phi,
                following_psi,
                following_phi
            ]
            
            feature_vectors.append(feature_vector)

    return feature_vectors


def save_angles(directory, clean_df):
    for index, row in clean_df.iterrows():
        bmrb_id = clean_df.loc[index, 'Entry_ID']
        print(f"index: {index}, bmrb_id: {bmrb_id}")
        filename = os.path.join(directory, f'bmr{bmrb_id}', 'features.pkl')
        
        # Check if features.pkl already exists
        if os.path.exists(filename):
            print(f"File {filename} already exists, skipping...")
            continue
        
        # Calculate features only if the file does not exist
        features = read_pdb_and_create_feature_vectors(directory, bmrb_id)
        
        # Skip saving if features is empty
        if not features:
            print(f"No features to save for BMRB ID {bmrb_id}")
            continue
        
        # Save the features to the file
        save_list_to_file(data=features, filename=filename)
        print(f"Saved features for BMRB ID {bmrb_id} to {filename}")

def put_data_together(directory, clean_df):
    rows = []  # Use a list to collect rows

    for index, row in clean_df.iterrows():
        bmrb_id = clean_df.loc[index, 'Entry_ID']
        print(f"index: {index}, bmrb ID: {bmrb_id}")
        
        try:
            features = read_features(directory=directory, bmrb_id=bmrb_id)
        except FileNotFoundError:
            print(f"No features file found for BMRB ID {bmrb_id}, skipping...")
            continue
        
        try:
            shifts = read_shifts(directory=directory, bmrb_id=bmrb_id)
        except FileNotFoundError:
            print(f"No shifts file found for BMRB ID {bmrb_id}, skipping...")
            continue
        
        # Convert features and shifts to DataFrames if they are not already
        features_df = pd.DataFrame(features)
        shifts_df = pd.DataFrame(shifts)
        
        # Check if the DataFrames are empty
        if features_df.empty:
            print(f"Features DataFrame is empty for BMRB ID {bmrb_id}, skipping...")
            continue
        
        if shifts_df.empty:
            print(f"Shifts DataFrame is empty for BMRB ID {bmrb_id}, skipping...")
            continue
        
        # Drop rows with any None values
        features_df.dropna(inplace=True)
        shifts_df.dropna(inplace=True)

        # Merge data based on conditions
        for _, shift_row in shifts_df.iterrows():
            residue_number = shift_row[0]
            residue_name = shift_row[1]
            if residue_name in amino_acids:
                matched_rows = features_df[(features_df.iloc[:, 0] == residue_number) & (features_df.iloc[:, 1] == residue_name)]
                for _, feature_row in matched_rows.iterrows():
                    new_row = [
                        bmrb_id,
                        residue_number,
                        residue_name,
                        shift_row[2],  # H_shift
                        shift_row[3],  # N_shift
                        shift_row[4],  # CO_shift
                        shift_row[5],  # CA_shift
                        feature_row[1],  # current_residue_type
                        feature_row[2],  # preceding_residue_type
                        feature_row[3],  # following_residue_type
                        feature_row[4],  # current_psi
                        feature_row[5],  # current_phi
                        feature_row[6],  # preceding_psi
                        feature_row[7],  # preceding_phi
                        feature_row[8],  # following_psi
                        feature_row[9]   # following_phi
                    ]
                    rows.append(new_row)
    
    # Create global_df from collected rows
    global_df = pd.DataFrame(rows, columns=[
        'bmrb_id', 'residue_number', 'residue_name', 'H_shift', 'N_shift', 'CO_shift', 'CA_shift', 
        'current_residue_type', 'preceding_residue_type', 'following_residue_type', 
        'current_psi', 'current_phi', 'preceding_psi', 'preceding_phi', 'following_psi', 'following_phi'
    ])

    # Save the global_df to a CSV file on the hard drive
    output_file = 'global_features.csv'
    global_df.to_csv(output_file, index=False)
    print(f"Saved combined data to {output_file}")