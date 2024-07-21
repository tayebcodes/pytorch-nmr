
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder

import numpy as np
from sklearn.preprocessing import OneHotEncoder
from pprint import pprint
import pickle

# List of standard three-letter amino acid codes
amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
               'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
               'THR', 'TRP', 'TYR', 'VAL']

def amino_acid_to_onehot(amino_acid_code):
    """
    Convert a single three-letter amino acid code to a one-hot encoded vector.

    Parameters:
    amino_acid_code (str): A single three-letter amino acid code.

    Returns:
    np.ndarray: One-hot encoded vector of the amino acid.
    """
    # Ensure the amino acid code is in the standard list
    if amino_acid_code not in amino_acids:
        raise ValueError(f"{amino_acid_code} is not a valid three-letter amino acid code")

    # Initialize OneHotEncoder
    onehot_encoder = OneHotEncoder(categories=[amino_acids], sparse_output=False)

    # Fit and transform the data
    onehot_encoded = onehot_encoder.fit_transform([[amino_acid_code]])

    return onehot_encoded[0]


def convert_angle_to_features(angle):
    """
    Convert a single angle in radians to sine and cosine values for machine learning features.

    Parameters:
    angle (float): An angle in radians.

    Returns:
    np.ndarray: A 1D array with two elements, where the first element is the sine value
                and the second element is the cosine value of the input angle.
    """
    sin_value = np.sin(angle)
    cos_value = np.cos(angle)
    features = np.array([sin_value, cos_value])
    return features


def prep_data(csv_filename):
    X_data = []
    y_data = []
    df = pd.read_csv(csv_filename)
    for index, row in df.iterrows():
        print(index)
        y_current = [
            df.loc[index, 'H_shift'],
            df.loc[index, 'N_shift'],
            df.loc[index, 'CO_shift'],
            df.loc[index, 'CA_shift'],
        ]
        
        current_psi_feature = convert_angle_to_features(df.loc[index, 'current_psi'])
        current_phi_feature = convert_angle_to_features(df.loc[index, 'current_phi'])
        preceding_psi_feature = convert_angle_to_features(df.loc[index, 'preceding_psi'])
        preceding_phi_feature = convert_angle_to_features(df.loc[index, 'preceding_phi'])
        following_psi_feature = convert_angle_to_features(df.loc[index, 'following_psi'])
        following_phi_feature = convert_angle_to_features(df.loc[index, 'following_phi'])
        
        X_current = np.concatenate([
            amino_acid_to_onehot(df.loc[index, 'current_residue_type']),
            amino_acid_to_onehot(df.loc[index, 'preceding_residue_type']),
            amino_acid_to_onehot(df.loc[index, 'following_residue_type']),
            current_psi_feature,
            current_phi_feature,
            preceding_psi_feature,
            preceding_phi_feature,
            following_psi_feature,
            following_phi_feature
        ])
        
        X_data.append(X_current)
        y_data.append(np.array(y_current))

    return np.array(X_data), np.array(y_data)




def prep_and_save_data(csv_filename):
    X_data, y_data = prep_data(csv_filename=csv_filename)
    with open('X_data.pkl', 'wb') as f:
        pickle.dump(X_data, f)
    with open('y_data.pkl', 'wb') as f:
        pickle.dump(y_data, f)