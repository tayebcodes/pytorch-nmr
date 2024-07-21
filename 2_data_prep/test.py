from helpers import *
from pprint import pprint
import pickle
import numpy as np

# Configuration dictionary to control which parts of the script should run
config = {
    "run_prep_and_save_data": False,
    "run_load_X_y_data": True,
    "run_save_npz": True
}

# Copy global_features.csv to this directory for further processing
# This results in X_data.pkl and y_data.pkl
if config["run_prep_and_save_data"]:
    csv_filename = 'global_features.csv'
    prep_and_save_data(csv_filename=csv_filename)

# Load X_data and y_data from pickle files
# Read the X_data.pkl and y_data.pkl data and organize them
if config["run_load_X_y_data"]:
    with open('X_data.pkl', 'rb') as f:
        X_data = pickle.load(f)

    with open('y_data.pkl', 'rb') as f:
        y_data = pickle.load(f)

# Save data in npz format and move them to Google Colab
if config["run_save_npz"]:
    np.savez('data.npz', X_data=X_data, y_data=y_data)
