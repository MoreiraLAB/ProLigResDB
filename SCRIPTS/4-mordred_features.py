# 4-mordred_features.py
# INPUT FILES: (path: "./DATA_ANALYSIS/") "ligands_dict_clean.csv" - dictionary with information on selected ligands
# OUTPUT FILES: (path: "./FEATURES/") - Folder with descriptors and feature information
#                                      - "mordred_descriptors.txt" - List of mordred descriptors
#                                      - "mordred_error.csv" - IDs without mordred features due to error
#              (path: "./H5_FILES/")- "mordred.hdf5" - H5 file with Mordred descriptors, keys correspond to PDB ligand ID
#                                   - "mordred_keys.txt" - txt with mordred keys
#              (path: "./FEATURES/MORDRED/") - Folder with Mordred descriptors for all ligand IDs, in .csv

__author__ = "C. Marques-Pereira, A.T. Gaspar"
__email__ = "amarques@cnc.uc.pt"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ProLigResDB: A Comprehensive Repository of Experimental Protein Residue-Ligand Interactions from Protein Data Bank"

from ProLigRes_variables import DATA_ANALYSIS_FOLDER, MORDRED_FOLDER, FEATURES_FOLDER, H5_FOLDER, write_txt
import pandas as pd
import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors
import csv
import os
import h5py
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

print("######## ProLigRes- Mordred Features ########")
# Create a new folder to save Mordred descriptors
if not os.path.exists(MORDRED_FOLDER):
    os.makedirs(MORDRED_FOLDER)
    print("- Directory Created: " + MORDRED_FOLDER)

# Function to calculate Mordred descriptors and save them in MORDRED folder
# Descriptors are also saved in an H5 file ("mordred.h5py") in FEATURES folder
def drug_feature_extraction(input_ids, input_smiles, path_to_ligands_features):
    calc = Calculator(descriptors, ignore_3D=True)

    mordred_error = open(FEATURES_FOLDER + "mordred_error.csv", "w")
    mordred_error = csv.writer(mordred_error, delimiter=";")
    keys_written = False
    counter = 0

    if not os.path.exists(str(path_to_ligands_features + input_ids[counter] + '.csv')):
        for smile in input_smiles:
            mol = Chem.MolFromSmiles(smile)

            if mol is None:
                mordred_error.writerow([input_ids[counter], smile])
            else:
                feature = calc(mol)
                feature_values = [str(value) if isinstance(value, (str, object)) else value for value in feature.values()]
                output_file = open(str(path_to_ligands_features + input_ids[counter] + '.csv'), "w")

                writer = csv.writer(output_file, delimiter=";")

                writer.writerow(feature.keys())
                writer.writerow(feature.values())
                output_file.close()
                h5_file.create_dataset(input_ids[counter], data=feature_values)
                counter +=1
                if not keys_written:
                    keys_written = True
                    with open(FEATURES_FOLDER + "mordred_descriptors.txt", "w") as keys_file:
                        keys_file.write("\n".join(map(str, feature.keys())))
                    print("- Total number of Mordred features: ", len(list(feature.keys())))
    write_txt(H5_FOLDER + "mordred_keys.txt", list(h5_file.keys()))

if os.path.exists(H5_FOLDER + "mordred.hdf5"):
    h5_file = h5py.File(H5_FOLDER + "mordred.hdf5", "a")
else:
    h5_file = h5py.File(H5_FOLDER + "mordred.hdf5", "w")

# ligand ids and smile list from clean dictionary with ligands of interest
ids = list(pd.read_csv(DATA_ANALYSIS_FOLDER + "ligands_dict_clean.csv", sep=";")["ID"])
smiles = list(pd.read_csv(DATA_ANALYSIS_FOLDER + "ligands_dict_clean.csv", sep=";")["smiles"])

# extract Mordred features for each ligand id
drug_feature_extraction(ids, smiles, MORDRED_FOLDER)

print("- H5 file with Mordred descriptors saved in " + FEATURES_FOLDER + "mordred.hdf5")
print("- Total Ligands: " + str(len(ids)))
print("- Mordred keys are the ligand PDB ID")
print("#######################################")
