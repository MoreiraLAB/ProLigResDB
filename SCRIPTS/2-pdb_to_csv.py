#2-pdb_to_csv.py
#INPUT FILES: (path: "./INPUT/PDB/") - Folder with downloaded PDB files
#             (path: "./INPUT/") "ligand.txt" - File with PDB compound dictionary information
#OUTPUT FILES: (path: "./INPUT/CSV/") - Folder with important information from PDB files, saved in CSV format
#              (path: "./DATA_ANALYSIS/") - Folder with important informations about the dataset.
#                                         - "ligands_dict_clean.csv" - dictionary with information on selected ligands
#                                         - "ligand_id.csv" - selected ligands PDB ID
#                                         - "pdb_sum_info.csv" - summarized information for each PDB ID entry
#                                         - "pdb_protein_sequences.csv" - chain sequence information and equal chains within the same PDB ID
#                                         - "pdb_with_models.csv" - PDB IDS excluded due to several models

__author__ = "C. Marques-Pereira, A.T. Gaspar"
__email__ = "amarques@cnc.uc.pt"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ProLigResDB: A Comprehensive Repository of Experimental Protein Residue-Ligand Interactions from Protein Data Bank"

import csv
import os
from ProLigRes_variables import PDB_FOLDER, CSV_FOLDER, DATA_ANALYSIS_FOLDER, HEADER_PDB_SUM, HEADER_PDB_SEQ
import functions_pdb_to_csv as funct

# Create a new folder to save CSV files
if not os.path.exists(CSV_FOLDER):
    os.makedirs(CSV_FOLDER)
    print("- Directory Created: "+CSV_FOLDER)
# Create a new folder to save Data Analysis files
if not os.path.exists(DATA_ANALYSIS_FOLDER):
    os.makedirs(DATA_ANALYSIS_FOLDER)
    print("- Directory Created: "+DATA_ANALYSIS_FOLDER)

# Write Data informatio in files "pdb_sum_info.csv", "pdb_protein_sequences.csv" and "pdb_with_models.csv"
pdb_sum = open(DATA_ANALYSIS_FOLDER+"pdb_sum_info.csv", "w")
pdb_sum_writer = csv.writer(pdb_sum, delimiter=";")
pdb_sum_writer.writerow(HEADER_PDB_SUM)

sequence = open(DATA_ANALYSIS_FOLDER+"pdb_protein_sequences.csv", "w")
sequence_writer=csv.writer(sequence, delimiter=";")
sequence_writer.writerow(HEADER_PDB_SEQ)

models = open(DATA_ANALYSIS_FOLDER+"pdb_with_models.csv", "w")

# empty dictionaries and lists
ligands_dict, dna_rna_ligands, equal_chains = {}, {}, []
percentage =["0%"]*26

# list files in PDB folder
pdb_folder = os.listdir(PDB_FOLDER)

for file in pdb_folder:
    id = file.replace("pdb","").replace(".ent", "")
    # file is the PDB file
    file = open(PDB_FOLDER+file, "r", encoding = "utf-8").readlines()

    # result is the important information retrieved from PDB file: [0]:pdb_information, [1]:sequence, [2]:ligands, [3]:ligands_dna_rna
    result = funct.pdb_to_csv(file, id, models, ligands_dict, dna_rna_ligands)
    # ligands_dict is updated with ligands information
    ligands_dict.update(result[2])
    # dna_rna_ligands is updated with dna and rna ligands
    dna_rna_ligands.update(result[3])
    # equal_chains has the information of same sequence chains within each PDB
    equal_chains = funct.equal_chains_results(result[1], result[0][2])
    # "pdb_protein_sequences.csv" will have the PDB ID, protein sequence and equal chain information
    sequence_writer.writerow([id, result[1], equal_chains])
    # "pdb_sum_info.csv" will have a summary of PDB information
    pdb_sum_writer.writerow(result[0])

print("- Summarized information on PDB dataset stored in "+DATA_ANALYSIS_FOLDER+"pdb_sum_info.csv")
print("- PDB protein sequence and equal chain information stored in "+DATA_ANALYSIS_FOLDER+"pdb_protein_sequences.csv")
print("- PDB files information with several models stored in "+DATA_ANALYSIS_FOLDER+"pdb_with_models.csv")
print("#######################################")
