#1-pdb_download.py
#INPUT FILES: (path: "./INPUT/") "pdb_ids.txt" - list with PDB IDS
#OUTPUT FILES: (path: "./INPUT/PDB/") - Folder with downloaded PDB files

__author__ = "C. Marques-Pereira, A.T. Gaspar"
__email__ = "amarques@cnc.uc.pt"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ProLigResDB: A Comprehensive Repository of Experimental Protein Residue-Ligand Interactions from Protein Data Bank"

from Bio.PDB import *
from ProLigRes_variables import INPUT_FOLDER, PDB_FOLDER
import os

print("######## ProLigRes- PDB download ########")

# Create a new folder to save PDB files
if not os.path.exists(PDB_FOLDER):
    os.makedirs(PDB_FOLDER)
    print("- Directory Created: "+PDB_FOLDER)

# List IDS from a txt file
def open_file(file_path):
    data_list = []
    with open(file_path, "r") as file:
        for line in file:
        # Use strip() to remove newline characters
            line = line.strip()
        # Use split(',') to separate values by commas
            values = line.split(',')
            data_list.append(values)

    return data_list

# Download all pdb files from an id list
def download_pdb(pdb_list):
    pdbl=PDBList()
    pdbl.download_pdb_files(pdb_list, pdir=PDB_FOLDER, obsolete=False, file_format="pdb")

# Get PDB IDS as a list fromo file "pdb_ids.txt"
ids_list = open_file(INPUT_FOLDER+"pdb_ids.txt")[0]
# Dowload all PDB files to PDB folder
pdb_download=download_pdb(ids_list)

print("- PDB ids from "+INPUT_FOLDER+"pdb_ids.txt \ndownloaded to "+PDB_FOLDER +" in .ent format.")
print("#########################################")
