#3-build_class.py
#INPUT FILES: (path: "./DATA_ANALYSIS/") "pdb_sum_info.csv" - summarized information for each PDB ID entry
#                                        "ligands_dict_clean.csv" - dictionary with information on selected ligands
#             (path: "./INPUT/CSV/") - Folder with important information from PDB files, saved in CSV format
#OUTPUT FILES: (path: "./H5_FILES/") - Folder with class H5 file
#                                    "class_keys.txt" - txt with class keys
#              (path: "./DATA_ANALYSIS/") "chains_for_fasta.csv" - Interacting Chains

__author__ = "C. Marques-Pereira, A.T. Gaspar"
__email__ = "amarques@cnc.uc.pt"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ProLigResDB: A Comprehensive Repository of Experimental Protein Residue-Ligand Interactions from Protein Data Bank"

from ProLigRes_variables import DATA_ANALYSIS_FOLDER, CSV_FOLDER, H5_FOLDER, write_txt
from scipy.spatial import distance_matrix
import pandas as pd
import os
import h5py
import csv
import numpy as np
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

# Function to select PDB IDS of interest:
# PDBs with at least one protein chain, no DNA or RNA chains, more than one ligand
# other than DNA or RNA ligands, are selected
def ids_of_interest(path_to_results):
    ids_list=[]

    pdb_id=list(pd.read_csv(DATA_ANALYSIS_FOLDER+"pdb_sum_info.csv", sep=";")["pdb"])
    num_prot_chains=pdb_id["num_prot_chains"]
    num_rna_chains=pdb_id["num_rna_chains"]
    num_dna_chains=pdb_id["num_dna_chains"]
    num_ligands=pdb_id["num_ligands"]
    num_dna_ligands=pdb_id["num_dna_ligands"]
    num_rna_ligands=pdb_id["num_rna_ligands"]

    for num in range(len(pdb_id)):
        if int(num_prot_chains[num])>0 and int(num_dna_chains[num])==0 and int(num_rna_chains[num])==0 \
            and int(num_ligands[num])!=0 and int(num_rna_ligands[num])==0 and int(num_dna_ligands[num])==0:
            ids_list.append(pdb_id[num])
        else:
            continue

    return ids_list

# Function to transform csv information into pandas, according to PDB file positions
# pandas columns: "atom_type", "atom_number", "atom", "aa/ligand", "chain", "amino acid_num","x","y","z"
def csv_to_pd(pdb_id):

    with open(CSV_FOLDER+pdb_id.lower()+".csv", 'r') as file:
        csv_id = []
        lines = file.readlines()
        for line in lines:
            line = line.replace("\n","")
            csv_id.append([line[0:4], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46],
            line[46:54]])

    csv_id_pd= pd.DataFrame([i for i in csv_id], columns=["atom_type", "atom_number", "atom", "aa/ligand", \
    "chain", "amino acid_num","x","y","z"])

    return csv_id_pd

# Function to calculate a distance matrix between every protein chain and ligand atom
# Binary classification, according to distance: class 1 if distance <= 5.0, class 0 if distance > 5.0
def subset_pd (pd_pdbfile, interest, id, chains_list, h5_class_w_hetatoms, chain_entries, neg, pos):
    pd_pdbfile[["x", "y", "z"]] = pd_pdbfile[["x", "y", "z"]].apply(pd.to_numeric)
    atoms = pd_pdbfile[pd_pdbfile["atom_type"] == "ATOM"]
    atoms_coordinates = atoms.loc[:,["x","y","z"]]

    ligands = pd_pdbfile[pd_pdbfile["atom_type"].str.startswith("HETA")]
    ligands_coordinates = ligands.loc[:,["x","y","z"]]
    distances = pd.DataFrame(distance_matrix(atoms_coordinates.values, ligands_coordinates.values), \
        index = atoms_coordinates.index, columns = ligands_coordinates.index)

    distances = distances.mask(distances <= 5.0, 1)
    distances = distances.mask(distances > 5.0, 0)

    table = pd.concat((distances, atoms.loc[:,["amino acid_num", "chain"]]), ignore_index=False, axis=1)
    amino_acid_num = table.groupby(table["amino acid_num"], axis=0, as_index=False).max()
    amino_acid_num = amino_acid_num.loc[:,["amino acid_num", "chain"]]
    amino_acid_num.index = amino_acid_num["amino acid_num"].astype("int64")
    
    table = table.groupby(table["amino acid_num"], axis=0).max()
    table.index = table.index.astype("int64")
    new_table = table.iloc[:,:-1].sort_index().transpose()
    new_table = pd.concat((new_table, ligands.loc[:,["aa/ligand"]]), ignore_index=False, axis=1)
    ligand_entries= list(new_table["aa/ligand"].unique())

    chains_interacting= []
    for lig in ligand_entries:
        if lig in interest:
            table_lig = new_table.loc[new_table["aa/ligand"]==lig,:].iloc[:,:-1].transpose()
            table_lig = pd.concat((table.loc[:,["chain"]],table_lig),ignore_index=False, axis=1)

            for chain in chain_entries:
                table_final= table_lig.loc[table_lig["chain"]==chain,:].drop(labels="chain", axis=1)
                amino_acids = amino_acid_num.loc[amino_acid_num["chain"]==chain,:]

                if 1 in table_final.values:
                    chains_interacting.append(chain)

                    class_1column = table_final.transpose()
                    class_1column = pd.concat((class_1column, new_table.loc[new_table["aa/ligand"]==lig,["aa/ligand"]]),\
                    ignore_index=False, axis=1)


                    class_1column = class_1column.groupby(class_1column["aa/ligand"], axis=0).max().transpose()
                    class_1column = pd.concat((amino_acids["amino acid_num"].astype("int64"),class_1column), ignore_index=False, axis=1)
                    class_1column_array = np.array(class_1column)

                    key_id = str(id)+"_"+ str(chain) +":" + str(lig)
                    h5_class_w_hetatoms.create_dataset(key_id, data = class_1column_array)

                    pos+=len(class_1column.loc[class_1column[str(lig)]==1.0])
                    neg+=len(class_1column.loc[class_1column[str(lig)]==0.0])

    if len(chains_interacting) != 0:
        chains_list[id] = list(set(chains_interacting))

    return chains_list, neg, pos, list(h5_class_w_hetatoms.keys())
                
# ids list from CSV file
ids_folder = os.listdir(CSV_FOLDER)
ids_list = [i.replace(".csv","") for i in ids_folder]

print("######## ProLigRes- Build Class ########")
print("- Number of CSV files processed from "+CSV_FOLDER+": "+str(len(ids_list)))

# Create a new folder to save Class H5 file
if not os.path.exists(H5_FOLDER):
    os.makedirs(H5_FOLDER)
    print("- Directory Created: "+H5_FOLDER)
# Create class_5A.hdf5 file
try:
    h5_file = h5py.File(H5_FOLDER+"class_5A.hdf5", "w")
except:
    os.remove(H5_FOLDER+"class_5A.hdf5")
    h5_file = h5py.File(H5_FOLDER+"class_5A.hdf5", "w")

# Information on ligands of interest in "ligands_dict_clean.csv"
ligands_of_interest= pd.read_csv(DATA_ANALYSIS_FOLDER+"ligands_dict_clean.csv", sep=";")
ligands_of_interest= list(ligands_of_interest["ID"])

# empty dictionary and counts
chains_for_fastafile = {}
neg, pos = 0, 0

for id in ids_list:
    # from csv to pandas
    panda = csv_to_pd(id)
    # get PDB chains interacting with ligands of interest
    chains_for_fastafile, neg, pos, class_keys=subset_pd(panda, ligands_of_interest, id, chains_for_fastafile, h5_file, list(set(panda["chain"])), neg, pos)
write_txt(H5_FOLDER+"class_keys.txt", class_keys)

# Save interacting chain information in "chains_for_fasta.csv"
try:
    writer = csv.writer(open(DATA_ANALYSIS_FOLDER+"chains_for_fasta.csv", "w"), delimiter=";")
except:
    os.remove(DATA_ANALYSIS_FOLDER+"chains_for_fasta.csv")
    writer = csv.writer(open(DATA_ANALYSIS_FOLDER+"chains_for_fasta.csv", "w"), delimiter=";")

for key, value in chains_for_fastafile.items():
    writer.writerow([key, value])

print("- H5 file with PDB class saved in "+H5_FOLDER+"class_5A.hdf5")
print("- Total Residues: "+ str(neg+pos))
print("- Negative Residues: "+ str(neg))
print("- Positive Residues: "+ str(pos))
print("- Key classes have PDB ID, chain and ligand information: PDB_CHAIN:LIG")
print("- Interacting chains are saved in "+DATA_ANALYSIS_FOLDER+"chains_for_fasta.csv")
print("- Chain keys saved in "+H5_FOLDER+"class_keys.txt")
print("#######################################")
