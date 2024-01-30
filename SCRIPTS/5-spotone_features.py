#5-spotone_features.py
#INPUT FILES: (path: "./DATA_ANALYSIS/") "pdb_sum_info.csv" - summarized information for each PDB ID entry
#                                        "pdb_protein_sequences.csv" - chain sequence information and equal chains within the same PDB ID
#OUTPUT FILES: (path: "./FEATURES/") - "proteins.fasta" - FASTA file with interacting protein chains
#              (path: "./FEATURES/SPOTONE/") - Creates SPOTONE folder with individual fasta for each chain and spotone features
#              (path: "./H5_FILES/") - "spotone.hdf5" - H5 file with spotone descriptors, key is PDB:Chain
#                                    - "spotone_keys.txt" - txt with spotone keys

__author__ = "C. Marques-Pereira, A.T. Gaspar"
__email__ = "amarques@cnc.uc.pt"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ProLigResDB: A Comprehensive Repository of Experimental Protein Residue-Ligand Interactions from Protein Data Bank"

from ProLigRes_variables import DATA_ANALYSIS_FOLDER, FEATURES_FOLDER, SPOTONE_FOLDER, H5_FOLDER, write_txt
import pandas as pd
import numpy as np
import time
from ast import literal_eval
import os
import csv
import h5py

print("######## ProLigRes- SPOTONE Features ########")

# Create a new folder to save Mordred descriptors
if not os.path.exists(SPOTONE_FOLDER):
    os.makedirs(SPOTONE_FOLDER)
    print("- Directory Created: "+SPOTONE_FOLDER)

def split_fasta(fasta_file, path_to_fasta):
    with open(fasta_file, "r") as file:
        fasta_content = file.read().splitlines()
    sequence = ""
    for line in fasta_content:
        line = line.strip()
        if line.startswith(">"):
            if sequence:
                sequence_name = line[1:] # Remove ">" from the header
                output_file_path = os.path.join(path_to_fasta, f"fasta_{sequence_name}.fasta")
                with open(output_file_path, "w") as output:
                    output.write(f">{sequence_name}\n")
                    output.write(sequence + "\n")
                sequence = ""  # Reset the sequence for the next header
        else:
            sequence += line

def get_entries(input):
    proteins_fasta=open(input,"r").readlines()
    output = [element.replace(">","").replace("\n","") for element in proteins_fasta if ">" in element]
    return output

def join_spotone_results(file, list, input_names, h5_file, keys_written):
	for element in list:
		dataset = file.loc[file["Input ID"] == element.split(":")[0]].loc[file["Chain"] == element.split(":")[1]]
		dataset = dataset.drop(columns="Chain")
		dataset_column= dataset.columns
		dataset = np.array(dataset.iloc[:,4:].astype("float64").round(10))
		try:
			index = input_names.index(element)
			dataset_name = input_names[index]
		except:
			dataset_name = [i for i in input_names if element in i][0]
		if dataset_name not in h5_file.keys():
			h5_file.create_dataset(dataset_name, data = dataset)
		if keys_written == False:
			keys_written = True
			with open(FEATURES_FOLDER + "spotone_descriptors.txt", "w") as keys_file:
				keys_file.write("\n".join(map(str, dataset_column)))
			print("- Total number of Spotone features: ", len(dataset_column))
		write_txt(H5_FOLDER+"spotone_keys.txt", h5_file.keys())
	return keys_written

sequences = pd.read_csv(DATA_ANALYSIS_FOLDER+"pdb_protein_sequences.csv", sep=";").loc[:,["pdb","chain_sequences"]]
pdb_sum = pd.read_csv(DATA_ANALYSIS_FOLDER+"pdb_sum_info.csv", sep=";").loc[:,["pdb","prot_chains"]]
equal_chains = pd.read_csv(DATA_ANALYSIS_FOLDER+"pdb_protein_sequences.csv", sep=";").loc[:,["equal_chains"]]

fasta= open(FEATURES_FOLDER+"proteins.fasta", "w")

#chains_for_fasta comes from build_classes.py and has only the chains that interact with ligands of interest.
chains_for_fasta= pd.read_csv(DATA_ANALYSIS_FOLDER+"chains_for_fasta.csv", sep=";", names=["id","chains"])
ids_list = list(chains_for_fasta["id"].unique())

for i in ids_list:
    index_seq = sequences[sequences["pdb"] == i].index
    index_pdb_sum = pdb_sum[pdb_sum["pdb"] == i].index
    chains_interacting = literal_eval(np.array(chains_for_fasta.loc[chains_for_fasta["id"]==i,"chains"])[0])
    sequence = literal_eval(list(sequences.loc[index_seq, "chain_sequences"])[0])
    protein_chains = literal_eval(list(pdb_sum.loc[index_pdb_sum,"prot_chains"])[0])
    pdb_equal_chains = literal_eval(np.array(equal_chains.loc[index_seq, "equal_chains"])[0])

    indexes = []
    elements_to_remove = []
    for element in chains_interacting:
        for pair in pdb_equal_chains:
            if element in pair:
                indexes.append(pdb_equal_chains.index(pair))
                elements_to_remove.append(element)

    diferent_chains = list(set(chains_interacting) - set(elements_to_remove))
    indexes = list(set(indexes))

    if len(diferent_chains) != 0:
        for num in range(len(diferent_chains)):
            index = protein_chains.index(diferent_chains[num])
            if ("X" not in sequence[index]) and (len(sequence[index])>=32):
                    fasta.write(">"+str(i)+":"+str(diferent_chains[num])+"\n")
                    fasta.write(sequence[index]+"\n")

    if len(indexes) != 0:
        for index in indexes:
            equals = pdb_equal_chains[index][0]
            chain_index = protein_chains.index(equals)
            if ("X" not in sequence[chain_index]) and (len(sequence[chain_index])>=32):
                fasta.write(">"+str(i)+":"+str(pdb_equal_chains[index]).replace(",","_").replace("'","").replace("[","").replace("]","").replace(" ","")+"\n")
                fasta.write(sequence[chain_index]+"\n")

fasta.close()

split_fasta(FEATURES_FOLDER+"proteins.fasta", SPOTONE_FOLDER)

if os.path.exists(H5_FOLDER + "spotone.hdf5"):
    h5_spotone = h5py.File(H5_FOLDER+"spotone.hdf5", "a")
else:
	h5_spotone= h5py.File(H5_FOLDER+"spotone.hdf5", "w")

list_of_files = [i for i in os.listdir(SPOTONE_FOLDER) if i.endswith(".fasta")]
python_raw_command = "nohup python3 -u spotone.py "

for file in list_of_files:
	if not os.path.exists(SPOTONE_FOLDER+file+"_features.csv"):
	    current_command = python_raw_command + str(file) + " &"
	    os.system(current_command)

files = os.listdir(SPOTONE_FOLDER)
names = get_entries(FEATURES_FOLDER+"proteins.fasta")

keys_written = False
for file in files:
    if "fasta_features.csv" in file:
        file = pd.read_csv(SPOTONE_FOLDER+file, sep=",")
        list = file['Input ID'] + ':' + file['Chain']
        keys_written = join_spotone_results(file, list.unique(), names, h5_spotone, keys_written)

print("- H5 file with SPOTONE descriptors saved in "+FEATURES_FOLDER+"spotone.hdf5")
print("- Total Chains: "+ str(len(list_of_files)))
print("- SPOTONE keys are the PDB ID and protein chains: PDB:Chains")
print("#######################################")
