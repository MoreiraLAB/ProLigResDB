import csv

"""
Variables useful for the bulk of ProLigRes project
"""

__author__ = "C. Marques-Pereira"
__email__ = "amarques@cnc.uc.pt"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "ProLigResDB: A Comprehensive Repository of Experimental Protein Residue-Ligand Interactions from Protein Data Bank"

"""
Folder locations
"""

DEFAULT_LOCATION = "Change/to/user/location/ProLigResDB/"
SCRIPT_FOLDER = DEFAULT_LOCATION + "SCRIPTS/"
INPUT_FOLDER = DEFAULT_LOCATION + "INPUT/"
PDB_FOLDER = INPUT_FOLDER + "PDB/"
CSV_FOLDER = INPUT_FOLDER + "CSV/"
DATA_ANALYSIS_FOLDER = DEFAULT_LOCATION + "DATA_ANALYSIS/"
H5_FOLDER = DEFAULT_LOCATION + "H5_FILES/"
FEATURES_FOLDER = DEFAULT_LOCATION + "FEATURES/"
MORDRED_FOLDER = FEATURES_FOLDER + "MORDRED/"
SPOTONE_FOLDER = FEATURES_FOLDER + "SPOTONE/"
RESOURCES_FOLDER = INPUT_FOLDER + "RESOURCES/"

"""
Header of csv files
"""

HEADER_PDB_SUM = ["pdb", "num_prot_chains", "prot_chains","len_prt_seq", "num_rna_chains", "rna_chains", "num_dna_chains", "dna_chains", "num_ligands", "ligands",
"ligand_type", "hetatm_type", "chains:ligands", "num_nonpolymer_ligands", "nonpolymer_ligands", "num_other_ligands", "other_ligands", "num_peptide_ligands",
"peptide_ligands", "num_saccharide_ligands", "saccharide_ligands", "num_dna_ligands", "dna_ligands", "num_rna_ligands", "rna_ligands"]

HEADER_PDB_INFO = ["Proteins w/ 1 Chain", "Proteins w/ 2 Chains", "Proteins w/ 3 Chains", "Proteins w/ 4 Chains",
"Proteins w/ more than 4 Chains", "Only proteins", "Only DNA", "Only RNA", "Protein + DNA", "Protein + RNA",
"Protein + DNA + RNA", "DNA + RNA", "1 Chain : 1 Ligands", "Chains = Ligands", "Chains < Ligands", "Chains > Ligands",
"Number of Dimers","Number of Tetramers", "Pair Equals", "Ligands of Interest","Nonpolymer Ligands", "Other Ligands", "Peptide Ligands",
"Saccharide Ligands", "DNA Ligands", "RNA Ligands", "Total Ligands", "Structures except Models"]#, "Protein Structures"]

HEADER_PDB_SEQ = ["pdb","chain_sequences", "equal_chains"]

polymer_type_dict = {"non-polymer":1, "other":2,
"d-peptide cooh carboxy terminus":3, "d-peptide nh3 amino terminus":4, "d-peptide linking":5, "d-beta-peptide, c-gamma linking":6, "d-gamma-peptide, c-delta linking":7,
"l-peptide cooh carboxy terminus":8, "l-peptide nh3 amino terminus":9, "l-peptide linking":10, 'l-beta-peptide, c-gamma linking':11, 'l-gamma-peptide, c-delta linking':12,
"peptide-like":13, "peptide linking":14,
"d-saccharide":15, "d-saccharide 1,4 and 1,4 linking":16, "d-saccharide 1,4 and 1,6 linking":17, 'd-saccharide, alpha linking':18, 'd-saccharide, beta linking':19,
"l-saccharide":20, "l-saccharide 1,4 and 1,4 linking":21,"l-saccharide 1,4 and 1,6 linking":22,'l-saccharide, alpha linking':23,'l-saccharide, beta linking':24,
"saccharide":25,
"dna oh 3 prime terminus":26, "dna oh 5 prime terminus":27, "dna linking":28, "l-dna linking":29,
"rna oh 3 prime terminus":30, "rna oh 5 prime terminus":31, "rna linking":32, "l-rna linking":33}

hetatm_type_dict= {"hetain":1, "hetad":2, "atomp":3, "atoms":4, "atomn":5, "?":6}

# SPOTONE variables
HEADER_COLUMNS = ["Input ID","Chain","Residue number","Residue name"]
INTERMEDIATE_SEP = "_"
SYSTEM_SEP = "/"
SEP = ","
COMPLEX_NAME_COL = "Input ID"
CHAIN_NAME_COL = "Chain"
RES_NUMBER_COL = "Residue number"
ENCODING_FILE = RESOURCES_FOLDER + "encoding.csv"
ENCODING_RESIDUE_NAME = "res_letter"
RES_NAME_COL = "Residue name"
THREE_TO_ONE_CODE_CONVERTER = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
ONE_TO_THRE_CODE_CONVERTER = {value:key for key, value in THREE_TO_ONE_CODE_CONVERTER.items()}
AMINO_PROPERTIES_FILE = RESOURCES_FOLDER + "amino_properties.csv"
AMINO_PROPERTIES_THREE_CODE_COL = "three_code"
FEATURES_HOLDER_COLS = ['Helix_propensity', 'Sheet_propensity', 'Helix_propensity_values', 'Sheet_propensity_values', 'MW', 'pka_carboxylate', 'pka_amine', 'pka_side_chain', 'C_number', 'H_number', 'N_number', 'O_number', 'S_number', 'Standard_area_free', 'Standard_area_protein', 'folded_buried_area', 'mean_fractional_area_loss', 'residue_mass', 'monoisotopic_mass']
CENTRED_FEATURES_WINDOWS = [2, 5, 7, 10, 25, 50, 75]
SECTIONS_FEATURES_SPLITS = [100]
PROCESSED_TERMINATION = INTERMEDIATE_SEP + "processed.csv" 
FEATURES_TERMINATION = INTERMEDIATE_SEP + "features.csv"


# Used functions
def open_txt_file(txt_file):
    file = open(txt_file, "r")
    lines = file.readlines()
    list = []
    for line in lines:
        line = line.replace("\n", "")
        list.append(str(line))
    file.close()
    return list

def write_csv(file_name, info_list, header_list):
    with open(file_name, "w", newline = "") as csvfile:
        writer = csv.writer(csvfile, delimiter = ";")
        if info_list[0] not in header_list:
            writer.writerow(header_list)
        for line in info_list:
            if type(line) is not list:
                line = line.split(";")
            else:
                line = [element.replace('"', "") for element in line]
            writer.writerow(line)

def write_txt(path, keys):
    with open(path, "w") as keys_file:
        keys_file.write("\n".join(map(str, keys)))
