from ProLigRes_variables import open_txt_file, write_csv, polymer_type_dict, hetatm_type_dict, INPUT_FOLDER, DATA_ANALYSIS_FOLDER, LIGAND_INPUT
import os

if not os.path.exists(DATA_ANALYSIS_FOLDER):
    os.makedirs(DATA_ANALYSIS_FOLDER)
    print("- Directory Created: "+DATA_ANALYSIS_FOLDER)

"""
info_ligands_file function takes a list from ligands lines and returns _chem_comp parameters that are useful, explicit in conditions
Types of ligands (column _type from ligands[2])
    D-peptide COOH carboxy terminus (1), D-peptide NH3 amino terminus (2), D-peptide linking (3),
    D-saccharide (4), D-saccharide 1,4 and 1,4 linking (5), D-saccharide 1,4 and 1,6 linking (6),
    DNA OH 3 prime terminus, DNA OH 5 prime terminus, DNA linking
    L-peptide COOH carboxy terminus (7), L-peptide NH3 amino terminus (8), L-peptide linking (9),
    L-saccharide (10), L-saccharide 1,4 and 1,4 linking (11), L-saccharide 1,4 and 1,6 linking (12)
    RNA OH 3 prime terminus, RNA OH 5 prime terminus, RNA linking
    non-polymer (13)
    other (14)
    saccharide (15)
"""
def info_ligands_file(file_lines, conditions):
    new_conditions = conditions.split(",")
    ligands_list = []
    comp_id = ""
    ligand_name = ""
    comp_type = ""
    comp_pdbx = ""
    smiles = ""
    counter = 0
    for line in file_lines:
        counter+=1
        if new_conditions[0] in line:
            comp_id= str(line[49:].strip())
            smiles_index = comp_id+" "+new_conditions[4]

        elif new_conditions[1] in line:
            if " ION" or " ion" not in line[49:]:
                ligand_name = str(line[49:].lower().strip())
                if ligand_name=="":
                    for num in range (counter, len(file_lines)):
                        if new_conditions[2] not in file_lines[num]:
                            ligand_name += file_lines[num]
                            continue
                        else:
                            break

        elif new_conditions[2] in line:
            comp_type = str(line[49:].lower().strip())


        elif new_conditions[3] in line:
            comp_pdbx = line[49:].lower().strip()

        elif new_conditions[4] in line:
            smiles = line[len(smiles_index)-1:].replace("\n", " ")
            for num in range(counter, len(file_lines)-1):
                if "SMILES_CANONICAL CACTVS" in file_lines[num]:
                    break
                elif 'SMILES_CANONICAL "OpenEye OEToolkits"' not in file_lines[num] and ('SMILES' not in file_lines[num]):
                    smiles += file_lines[num].replace("\n", " ")
                else:
                    break

            smiles = smiles.replace(" ", ",").replace('"', "").replace(";","").split(",")
            smiles = [i for i in smiles if i != ""]
            smiles = smiles[-1]

        elif comp_id != "" and ligand_name != "" and comp_type != "" and comp_pdbx != "" and smiles != "":
            ligands_list.append([comp_id, ligand_name, comp_type, comp_pdbx, smiles])
            comp_id = ""
            ligand_name = ""
            comp_type = ""
            comp_pdbx = ""
            smiles = ""

    return ligands_list

"""
add_info_in_list function takes a list, an index and a condition and returns the
ID information from lines with the detailed conditions
"""
def add_info_in_list (ligands_list, index, conditions):
    new_conditions = conditions.split(",")
    #black_list = ["dna oh 3 prime terminus", "dna oh 5 prime terminus", "dna linking", "rna oh 3 prime terminus", "rna oh 5 prime terminus", "rna linking" ]
    new_list = []
    if index != 0:
        for line in ligands_list:
            polymer_type=line[index-1].replace('"', "")
            hetatm_type=line[index].replace('"', "")
            if (hetatm_type in new_conditions): #and (polymer_type not in black_list):
                new_list.append((line[0], polymer_type_dict[polymer_type],hetatm_type_dict[hetatm_type]))
    else:
        new_list = [element[index] for element in ligands_list]

    return new_list

"""
ligands is a variable with all lines from ligand dictionary (= ligands_txt), opened with open_txt_file function from functions_clean_ligands_list script
ligand_clean saves the information (ID, name, type and pdbx_type) from compounds on ligands variable using the info_ligands_file function from functions_clean_ligands_list script
ligand_clean_dict.csv will be saved with the information from ligand_clean using the function write_csv from functions_clean_ligands_list
"""
ligands = open_txt_file(INPUT_FOLDER+LIGAND_INPUT)
ligand_clean = info_ligands_file(ligands, '_chem_comp.id,_chem_comp.name,_chem_comp.type,_chem_comp.pdbx_type,SMILES_CANONICAL "OpenEye OEToolkits"')
write_csv(DATA_ANALYSIS_FOLDER+"ligands_dict_clean.csv", ligand_clean, ["ID", "Name", "polymer_type", "hetatm_type", "smiles"])

"""
lig_residue_id is a list with the IDs of all ligands, using the function add_info_in_list from functions_clean_ligands_list
ligand_id.csv will be saved with the information of ligands IDs using write_csv function from functions_clean_ligands_list
good_drug_id saves the IDs of ligands that are hetain and hetad in pdbx_type column, using add_info_in_list function from functions_clean_ligands_list
other_ids saves the IDs of ligands that are atomp and atoms in pdbx_type column, using add_info_in_list function from functions_clean_ligands_list
ids_unidentified saves the IDs of ligands that are "?" in pdbx_type column, using add_info_in_list function from functions_clean_ligands_list
bad_drugs_id is a list with 275 ligands of no interest. This list will be used in pdb_to_csv script to select the lines of interest from each PDB protein
"""
lig_residue_id = add_info_in_list(ligand_clean, 0, "")
write_csv(DATA_ANALYSIS_FOLDER+"ligand_id.csv", lig_residue_id, ["ID"])
comp_interest_id = add_info_in_list(ligand_clean, 3, "hetain,hetad,atomp,atoms,?,atomn")
comp_not_interest_id = [i for i in lig_residue_id if i not in [j[0] for j in comp_interest_id]]

print("######## ProLigRes- PDB to CSV ########")
print("- Clean Ligand Dictionary with Ligands of interest stored in "+DATA_ANALYSIS_FOLDER+"ligands_dict_clean.csv")
print("- Dictionary IDS of all Ligands of interest stored in "+DATA_ANALYSIS_FOLDER+"ligand_id.csv")
