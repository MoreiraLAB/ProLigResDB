from ProLigRes_variables import CSV_FOLDER
from compound_dictionary import comp_interest_id

"""
pdb_to_csv function extract the lines beginning with ATOM/HETATM of interest and writes them in a .csv.
Saves all information in pdb_sum.csv (if it has a DNA/RNA/PRT in ATOM, how many chains the protein has, which chains,
how many and which ligands, and if ligands are nonpolymer, other, peptide, saccharide, DNA or RNA)
Saves in sequences.csv a sequence for each protein chain and compares chains
It returns a dictionary that takes a PDB number for key and two lists for values: the first list has the pdb chains
and the second each ligand_namein the pdb.
"""
def pdb_to_csv(input_file, id, model_file, ligands, ligands_dna_rna):
    """
    drugs catches all information of good drugs (type and pdbx_type of ligand)
    pdb_dict saves a dictionary with simple information of each pdb
    pdb_sum.csv will save all informations in each .pdb
    sequence.csv saves all sequences of all protein chains and equal chains for each protein
    big_list will append all important numbers in within every PDBs and percentage list saves the same numbers in percentage
    total_ligands will be how many ligands (in good_drug) are in every PDBs
    aa is a dictionary that correlates the 3 letter code of amino acids with an 1 letter code to retrieve sequences
    identifier keeps all codes for DNA and RNA in ATOM lines
    """
    drugs = [d[0] for d in comp_interest_id]

    aa = {"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "ASX":"B", "CYS":"C", "GLU":"E", "GLN":"Q", "GLX":"Z", "GLY":"G", "HIS":"H", "ILE":"I", "LEU":"L",
    "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V"}
    identifier = ["DT", "DA", "DC", "DG", "DI", "G", "U", "A", "C", "I"]

    """
    each pdb in input_list will be opened and read line by line
    list_of_lines will save every line of interest: ATOM (DNA/RNA/PRT) HETATM(lines with good_drug - from 3 letter code)
    prot_chains will append every protein chains in PDBs with proteins
    dna_chains will append every DNA chain
    rna_chains will append every RNA chain
    ligands will append every 3 letter code from good_drug only
    ligand_types will catch the type of a specific ligand_namefrom drugs information and keep it
    nonpolymer_ligands, other_ligands, peptide_ligands, saccharide_ligands, dna_ligands and rna_ligands will save the information of ligands of these types and keep the code for pdbx_type
    seq, residue and chain begin empty. seq will append the sequence for each chain, residue will correspond to the previous aminoacid residue and the chain to the chain of the previous aminoacid
    sequence will append all sequences of all protein chains
    equal_chains will append all equal chains for each protein
    writer will write every omportant line in a new .csv for each interesting PDB
    """

    csvfile = open(CSV_FOLDER+id+".csv", "w")

    ligands_list, prot_chains, dna_chains, rna_chains =[], [], [], []
    ligand_types, nonpolymer_ligands, other_ligands, peptide_ligands=[], [], [], []
    saccharide_ligands, dna_ligands, rna_ligands, sequence =[], [], [], []
    seq, residue_number, chain ="","",""

    """
    each line of PDB will be read and only the lines starting with either "ATOM" or "HETATM" will be analysed and appended to list_of_lines
    if the line begins with ATOM, residue_name can be a PRT, DNA or RNA atom. If residue_name is in identifier list, the atom is DNA or RNA and will append DNA or RNA chains in dna_chains and rna_chains, respectively
    if residue_name is a PRT, the program will save the chain in prot_chains and for each chain it will be saving its sequence, using aa dictionary.
    """
    for line in input_file:
        if line[0:4] == "ATOM":
            residue_name = line[17:20].replace(" ","")
            insertion = line[26]

            if insertion.isalpha() == False:
                if residue_name not in identifier:
                    if (line[21] == chain) and (line[22:26]!= residue_number):
                        try:
                            seq += aa[residue_name]
                        except:
                            seq += "X"
                        residue_number = line[22:26]
                    elif (line[21] != chain) and (seq != ""):
                        sequence.append(seq)
                        try:
                            seq = aa[residue_name]
                        except:
                            seq = "X"
                        residue_number= line[22:26]
                    elif (line[21]!= chain) and (seq == ""):
                        try:
                            seq = aa[residue_name]
                        except:
                            seq = "X"
                        residue_number= line[22:26]
                    if line[21] not in prot_chains:
                        prot_chains.append(line[21])
                        chain= prot_chains[-1]
                elif residue_name in identifier[:5]:
                    if line[21] not in dna_chains:
                        dna_chains.append(line[21])
                elif residue_name in identifier[5:]:
                    if line[21] not in rna_chains:
                        rna_chains.append(line[21])

                if "\n" in line:
                    csvfile.write(line)
                else:
                    csvfile.write(line+"\n")


        """
        if the line begins with HETATM, ligand can be a drug or other compounds not of interest
        if the ligand_nameis not a bad_drug, the line will be saved in list_of_lines and the 3 letter code of ligands will be appended to ligands
        type information of ligands in this pdb will be saved in ligand_types and pdbx_type code will be saved in the respective variable (nonpolymer, other, peptide, saccharide, dna or rna)
        when the program sees the type of a ligand, it increments in big_list within the correspondent index and its percentage in percentage list (NOTE: this numbers and percentage are made with every ligands! -good_drug)
        [18]:"sum_ligands_pdbs", [19]:"nonpolymer", [20]:"other", [21]:"peptide", [22]:"saccharide", [23]:"DNA", [24]:"RNA", [25]:"total_interest_ligands
        """

        if line[0:6] == "HETATM":
            ligand_name = line[17:20].strip()
            try:
                ligand_info = drugs.index(ligand_name)
                #print(ligand_info)
                ligand_info = comp_interest_id[ligand_info]

            except:
                continue

            if "\n" in line:
                csvfile.write(line)

            elif "\n" not in line:
                csvfile.write(line+"\n")

            if ligand_info in ligands_list:
                continue

            else:
                ligands_list.append(ligand_info)

        if line[0:5]== "MODEL":
            model_file.write(id+"\n")
            pdb_information = None
            return pdb_information


    if seq!="":
        sequence.append(seq)

    for tuple in ligands_list:
        ligand_types.append((tuple[1],tuple[2]))

        if tuple[1] < 26:
            if tuple[0] in ligands.keys():
                ligands[tuple[0]] += 1
            else:
                ligands[tuple[0]] = 1

        if tuple[1] == 1:
            nonpolymer_ligands.append((tuple[0],tuple[1]))

        elif tuple[1] == 2:
            other_ligands.append((tuple[0],tuple[1]))

        elif tuple[1] >= 3 and tuple[1] <= 14:
            peptide_ligands.append((tuple[0],tuple[1]))

        elif tuple[1] >= 15 and tuple[1] <= 25:
            saccharide_ligands.append((tuple[0],tuple[1]))

        if tuple[1] >= 26 and tuple[1] <= 28:
            if tuple[0] not in ligands_dna_rna.keys():
                dna_ligands.append((tuple[0],tuple[1]))

        elif tuple[1] >= 29 and tuple[1] <= 33:
            if tuple[0] not in ligands_dna_rna.keys():
                rna_ligands.append((tuple[0],tuple[1]))

        if tuple[1] >= 26:
            if tuple[0] in ligands_dna_rna.keys():
                ligands_dna_rna[tuple[0]] += 1
            else:
                ligands_dna_rna[tuple[0]] = 1

    pdb_information =[id, len(prot_chains), [p for p in prot_chains], [len(s) for s in sequence], len(rna_chains),
    [r for r in rna_chains], len(dna_chains), [d for d in dna_chains], len(ligands_list), [l[0] for l in ligands_list],
    [tuple[0] for tuple in ligand_types], [tuple[1] for tuple in ligand_types], str(len(prot_chains))+":"+str(len(ligands_list)),
    len(nonpolymer_ligands), [n for n in nonpolymer_ligands], len(other_ligands), [o for o in other_ligands],
    len(peptide_ligands), [pep for pep in peptide_ligands], len(saccharide_ligands), [s for s in saccharide_ligands],
    len(dna_ligands), [dna for dna in dna_ligands],len(rna_ligands), [rna for rna in rna_ligands]]


    """
    pdb_information saves all pdb information with pdb name as key: protein, rna and dna chains, ligands, ligand type and pdbx_type.
    Only for pdbs with protein, the last chain sequence will be saved in sequence and sequences will be compared to discover equal chains
    only PDBs with a protein will save list_of_lines in a PDB.
    This .csv has only the important lines of a PDB
    chain sequences and equal chains will be saved in sequence.csv
    """
    return pdb_information, sequence, ligands, ligands_dna_rna

"""
equal_chains function discovers all equal_chains in a protein.
is_equal begins with FALSE in every chains and the program will compare every chain except the last with all the others, only if both chains are FALSE.
When two chains are equal, the is_equal value changes to TRUE and the second chains is saved in other list
if the first chain is true, seqs_equal will append all equal chains and save them in equal_chains
sequence_list will append: pdb, sequence of protein chains and witch chains are equal
"""
def equal_chains_results(sequence, prot_chains):
    equal_chains = []
    num_chain=len(sequence)
    is_equal = [False]*num_chain
    if num_chain > 1:
        for c in range(num_chain-1):
            if not is_equal[c]:
                other =[]
                for k in range(c+1,num_chain):
                    if (sequence[c] == sequence[k]) and (not is_equal[k]):
                        other.append(prot_chains[k])
                        is_equal[k] = True
                        is_equal[c] = True
                if is_equal[c]:
                    seqs_equals = [prot_chains[c]]
                    for i in other:
                        seqs_equals.append(i)
                    equal_chains.append(seqs_equals)

    return equal_chains

