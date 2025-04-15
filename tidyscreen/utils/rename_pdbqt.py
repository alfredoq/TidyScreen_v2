

mol2_file = "/home/fredy/Desktop/test_tidyscreen/test3/chemspace/misc/yulia_lig_files/AAUCXXCMBFZYRT-UHFFFAOYSA-N_sybyl.mol2"
mol2_file_out = "/home/fredy/Desktop/test_tidyscreen/test3/chemspace/misc/yulia_lig_files/AAUCXXCMBFZYRT-UHFFFAOYSA-N_sybyl-RENAMED.mol2"
pdqbt_file = "/home/fredy/Desktop/test_tidyscreen/test3/chemspace/misc/yulia_lig_files/AAUCXXCMBFZYRT-UHFFFAOYSA-N_sybyl.pdbqt"

def get_atom_names_from_mol2(mol2_file):
    """Extracts atom names from a Mol2 file manually."""
    atom_names = []
    atom_section = False
    with open(mol2_file, "r") as file:
       
       for line in file:
            line_split = line.split()
            if len(line_split) > 0:
                if line_split[0] == "@<TRIPOS>ATOM":
                    atom_section = True
                    continue
                elif line_split[0] == "@<TRIPOS>BOND":
                    break  # End of ATOM section

            if  atom_section and len(line_split) > 1:
                atom_names.append(line_split[1])  # The second column is the atom name
    
    return atom_names

def rename_pdbqt_file(pdbqt_file,atom_names,output_file):
    with open(pdbqt_file, 'r') as file:
        id = 0
        for line in file:
            line_split = line.rstrip().split()
            if line_split[0] == "ATOM":
                print(f"replace {line_split[2]} by {atom_names[id]}")
                id += 1

atom_names = get_atom_names_from_mol2(mol2_file)
rename_pdbqt_file(pdqbt_file,atom_names,mol2_file_out)