import sqlite3
from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage, rdMolDraw2D
import sys
from pandarallel import pandarallel
import pandas as pd
import os
import shutil
from pathlib import Path
from rdkit.Chem import AllChem
import tarfile
import io
import subprocess
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

mol2_file = "/tmp/AAUCXXCMBFZYRT-UHFFFAOYSA-N_sybyl.mol2"
pdb_file = "/tmp/AAUCXXCMBFZYRT-UHFFFAOYSA-N.pdb"
target_pdqbt_file = "/tmp/AAUCXXCMBFZYRT-UHFFFAOYSA-N.pdbqt"

def get_atom_names_from_mol2(mol2_file):
    """Extracts atom names from a Mol2 file manually."""
    atom_names = []
    atom_ref_coord = []
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
                atom_ref_coord.append(round(float(line_split[2]),3)) # The third column is the 'x' coord of the atom name
  
    return atom_names, atom_ref_coord

def prepare_pdqbt_with_meeko(mol2_file):
    mol = Chem.MolFromMol2File(mol2_file,removeHs=False)
    mk_prep = MoleculePreparation()
    mol_setup_list = mk_prep(mol)
    molsetup = mol_setup_list[0]
    
    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)

    with open(target_pdqbt_file,'w') as pdbqt_file:
        pdbqt_file.write(pdbqt_string[0])

def rename_pdbqt_file(target_pdqbt_file,atom_names, atom_ref_coords):
    output_file = target_pdqbt_file.replace(".pdbqt", "-renamed.pdbqt")
    with open(target_pdqbt_file,'r') as readfile, open(output_file,'w') as writefile:
        for line in readfile:
            line_split = line.rstrip().split()
            if line_split[0] == "ATOM":
                if float(line_split[5]) in atom_ref_coords:
                    # Get the index of the corresponding coordinate
                    index = atom_ref_coords.index(float(line_split[5]))
                    # Replace the atom name
                    line_split[2] = atom_names[index]
                    # Parse a new line be written to the renamed .pdbqt file
                    column_widths = [9, 4, 4, 8, 6, 7, 8, 8, 6, 6, 10, 3]
                    # Format the line in accordance to .pdbqt 
                    formated_line = line_split[0].ljust(column_widths[0]) + line_split[1].ljust(column_widths[1]) + line_split[2].ljust(column_widths[2]) + line_split[3].ljust(column_widths[3]) + line_split[4].ljust(column_widths[4]) + line_split[5].rjust(column_widths[5]) + line_split[6].rjust(column_widths[6]) + line_split[7].rjust(column_widths[7]) + line_split[8].rjust(column_widths[8]) + line_split[9].rjust(column_widths[9]) + line_split[10].rjust(column_widths[10]) + " " + line_split[11].ljust(column_widths[11])
                    # Write the formatted line to the output file
                    writefile.write(f"{formated_line} \n")
            else:
                writefile.write(line)

    return writefile

atom_names, atom_ref_coords = get_atom_names_from_mol2(mol2_file)
prepare_pdqbt_with_meeko(mol2_file)
rename_pdbqt_file(target_pdqbt_file,atom_names,atom_ref_coords)




# mol = Chem.MolFromMol2File(mol2_file,removeHs=False)

# # Assign the extracted atom names to the molecule
# if mol:
#     for atom, name in zip(mol.GetAtoms(), atom_names):
#         atom.SetProp("atomName", name)  # Assign name as a property

# # # Show atom numbering
# # for atom in mol.GetAtoms():
# #     print(f"Atom Index: {atom.GetIdx()}, Name: {atom.GetProp('atomName')}")

# mk_prep = MoleculePreparation()
# mol_setup_list = mk_prep(mol)
# molsetup = mol_setup_list[0]

# pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)

# # Show atom numbering
# for atom in mol.GetAtoms():
#     print(f"Atom Index: {atom.GetIdx()}, Name: {atom.GetProp('atomName')}")

