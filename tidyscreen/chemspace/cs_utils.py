import sqlite3
from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage, rdMolDraw2D
import sys
from tidyscreen import tidyscreen as tidyscreen
from pandarallel import pandarallel
import pandas as pd
import os
import glob
import shutil
from pathlib import Path
from rdkit.Chem import AllChem
import tarfile
import io
import subprocess
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

def check_smiles(smiles):
    
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Problem reading SMILES columns - Example {smiles} \n Stopping...")
        sys.exit()
    else:
        print("SMILES column valid...")

def process_input_df(df):
    """
    Will rename the columns of the df constructed from a .csv input. If 'name' or 'flag' columns does not exist they will be created.
    """
    rename_dict = {
    'index': 'id', # Renaming the index column
    0: 'SMILES',  # Renaming 'A' to 'New_A'
    1: 'name',  # Renaming 'B' to 'New_B'
    2: 'flag'   # 'C' doesn't exist, it will be created
    }
    
    for old_col, new_col in rename_dict.items():
        if old_col in df.columns:
            # Rename existing columns
            df = df.rename(columns={old_col: new_col})
        else:
            # Create new columns if they don't exist
            if new_col == "name":
                df[new_col] = "nn" 
            elif new_col == "flag":
                df[new_col] = 0

    df = compute_inchi_key_for_whole_df(df)

    return df

def compute_inchi_key_for_whole_df(df):
    pandarallel.initialize(progress_bar=False)
    df["inchi_key"] = df["SMILES"].parallel_apply(lambda smiles: compute_inchi_key(smiles))
    return df

def compute_inchi_key(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchi_key = Chem.MolToInchiKey(mol)
        return inchi_key
    except:
        pass

def check_table_presence(conn,table_name):
    "Will return 1 if exists, otherwise returns 0"
    cursor = conn.cursor()
    cursor.execute("""SELECT name FROM sqlite_master WHERE type='table' AND name=?;""", (table_name,))

    table_exists = cursor.fetchone()

    if table_exists:
        print(f"Table {table_name} exists.")
        return 1
    else:
        #print(f"Table '{table_name}' DOES NOT exists.")
        return 0

def save_df_to_db(db,df,table_name):
    conn = tidyscreen.connect_to_db(db)
    exists = check_table_presence(conn,table_name)

    if exists == 1:
        replace_table_action = input(f"The table named '{table_name}' already exists in the database. I will replace it, are you ok with that? (y/n): ")

        if replace_table_action == 'y':
            print(f"I will overwrite table '{table_name}' as indicated.'")
            pass
        else:
            print(f"Quiting to safely retain table {table_name}. Stopping...")
            sys.exit()

    try:
        df.to_sql(con=conn, name=table_name,if_exists="replace",index=None)
        
    except Exception as error:
        print(error)

def list_ligands_tables(db):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    # Query to get all table names
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    # Fetch all table names
    tables = cursor.fetchall()
    # Print table names
    print("Tables in the database:")
    for table in tables:
        print(table[0])

    conn.close()

def delete_ligands_table(db,table_name):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    # Query to get all table names
    cursor.execute(f"DROP TABLE IF EXISTS {table_name};")
    # Fetch all table names
    conn.commit()
    conn.close
    # Print table names
    print(f"Table '{table_name}' has been deleted from '{db}'")

def depict_ligands_table(db,table_name,output_path,max_mols_ppage=25):
    conn = tidyscreen.connect_to_db(db)
    sql=f"SELECT SMILES, inchi_key, name FROM {table_name};"
    molecules_df = pd.read_sql_query(sql,conn)

    # Since the molecules might be processed in chunks, it is usefull to create lists with the information in order to process them
    smiles_list, inchi_key_list, names_list = molecules_df["SMILES"].to_list(), molecules_df["inchi_key"].to_list(),molecules_df["name"].to_list()

    # Generate a naming list
    lengeds_list = combine_strings_in_lists([inchi_key_list,names_list])

    # Create the chunks of objects to process according to the max_mols_ppage parameter provided
    smiles_list_chunks = split_list_in_chunks(smiles_list, max_mols_ppage)
    lengeds_list_chunks = split_list_in_chunks(lengeds_list, max_mols_ppage)

    counter = 0 # used to number figures with chunks
    for index, list_item in enumerate(smiles_list_chunks):
        df_chunk  = pd.DataFrame({'SMILES':list_item})
        mols_item = [Chem.MolFromSmiles(smile) for smile in list_item]
        img = MolsToGridImage(mols=mols_item, legends=lengeds_list_chunks[index], molsPerRow=5,subImgSize=(400,250))
        img.save(f"{output_path}/{table_name}_{counter}.png")
        counter+=1

def combine_strings_in_lists(list_of_lists):
    # Take the first list as reference
    reference_list = list_of_lists[0]

    # Check if some list in the list of list is not equal in length. In that case, inform and exit
    for item in list_of_lists:
        if len(reference_list) != len(item):
            print(f"The list: \n '{item}' \n is not equal in length to the rest of lists provided. Stopping.")
            sys.exit()

    # Append the string values to the resulting list
    appended_list = [' \n '.join(strings) for strings in zip(*list_of_lists)]

    return appended_list

def split_list_in_chunks(list,nbr_items):
    """"
    Will receive a list of items, and will return a list of lists each one containing the corresponding 'nbr_items'
    ------
    Parameters
    ------
    -list: input list to be processed in chunks
    -nbr_item: number of chunks included in each component list
    """
    list_of_chunks = [list[i:i+nbr_items] for i in range(0, len(list), nbr_items)]

    return list_of_chunks

def check_folder_presence(folder,create=0):
    if os.path.isdir(folder):
        exists = input(f"Folder {folder} already exists, delete and continue? (y/n): ")

        if exists == 'y':
            print(folder)
            shutil.rmtree(folder)
            print(f"Folder: '{folder}' deleted. Continuing...")
            if create == 1:
                # Create the corresponding folder
                Path(f"{folder}").mkdir(parents=True, exist_ok=False)
        else:
            print(f"Preventing overwritting of '{folder}'. Stopping...")
            sys.exit()
    else:
        if create == 1:
            # Create the corresponding folder
            Path(f"{folder}").mkdir(parents=True, exist_ok=False)

def process_all_mols_in_table(db,table_name):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    try: 
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN pdb_file BLOB;") # Will create a 'pdb_file' column of type BLOB
        conn.commit()
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN mol2_file_sybyl BLOB;") # Will create a 'mol2_file' column of type BLOB
        conn.commit()
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN mol2_file_gaff BLOB;") # Will create a 'mol2_file' column of type BLOB
        conn.commit()
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN frcmod_file BLOB;") # Will create a 'frcmod_file' column of type BLOB
        conn.commit()
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN pdbqt_file BLOB;") # Will create a 'pdbqt_file' column of type BLOB
        conn.commit()
    except Exception as error:
        print(f"Error creating colums in: '{table_name}'")
        print(error)
        sys.exit()
    
    sql = f"""SELECT id, SMILES, Name FROM '{table_name}';"""
    df = pd.read_sql_query(sql,conn)
    pandarallel.initialize(progress_bar=True)
    df.parallel_apply(lambda row: append_ligand_mols_blob_object_to_table(db,table_name,row), axis=1)

    # Once all files have been generated and stored in the db, clean '/tmp' dir
    clean_dir("/tmp")

def clean_dir(directory):
    extensions = ["mol2", "pdb", "pdbqt", "frcmod"]  # Add extensions to delete
    for ext in extensions:
        pattern = os.path.join(directory, f"*.{ext}")
        for file in glob.glob(pattern):
            try:
                os.remove(file)
            except Exception as e:
                pass

def append_ligand_mols_blob_object_to_table(db,table_name,row):
    # Prepare and store the corresponding .pdb file
    pdb_file, tar_pdb_file, net_charge = pdb_from_smiles(row["SMILES"])
    store_file_as_blob(db,table_name,'pdb_file',tar_pdb_file,row)
    
    # Prepare and store the corresponding .mol2 file
    mol2_sybyl_file, tar_mol2_sybyl_file, tar_mol2_gaff_file, tar_frcmod_file = mol2_from_pdb(pdb_file,net_charge)
    # Store the .mol2 file - atom type: Sybyl
    store_file_as_blob(db,table_name,'mol2_file_sybyl',tar_mol2_sybyl_file,row)
    
    # Store the .mol2 file - atom type: gaff
    store_file_as_blob(db,table_name,'mol2_file_gaff',tar_mol2_gaff_file,row)
    # Store the .frcmod file
    store_file_as_blob(db,table_name,'frcmod_file',tar_frcmod_file,row)
    
    #Prepare and store the corresponding .pdbqt file
    tar_pdbqt_file = pdbqt_from_mol2(mol2_sybyl_file)
    # Store the .pdbqt file
    store_file_as_blob(db,table_name,'pdbqt_file',tar_pdbqt_file,row)

def store_file_as_blob(db,table_name,colname,tar_pdb_file,row):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    id = row['id']
    # Will store the '.pdb' file as a blob object into the database
    cursor.execute(f"UPDATE {table_name} SET {colname} = ? WHERE id = {id};",(tar_pdb_file,))
    conn.commit()

def pdb_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol_hs = Chem.AddHs(mol)
    confs, mol_hs, ps = generate_ligand_conformers(mol_hs)
    selected_mol = get_conformer_rank(confs, mol_hs, ps)
    pdb_file, inchi_key = save_ligand_pdb_file(selected_mol)

    net_charge = compute_molecule_net_charge(mol)
        
    # Compress the pdb_file
    tar_pdb_file = generate_tar_file(pdb_file)
    
    return pdb_file, tar_pdb_file, net_charge

def compute_molecule_net_charge(mol):
    molecule_charge = 0
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        molecule_charge = molecule_charge + charge

    return molecule_charge

def mol2_from_pdb(pdb_file,net_charge):
    # Get the prefixes for the file
    file_prefix = pdb_file.split('/')[-1].replace(".pdb","")
    # Compute .mol2 file
    antechamber_path = shutil.which('antechamber')
    
    # Compute mol2 with Sybyl Atom Types - for compatibility with RDKit and Meeko
    command1 = f'{antechamber_path} -i {pdb_file} -fi pdb -o /tmp/{file_prefix}_sybyl.mol2 -fo mol2 -c bcc -nc {net_charge} -at sybyl -pf y' # The 'sybyl' atom type convention is used for compatibility with RDKit
    subprocess.run(command1, shell=True, capture_output=True, text=True)
    
    # Compress the sybyl.mol2 file for storage
    tar_mol2_sybyl_file = generate_tar_file(f'/tmp/{file_prefix}_sybyl.mol2')
    
    # Compute mol2 with gaff Atom Types - for compatibility with Amber
    command2 = f'{antechamber_path} -i {pdb_file} -fi pdb -o /tmp/{file_prefix}_gaff.mol2 -fo mol2 -c bcc -nc {net_charge} -at gaff -pf y' # The 'sybyl' atom type convention is used for compatibility with RDKit
    subprocess.run(command2, shell=True, capture_output=True, text=True)
    
    # Compress the gaff.mol2 file for storage
    tar_mol2_gaff_file = generate_tar_file(f'/tmp/{file_prefix}_gaff.mol2')
    
    # Compute .frcmod file consistent with the mol2
    parmchk_path = shutil.which('parmchk2')
    command2 = f'{parmchk_path} -i /tmp/{file_prefix}_gaff.mol2 -f mol2 -o /tmp/{file_prefix}.frcmod'
    subprocess.run(command2, shell=True, capture_output=True, text=True)
    
    # Compress the .frcmod file for storage
    tar_frcmod_file = generate_tar_file(f'/tmp/{file_prefix}.frcmod')

    return f'/tmp/{file_prefix}_sybyl.mol2', tar_mol2_sybyl_file, tar_mol2_gaff_file, tar_frcmod_file

def pdbqt_from_mol2(mol2_file):
    # Load the corresponding .mol2 file 
    file_prefix = mol2_file.split('/')[-1].replace('.mol2','')
    mol = Chem.MolFromMol2File(mol2_file,removeHs=False)
    
    atoms_dict = create_meeko_atoms_dict()
    mk_prep = MoleculePreparation(merge_these_atom_types=("H"),charge_model="read", charge_atom_prop="_TriposPartialCharge",add_atom_types=atoms_dict)
    
    mol_setup_list = mk_prep(mol)
    molsetup = mol_setup_list[0]

    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)
    
    with open(f'/tmp/{file_prefix}.pdbqt','w') as pdbqt_file:
        pdbqt_file.write(pdbqt_string[0])

    # In this section, .pdbqt atoms will be renamed to match the original .mol2 file

    atom_names, atom_ref_coords = get_atom_names_from_mol2(mol2_file)
    renamed_pdbqt_file = rename_pdbqt_file(f'/tmp/{file_prefix}.pdbqt',atom_names, atom_ref_coords)

    tar_pdbqt_file = generate_tar_file(renamed_pdbqt_file)
    
    return tar_pdbqt_file

def get_atom_names_from_mol2_REDESIGNED(mol2_file):
    """Extracts atom names from a Mol2 file manually.
    This function has been redesigned in order to extraxt the (x,y,z) coordinates with an atom name, so that the renaming is based on the three coordinates
    """
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
                x_coord = (round(float(line_split[2]),3))
                y_coord = (round(float(line_split[3]),3))
                z_coord = (round(float(line_split[4]),3))
                atom_ref_coord.append((x_coord,y_coord,z_coord)) # Construct a tuple containing the x,y,z coordinates
  
    return atom_names, atom_ref_coord

def rename_pdbqt_file(target_pdqbt_file,atom_names, atom_ref_coords):
    output_file = target_pdqbt_file.split('_')[0] + ".pdbqt"
    with open(target_pdqbt_file,'r') as readfile, open(output_file,'w') as writefile:
        for line in readfile:
            line_split = line.rstrip().split()
            if line_split[0] == "ATOM":
                # Construct a tuple containing (x,y,z) coordinates
                x_coord = (round(float(line_split[5]),3))
                y_coord = (round(float(line_split[6]),3))
                z_coord = (round(float(line_split[7]),3))
                coords_tuple = (x_coord,y_coord,z_coord)
                if coords_tuple in atom_ref_coords:
                    # Get the index of the corresponding coordinates tuple
                    index = atom_ref_coords.index(coords_tuple)
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
    
    return output_file

def create_meeko_atoms_dict():
    atoms_dict = [{"smarts": "[#1]", "atype": "H",},
    {"smarts": "[#1][#7,#8,#9,#15,#16]","atype": "HD"},
    {"smarts": "[#5]", "atype": "B"},
    {"smarts": "[C]", "atype": "C"},
    {"smarts": "[c]", "atype": "A"},
    {"smarts": "[#7]", "atype": "NA"},
    {"smarts": "[#8]", "atype": "OA"},
    {"smarts": "[#9]", "atype": "F"},
    {"smarts": "[#12]", "atype": "Mg"},
    {"smarts": "[#14]", "atype": "Si"},
    {"smarts": "[#15]", "atype": "P"},
    {"smarts": "[#16]", "atype": "S"},
    {"smarts": "[#17]", "atype": "Cl"},
    {"smarts": "[#20]", "atype": "Ca"},
    {"smarts": "[#25]", "atype": "Mn"},
    {"smarts": "[#26]", "atype": "Fe"},
    {"smarts": "[#30]", "atype": "Zn"},
    {"smarts": "[#35]", "atype": "Br"},
    {"smarts": "[#53]", "atype": "I"},
    {"smarts": "[#7X3v3][a]", "atype": "N", "comment": "pyrrole, aniline"},
    {"smarts": "[#7X3v3][#6X3v4]", "atype": "N", "comment": "amide"},
    {"smarts": "[#7+1]", "atype": "N", "comment": "ammonium, pyridinium"},
    {"smarts": "[SX2]", "atype": "SA", "comment": "sulfur acceptor"},
    {"smarts": "[#1][#6X3]:[#6X3]([#6X4])[#7]:[#7][#7][#6X4]", "atype": "HD", "comment": "4,5-H in 1,2,3-triazole"},
    ]

    return atoms_dict

def generate_tar_file(file):
    """Compress the pdb file into a TAR archive and return its binary content."""
    tar_buffer = io.BytesIO()  # Create an in-memory buffer
    with tarfile.open(fileobj=tar_buffer, mode="w") as tar:
        tar.add(file, arcname=file.split("/")[-1])  # Store only filename in TAR
    
    return tar_buffer.getvalue()  # Return TAR file as binary

def generate_ligand_conformers(mol,nbr_confs=50,mmff='MMFF94',maxIters=10, conformer_rank=0):
    props = AllChem.ETKDG()
    props.pruneRmsThresh = 0.25
    props.useRandomCoords = True
    props.numThreads = 1
    
    # Create the object containing the conformers
    confs = AllChem.EmbedMultipleConfs(mol,nbr_confs,props)
    ps = AllChem.MMFFGetMoleculeProperties(mol)
    #AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant=mmff, maxIters=maxIters)
    AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant='MMFF94s', maxIters=maxIters)
    
    return confs, mol, ps

def get_conformer_rank(confs, mol_hs, ps,rank=0):
    
    conformers_energies_dict={} # Empty dictionary to store conformers

    for conf in confs:
        try:
            ff = AllChem.MMFFGetMoleculeForceField(mol_hs,ps,confId=conf)
            ff.Minimize()
            energy_value = ff.CalcEnergy()
            conformers_energies_dict[conf] = energy_value
        except:
            continue

    # Sort the dictionary 
    conformers_energies_dict_sorted=sorted(conformers_energies_dict.items(), key=lambda x: x[1])

    # The following will store the lowest energy conformer as .pdbqt
    selected_conformer = Chem.MolToMolBlock(mol_hs,confId=conformers_energies_dict_sorted[rank][0])
    selected_mol = Chem.MolFromMolBlock(selected_conformer, removeHs=False)

    return selected_mol

def save_ligand_pdb_file(selected_mol):
    inchi_key = Chem.MolToInchiKey(selected_mol)
    pdb_file = f'/tmp/{inchi_key}.pdb'
    
    Chem.MolToPDBFile(selected_mol,pdb_file)

    clean_pdb_records(pdb_file)
    
    return pdb_file, inchi_key

def clean_pdb_records(pdb_file):
    with open(pdb_file, 'r') as file:
        lines_to_check = file.readlines()
    
    with open(pdb_file,'w') as file:
        for line in lines_to_check:
            if 'HETATM' in line: # Will mantain only the lines containing ligand atom
                file.write(line)
                
def create_blob_object(file):
    """
    This function will create a binary object from a file provided that is aimed to be stored as a B
LOB object within the database.
    ------
    Parameters:
    ------
    - filename: the filename in the corresponding format to be stored in the table containing BLOB objec
ts.
    ------
    Returns:
    ------
    A blob object containing the ligand .pdbqt file in a compressed (.tar.gz) format.
    """
    with open(file, 'rb') as file:
        blobData = file.read()
    
    file.close()
    
    return blobData

def retrieve_blob_ligfiles(db,table_name,output_path,ligname,blob_colname):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    
    if ligname == None:

        # Retrieve all BLOBs
        cursor.execute(f"SELECT inchi_key, {blob_colname} FROM {table_name}")  # Adjust table/column names
        rows = cursor.fetchall()
        
        for row in rows:
            inchi_key, file_blob = row
            output_file = os.path.join(output_path, f"{inchi_key}.tar")
            
            # Write raw BLOB data in compressed mode
            with open(output_file, "wb") as f:
                f.write(file_blob)
            
            # Extract the file
            with tarfile.open(output_file, 'r:*') as tar_ref:
                tar_ref.extractall(output_path)
                os.remove(output_file)


    else:
        # Retrieve custom ligand by 'ligname'
        
        try: 
    
            cursor.execute(f"SELECT inchi_key, {blob_colname} FROM {table_name}")  # Adjust table/column names
            rows = cursor.fetchall()
            
            for row in rows:
                inchi_key, file_blob = row
                output_file = os.path.join(output_path, f"{inchi_key}.tar")
                
                if inchi_key == ligname:
                    # Write raw BLOB data in compressed mode
                    with open(output_file, "wb") as f:
                        f.write(file_blob)
                    
                    # Extract the file
                    with tarfile.open(output_file, 'r:*') as tar_ref:
                        tar_ref.extractall(output_path)
                        os.remove(output_file)
                        
        except Exception as error:
            print("Error retrieving ligand. Stopping...")
            print(error)
            sys.exit()


