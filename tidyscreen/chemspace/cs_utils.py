import warnings
warnings.filterwarnings("ignore")
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
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from tidyscreen.GeneralFunctions import general_functions as general_functions

def check_smiles(smiles):
    
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Problem reading SMILES columns - Example {smiles} \n Stopping...")
        sys.exit()
    else:
        print("SMILES column valid...")

def process_input_df(df,db,file):
    """
    Will rename the columns of the df constructed from a .csv input. If 'name' or 'flag' columns does not exist they will be created.
    """
    # Rename the columns of the df to give representative naming
    rename_dict = {
    0: 'SMILES',  
    1: 'name',  
    2: 'flag',   
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

    ### The processing of the SMILES will be done in the following order:
    # Initialize pandarallel to process all steps in parallel
    pandarallel.initialize(progress_bar=True) 
    # Sanitization performed in parallel
    print("Sanitizing SMILES")
    df_sanitized = pd.DataFrame()
    df_sanitized[["SMILES","name","flag"]] = df.parallel_apply(lambda row: sanitize_smiles_single(row,db,file), axis=1, result_type="expand")
    # Drop rows excluded by sanitization
    df_sanitized = df_sanitized.dropna()
    # Enumerate stereoisomers in parallel
    pandarallel.initialize(progress_bar=True) 
    print("Enumerating stereoisomers")
    df_enumerated = pd.DataFrame() # Create the enunmerated dataframe to return values from pandarallel
    df_enumerated[["SMILES","name","flag","stereo_nbr","stereo_config"]] = df_sanitized.parallel_apply(lambda row: enumerate_stereoisomers_single(row,db,file), axis=1, result_type="expand")
    # Computation the InChI key for the whole dataframe in parallel
    print("Computing InChIKey")
    pandarallel.initialize(progress_bar=True)
    df_enumerated["inchi_key"] = df_enumerated.parallel_apply(lambda row: compute_inchi_key_refactored(row,db,file),axis=1)
    # Delete duplicated molecules based on inchi_key
    df_ready_checked = df_enumerated.drop_duplicates(subset='inchi_key', keep='first')
    # Create an 'id' column
    df_ready_checked = df_ready_checked.reset_index().rename(columns={'index':'id'})

    return df_ready_checked

def sanitize_smiles(df,db,file):
    # Create the output dataframe containing the sanitized SMILES
    df_sanitized = pd.DataFrame(columns=['SMILES','name','flag'])
    for index, row in df.iterrows():
        try: 
            mol = Chem.MolFromSmiles(row["SMILES"])
            
            if mol is None:
                general_functions.write_failed_smiles_to_db(row["SMILES"],db,file)
                continue  # Invalid SMILES
            
            # Split into fragments and keep the largest
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if not frags:
                general_functions.write_failed_smiles_to_db(row["SMILES"],db,file)
                continue  # No fragments found
            
            # The largest fragment is the one with the most atoms and that successfully passed the sanitization
            largest = max(frags, key=lambda m: m.GetNumAtoms())
        
        except Exception as error:
            general_functions.write_failed_smiles_to_db(row["SMILES"],db,file)
            continue  # Error in processing SMILES
        
        # Sanitize the largest fragment
        try:
            Chem.SanitizeMol(largest)
            can_smi = Chem.MolToSmiles(largest, isomericSmiles=True)
            new_row = pd.DataFrame({'SMILES': can_smi, 'name': [row["name"]], 'flag': [row["flag"]]})
            df_sanitized = pd.concat([df_sanitized, new_row], ignore_index=True)
            
        except Exception:
            general_functions.write_failed_smiles_to_db(row["SMILES"],db,file)
            continue  # Error in sanitization
    
    # return the dataframe containing the sanitized SMILES
    return df_sanitized

### Developing a function for sanitizing the SMILES from a single input
def sanitize_smiles_single(row,db,file):
    smiles = row["SMILES"]
    name = row["name"]
    flag = row["flag"]
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as error:
        general_functions.write_failed_smiles_to_db(smiles,db,file)
        return
        
    if mol is None:
        general_functions.write_failed_smiles_to_db(smiles,db,file)
        return
            
    # Split into fragments and keep the largest
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not frags:
        general_functions.write_failed_smiles_to_db(smiles,db,file)
        return
            
    # The largest fragment is the one with the most atoms and that successfully passed the sanitization
    largest = max(frags, key=lambda m: m.GetNumAtoms())
        
    # Sanitize the largest fragment
    try:
        Chem.SanitizeMol(largest)
        can_smi = Chem.MolToSmiles(largest, isomericSmiles=True)
        return can_smi, name, flag
    
    except Exception:
        general_functions.write_failed_smiles_to_db(smiles,db,file)
        return

def enumerate_stereoisomers(df,db,file):
    options = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=32)
    df_enumerated = pd.DataFrame(columns=['id','SMILES','name','flag',"stereo_nbr","stereo_config"])
    for index, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["SMILES"])
        mol_hs = Chem.AddHs(mol)
        isomers = list(EnumerateStereoisomers(mol_hs,options=options))
        for isomer in isomers:
            isomer_no_hs = Chem.RemoveHs(isomer)
            smiles = Chem.MolToSmiles(isomer_no_hs)
            # Get the stereo configuration of the corresponding isomer
            stereo_config_tuple = Chem.FindMolChiralCenters(isomer,force=True,includeUnassigned=True,useLegacyImplementation=True)
            # Wil convert the list containing stereo into string for storage in SQL database
            stereo_config_list = [item[1] for item in stereo_config_tuple]
            stereo_config = ','.join(stereo_config_list)
            # Append the corresponding registry to the dataframe
            df_enumerated = df_enumerated._append({"SMILES":smiles,"name":row["name"],"flag":row["flag"],"stereo_nbr":len(isomers),"stereo_config":stereo_config},ignore_index=True)
            
    # Create the if column by assigning the index of the enumerated df
    df_enumerated['id'] = df_enumerated.index

    return df_enumerated

def enumerate_stereoisomers_single(row,db,file):
    smiles = row["SMILES"]
    name = row["name"]
    flag = row["flag"]
    # Set the options for the enumeration
    options = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=32)
    mol = Chem.MolFromSmiles(smiles)
    mol_hs = Chem.AddHs(mol)
    isomers = list(EnumerateStereoisomers(mol_hs,options=options))
    # Loop through the isomers and store the corresponding information
    for isomer in isomers:
        isomer_no_hs = Chem.RemoveHs(isomer)
        smiles = Chem.MolToSmiles(isomer_no_hs)
        # Get the stereo configuration of the corresponding isomer
        stereo_config_tuple = Chem.FindMolChiralCenters(isomer,force=True,includeUnassigned=True,useLegacyImplementation=True)
        # Wil convert the list containing stereo into string for storage in SQL database
        stereo_config_list = [item[1] for item in stereo_config_tuple]
        stereo_config = ','.join(stereo_config_list)
        if stereo_config == "":
            stereo_config = "n/a"
        
        return smiles, name, flag, len(isomers), stereo_config

def compute_inchi_key_for_whole_df(df,db,file):
    pandarallel.initialize(progress_bar=False)
    df["inchi_key"] = df["SMILES"].parallel_apply(lambda smiles: compute_inchi_key(smiles,db,file))
    return df

def compute_inchi_key_refactored(row,db,file):
    smiles = row["SMILES"]
    mol = Chem.MolFromSmiles(smiles)
    inchi_key = Chem.MolToInchiKey(mol)
    return inchi_key
    
def compute_inchi_key(smiles,db,file):
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchi_key = Chem.MolToInchiKey(mol)
        return inchi_key
        
    except:
        general_functions.write_failed_smiles_to_db(smiles,db,file)

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
    
    # Try to delete the corresponding table
    try:
        cursor.execute(f"DROP TABLE {table_name};")
        # Fetch all table names
        conn.commit()
        conn.close
        # Print table names
        print(f"Table '{table_name}' has been deleted from '{db}'")
    
    except Exception as error:
        print(f"Error deleting table '{table_name}': \n {error}")
        print("Stopping...")
        sys.exit()

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

def process_all_mols_in_table(db,table_name,charge,pdb,mol2_sybyl,mol2_gaff2,pdbqt):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    
    # # Define the list of column name to be created to store mol objects based on output selection
    # list_mol_objects_colnames = []
    # list_mol_objects_colnames_types = []

    # if pdb == 1:
    #     list_mol_objects_colnames.append("pdb_file")
    #     list_mol_objects_colnames_types.append("BLOB")
    # if mol2_sybyl == 1:
    #     list_mol_objects_colnames.append("mol2_file_sybyl")
    #     list_mol_objects_colnames_types.append("BLOB")
    # if mol2_gaff2 == 1:
    #     list_mol_objects_colnames.append("mol2_file_gaff")
    #     list_mol_objects_colnames_types.append("BLOB")
    #     list_mol_objects_colnames.append("frcmod_file")
    #     list_mol_objects_colnames_types.append("BLOB")
    # if pdbqt == 1:
    #     list_mol_objects_colnames.append("pdbqt_file")
    #     list_mol_objects_colnames_types.append("BLOB")

    # # Finally add the charge model if mol2 are requested:
    
    # if mol2_sybyl == 1 or mol2_gaff2 == 1:
    #     list_mol_objects_colnames.append("charge_model")
    #     list_mol_objects_colnames_types.append("TEXT")
        
    # # Exit if no molecules is selected
    # if len(list_mol_objects_colnames) == 0:
    #     print("No computation of molecules was requested. Stopping...")
    #     sys.exit()
    
    list_mol_objects_colnames, list_mol_objects_colnames_types = create_mols_columns_from_selection(pdb,mol2_sybyl,mol2_gaff2,pdbqt)
        
    check_columns_existence_in_table(conn, table_name, list_mol_objects_colnames)

    # Create all the columns to store mol objects
    try:
        for index, column in enumerate(list_mol_objects_colnames):
            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column} {list_mol_objects_colnames_types[index]};")
    except Exception as error:
        print(error)
        print("Stopping")
        sys.exit()


    # Compute the corresponding molecules
    sql = f"""SELECT id, SMILES, Name FROM '{table_name}';"""
    df = pd.read_sql_query(sql,conn)
    pandarallel.initialize(progress_bar=True) # Activate Progress Bar
    df.parallel_apply(lambda row: append_ligand_mols_blob_object_to_table(db,table_name,row,charge,pdb,mol2_sybyl,mol2_gaff2,pdbqt), axis=1)

    try:
        pass
        #clean_temp_dir(db,table_name)
    except:
        pass

def create_mols_columns_from_selection(pdb,mol2_sybyl,mol2_gaff2,pdbqt):
    # Define the list of column name to be created to store mol objects based on output selection
    list_mol_objects_colnames = []
    list_mol_objects_colnames_types = []

    if pdb == 1:
        list_mol_objects_colnames.append("pdb_file")
        list_mol_objects_colnames_types.append("BLOB")
    if mol2_sybyl == 1:
        list_mol_objects_colnames.append("mol2_file_sybyl")
        list_mol_objects_colnames_types.append("BLOB")
    if mol2_gaff2 == 1:
        list_mol_objects_colnames.append("mol2_file_gaff")
        list_mol_objects_colnames_types.append("BLOB")
        list_mol_objects_colnames.append("frcmod_file")
        list_mol_objects_colnames_types.append("BLOB")
    if pdbqt == 1:
        list_mol_objects_colnames.append("pdbqt_file")
        list_mol_objects_colnames_types.append("BLOB")

    # Finally add the charge model if mol2 are requested:
    
    if mol2_sybyl == 1 or mol2_gaff2 == 1:
        list_mol_objects_colnames.append("charge_model")
        list_mol_objects_colnames_types.append("TEXT")
        
    # Exit if no molecules is selected
    if len(list_mol_objects_colnames) == 0:
        print("No computation of molecules was requested. Stopping...")
        sys.exit()


    return list_mol_objects_colnames, list_mol_objects_colnames_types

def check_columns_existence_in_table(conn,table_name,columns_list):
    cursor = conn.cursor()
    try:
        cursor.execute(f"PRAGMA table_info({table_name})")
    except Exception as error:
        print(error)
        
    columns = [row[1] for row in cursor.fetchall()]
    # Return True if all items in 'columns_list' already exist as columns in 'table_name'
    if all(item in columns for item in columns_list):
        print(f"The columns: \n '{columns_list} \n already exists in '{table_name}'. \n To process again delete the whole table and recompute molecules. Stoppping...")

        sys.exit()

def append_ligand_mols_blob_object_to_table(db,table_name,row,charge,pdb,mol2_sybyl,mol2_gaff2,pdbqt):
    
    action = 0 
    
    if pdb == 1:
        # Prepare and store the corresponding .pdb file
        pdb_file, tar_pdb_file, net_charge = pdb_from_smiles(row["SMILES"])
        store_file_as_blob(db,table_name,'pdb_file',tar_pdb_file,row)
        action = 1
    
    if mol2_sybyl == 1 and pdb == 1:
        # Prepare and store the corresponding .mol2 file
        output_file_sybyl, output_tar_file_sybyl = mol2_from_pdb(pdb_file,net_charge,charge,at="sybyl")
        # Store the .mol2 file - atom type: Sybyl
        store_file_as_blob(db,table_name,'mol2_file_sybyl',output_tar_file_sybyl,row)
    
    if mol2_gaff2 == 1 and pdb == 1:
        # Prepare and store the corresponding .mol2 file
        output_file_gaff, output_tar_file_gaff = mol2_from_pdb(pdb_file,net_charge,charge,at="gaff")
        # Store the .mol2 file - atom type: gaff
        store_file_as_blob(db,table_name,'mol2_file_gaff',output_tar_file_gaff,row)
        # Compute the .frcmod file
        output_tar_frcmod = compute_frcmod_file(output_file_gaff,at="gaff")
        # Store the .frcmod file
        store_file_as_blob(db,table_name,'frcmod_file',output_tar_frcmod,row)
    
    if pdbqt == 1 and mol2_sybyl == 1 and pdb == 1:
        #Prepare and store the corresponding .pdbqt file
        tar_pdbqt_file = pdbqt_from_mol2(output_file_sybyl)
        #print(tar_pdbqt_file)
        # # Store the .pdbqt file
        store_file_as_blob(db,table_name,'pdbqt_file',tar_pdbqt_file,row)

    if mol2_sybyl == 1 or mol2_gaff2 == 1 and pdb == 1:
        # Append the charge model to the row
        store_string_in_column(db,table_name,"charge_model",charge,row)

    if action == 0:
        print("No computation of molecules was requested. Stopping...")
        sys.exit()

def pdb_from_smiles(smiles):
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol_hs = Chem.AddHs(mol)
        confs, mol_hs, ps = generate_ligand_conformers(mol_hs)
        selected_mol = get_conformer_rank(confs, mol_hs, ps)
        pdb_file, inchi_key = save_ligand_pdb_file(selected_mol)
        # Compute the net charge of the molecule for potential use in antechamber
        net_charge = compute_molecule_net_charge(mol)
            
        # Compress the pdb_file
        tar_pdb_file = generate_tar_file(pdb_file)
        
        return pdb_file, tar_pdb_file, net_charge
    
    except Exception as error:
        print(f"Problem with SMILES: \n {smiles}")

def mol2_from_pdb(pdb_file,net_charge,charge,at):
    # Get the prefixes for the file
    file_prefix = pdb_file.split('/')[-1].replace(".pdb","")
    # Compute .mol2 file
    antechamber_path = shutil.which('antechamber')
    
    temp_folder = f"/tmp/{file_prefix}"
    os.makedirs(temp_folder, exist_ok=True)
    
    output_file = f"{temp_folder}/{file_prefix}_{at}.mol2"
    # Compute mol2 with Sybyl Atom Types - for compatibility with RDKit and Meeko
    antechamber_command = f'cd {temp_folder} && {antechamber_path} -i {pdb_file} -fi pdb -o {output_file} -fo mol2 -c {charge} -nc {net_charge} -at {at} -pf y' # The 'sybyl' atom type 
    #convention is used for compatibility with RDKit
    subprocess.run(antechamber_command, shell=True, capture_output=True, text=True)
        
    # Compress the output file for storage
    output_tar_file = generate_tar_file(output_file)
        
    return output_file, output_tar_file

def compute_frcmod_file(mol2_file,at):
    # Get the prefixes for the file
    file_prefix = mol2_file.split('/')[-1].replace(f"_{at}.mol2","")
    temp_folder = f"/tmp/{file_prefix}"
    output_file = f"{temp_folder}/{file_prefix}.frcmod"
    # Compute the frcmod file
    os.makedirs(temp_folder, exist_ok=True)
    parmchk2_path = shutil.which('parmchk2')
    
    parmchk2_command = f'{parmchk2_path} -i {mol2_file} -f mol2 -o {output_file}' 
    subprocess.run(parmchk2_command, shell=True, capture_output=True, text=True)
    
    # Compress the output file for storage
    output_tar_file = generate_tar_file(output_file)
    
    return output_tar_file
    
def pdbqt_from_mol2(mol2_file):
    # Load the corresponding .mol2 file 
    file_prefix = mol2_file.split('/')[-1].replace('_sybyl.mol2','')
    
    try:
        mol = Chem.MolFromMol2File(mol2_file,removeHs=False)
    except:
        print(mol2_file)
    
    atoms_dict = create_meeko_atoms_dict()
    mk_prep = MoleculePreparation(merge_these_atom_types=("H"),charge_model="read", charge_atom_prop="_TriposPartialCharge",add_atom_types=atoms_dict)
    
    mol_setup_list = mk_prep(mol)
    molsetup = mol_setup_list[0]

    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)
    
    pdbqt_outfile = f'/tmp/{file_prefix}/{file_prefix}_tmp.pdbqt'
    with open(pdbqt_outfile,'w') as pdbqt_file:
        pdbqt_file.write(pdbqt_string[0])

    # In this section, .pdbqt atoms will be renamed to match the original .mol2 file

    atom_names, atom_ref_coords = get_atom_names_from_mol2(mol2_file)
    renamed_pdbqt_file = rename_pdbqt_file(pdbqt_outfile,atom_names, atom_ref_coords,file_prefix)

    tar_pdbqt_file = generate_tar_file(renamed_pdbqt_file)
    
    return tar_pdbqt_file

def store_file_as_blob(db,table_name,colname,tar_pdb_file,row):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    id = row['id']
    # Will store the '.pdb' file as a blob object into the database
    cursor.execute(f"UPDATE {table_name} SET {colname} = ? WHERE id = {id};",(tar_pdb_file,))
    conn.commit()

def store_string_in_column(db,table_name,colname,string,row):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    id = row['id']
    # Will store the '.pdb' file as a blob object into the database
    cursor.execute(f"UPDATE {table_name} SET {colname} = ? WHERE id = {id};",(string,))
    conn.commit()

def compute_molecule_net_charge(mol):
    molecule_charge = 0
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        molecule_charge = molecule_charge + charge

    return molecule_charge

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

def rename_pdbqt_file(target_pdqbt_file,atom_names, atom_ref_coords,file_prefix):
    output_file = f"/tmp/{file_prefix}/{file_prefix}.pdbqt"
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
    filename = file.split("/")[-1]
    tar_buffer = io.BytesIO()  # Create an in-memory buffer
    with tarfile.open(fileobj=tar_buffer, mode="w") as tar:
        tar.add(file, arcname=filename)  # Store only filename in TAR
    
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
    temp_folder = f'/tmp/{inchi_key}'
    os.makedirs(temp_folder, exist_ok=True)
    pdb_file = f'{temp_folder}/{inchi_key}.pdb'
    
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
            with tarfile.open(output_file, 'r') as tar_ref:
                try:
                    tar_ref.extractall(path=output_path)
                except Exception as error:
                    print(error)
                
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

def clean_temp_dir(db,table_name):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    
    # Get the inchi key values to construct a list of directories to delete
    cursor.execute(f"SELECT inchi_key FROM {table_name}")
    rows = cursor.fetchall()
    
    # Extract the values into a list
    inchikey_values = [row[0] for row in rows]
    
    for inchi_key in inchikey_values:
        shutil.rmtree(f"/tmp/{inchi_key}")
    
    conn.close
    
    # extensions = ["mol2", "pdb", "pdbqt", "frcmod"]  # Add extensions to delete
    # for ext in extensions:
    #     pattern = os.path.join(directory, f"*.{ext}")
    #     for file in glob.glob(pattern):
    #         try:
    #             os.remove(file)
    #         except Exception as e:
    #             pass