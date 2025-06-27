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
import time
import pickle
import multiprocessing
from espaloma_charge import charge
from rdkit.Chem import Descriptors
import json


def check_smiles(smiles):
    
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Problem reading SMILES columns - Example {smiles} \n Stopping...")
        sys.exit()
    else:
        print("SMILES column valid...")

def process_input_df(df,db,file,stereo_enum):
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
    df[["SMILES","name","flag"]] = df.parallel_apply(lambda row: sanitize_smiles_single(row,db,file), axis=1, result_type="expand")
    # Drop rows excluded by sanitization
    df = df.dropna()
    # Enumerate stereoisomers in parallel
    pandarallel.initialize(progress_bar=True) 
    # Enumerate stereoisomers if requested
    if stereo_enum == 1:
        print("Enumerating stereoisomers")
        df = pd.DataFrame() # Create the enumerated dataframe to return values from pandarallel
        df[["SMILES","name","flag","stereo_nbr","stereo_config"]] = df.parallel_apply(lambda row: enumerate_stereoisomers_single(row,db,file), axis=1, result_type="expand")
        # Computation the InChI key for the whole dataframe in parallel
    
    # Compute the InChIKey for the whole dataframe
    print("Computing InChIKey")
    pandarallel.initialize(progress_bar=True)
    df["inchi_key"] = df.parallel_apply(lambda row: compute_inchi_key_refactored(row,db),axis=1)
    # Delete duplicated molecules based on inchi_key
    df = df.drop_duplicates(subset='inchi_key', keep='first')
    # Create an 'id' column
    df = df.reset_index().rename(columns={'index':'id'})

    return df

def sanitize_smiles(df,db,file):
    # Create the output dataframe containing the sanitized SMILES
    df_sanitized = pd.DataFrame(columns=['SMILES','name','flag'])
    for index, row in df.iterrows():
        try: 
            mol = Chem.MolFromSmiles(row["SMILES"])
            
            if mol is None:
                fail_message = "Failed at .csv input stage - mol read"
                general_functions.write_failed_smiles_to_db(row["SMILES"],db,file,fail_message)
                #continue  # Invalid SMILES
            
            # Split into fragments and keep the largest
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if not frags:
                fail_message = "Failed at .csv input stage - fragments generation"
                general_functions.write_failed_smiles_to_db(row["SMILES"],db,file,fail_message)
                #continue  # No fragments found
            
            # The largest fragment is the one with the most atoms and that successfully passed the sanitization
            largest = max(frags, key=lambda m: m.GetNumAtoms())
        
        except Exception as error:
            fail_message = "Failed at .csv input stage - general error"
            general_functions.write_failed_smiles_to_db(row["SMILES"],db,file,fail_message)
            #continue  # Error in processing SMILES
        
        # Sanitize the largest fragment
        try:
            Chem.SanitizeMol(largest)
            can_smi = Chem.MolToSmiles(largest, isomericSmiles=True)
            new_row = pd.DataFrame({'SMILES': can_smi, 'name': [row["name"]], 'flag': [row["flag"]]})
            df_sanitized = pd.concat([df_sanitized, new_row], ignore_index=True)
            
        except Exception:
            fail_message = "Failed at .csv input stage - sanitization error"
            general_functions.write_failed_smiles_to_db(row["SMILES"],db,file,fail_message)
            continue  # Error in sanitization
    
    # return the dataframe containing the sanitized SMILES
    return df_sanitized

def sanitize_smiles_single(row,db,file):
    smiles = row["SMILES"]
    name = row["name"]
    flag = row["flag"]
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as error:
        fail_message = "Failed at .csv input stage - sanitization error"
        general_functions.write_failed_smiles_to_db(smiles,db,file,fail_message)
        return
        
    if mol is None:
        fail_message = "Failed at .csv input stage - mol generation at sanitization error"
        general_functions.write_failed_smiles_to_db(smiles,db,file,fail_message)
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

def compute_inchi_key_refactored(row,db):
    try:
        smiles = row["SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        inchi_key = Chem.MolToInchiKey(mol)
        
        return inchi_key
    
    except:
        return None
    
def compute_inchi_key(smiles,db,file):
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchi_key = Chem.MolToInchiKey(mol)
        return inchi_key
        
    except:
        fail_message = "Failed at .csv input stage - InChI_Key computation error"
        general_functions.write_failed_smiles_to_db(smiles,db,file,fail_message)

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

def depict_ligands_table(db,table_name,output_path,max_mols_ppage,limit,random):
    conn = tidyscreen.connect_to_db(db)
    sql=f"SELECT SMILES, inchi_key, name FROM {table_name};"
    molecules_df = pd.read_sql_query(sql,conn)

    if limit > 0:
        # If a limit is provided, take only the first 'limit' rows
        if random == True:
            molecules_df = molecules_df.sample(limit)
        else:
            molecules_df = molecules_df.head(limit)

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

def process_all_mols_in_table(db,table_name,charge,pdb,mol2,pdbqt,temp_dir):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    
    list_mol_objects_colnames, list_mol_objects_colnames_types = create_mols_columns_from_selection(pdb,mol2,pdbqt)
        
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
    df.parallel_apply(lambda row: append_ligand_mols_blob_object_to_table(db,table_name,row,charge,pdb,mol2_sybyl,mol2_gaff2,pdbqt,temp_dir), axis=1)

    try:
        pass
        #clean_temp_dir(db,table_name)
    except:
        pass

def create_mols_columns_from_selection(pdb,mol2,pdbqt):
    # Define the list of column name to be created to store mol objects based on output selection
    list_mol_objects_colnames = []
    list_mol_objects_colnames_types = []

    if pdb == 1:
        list_mol_objects_colnames.append("pdb_file")
        list_mol_objects_colnames_types.append("BLOB")
    if mol2 == 1:
        list_mol_objects_colnames.append("mol2_file_sybyl")
        list_mol_objects_colnames_types.append("BLOB")
        list_mol_objects_colnames.append("mol2_file_gaff")
        list_mol_objects_colnames_types.append("BLOB")
        list_mol_objects_colnames.append("frcmod_file")
        list_mol_objects_colnames_types.append("BLOB")
    if pdbqt == 1:
        list_mol_objects_colnames.append("pdbqt_file")
        list_mol_objects_colnames_types.append("BLOB")

    # Finally add the charge model if mol2 are requested:
    
    if mol2 == 1:
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

def append_ligand_mols_blob_object_to_table(db,table_name,row,charge,pdb,mol2_sybyl,mol2_gaff2,pdbqt,temp_dir):
    
    action = 0 
    
    if pdb == 1:
        # Prepare and store the corresponding .pdb file
        pdb_file, tar_pdb_file, net_charge = pdb_from_smiles(row["SMILES"],temp_dir)
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

def pdb_from_smiles(smiles,dir,conf_rank):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol_hs = Chem.AddHs(mol)
        confs, mol_hs, ps = generate_ligand_conformers(mol_hs)
        selected_mol = get_conformer_rank(confs, mol_hs, ps, conf_rank)
        pdb_file, inchi_key = save_ligand_pdb_file(selected_mol,dir)
        # Compute the net charge of the molecule for potential use in antechamber
        net_charge = compute_molecule_net_charge(mol)
            
        # Compress the pdb_file
        tar_pdb_file = generate_tar_file(pdb_file)
        
        return pdb_file, tar_pdb_file, net_charge
    
    except Exception as error:
        print(f"Problem with SMILES: \n {smiles}")

def pdb_from_smiles_safetime(smiles,dir):
    start = time.time()
    try:
        while True:
            print("COMPUTANDO")
            mol = Chem.MolFromSmiles(smiles)
            mol_hs = Chem.AddHs(mol)
            confs, mol_hs, ps = generate_ligand_conformers(mol_hs)
            print("SALI")
            selected_mol = get_conformer_rank(confs, mol_hs, ps)
            pdb_file, inchi_key = save_ligand_pdb_file(selected_mol,dir)
            # Compute the net charge of the molecule for potential use in antechamber
            net_charge = compute_molecule_net_charge(mol)
                
            # Compress the pdb_file
            tar_pdb_file = generate_tar_file(pdb_file)
            
            if time.time() - start > 10:
                print("Timeout reached for row:", smiles)
                return None
            
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

def sybyl_mol2_from_pdb(inchi_key,charge,temp_dir):
    
    pdb_file = f"{temp_dir}/{inchi_key}.pdb"
    output_file = f"{temp_dir}/{inchi_key}_sybyl.mol2"
    # Since Antechamber creates temporary files for the calculation, and they should not overwrite each other when executing in parallel, create a temporary folder following the molecule InChIKey, 'cd' into it and run the computation.
    antechamber_temp_folder = f"{temp_dir}/{inchi_key}"
    os.makedirs(antechamber_temp_folder,exist_ok=True) # Create the corresponding temp directory
    # Compute the ligand net charge
    net_charge = compute_charge_from_pdb(pdb_file)
    # Compute mol2 with Sybyl Atom Types - for compatibility with RDKit and Meeko
    antechamber_path = shutil.which('antechamber')
    try: 
        antechamber_command = f'cd {antechamber_temp_folder} && {antechamber_path} -i {pdb_file} -fi pdb -o {output_file} -fo mol2 -c {charge} -nc {net_charge} -at sybyl -pf y' # The 'sybyl' atom type convention is used for compatibility with RDKit
        subprocess.run(antechamber_command, shell=True, capture_output=True, text=True)

        # Once computation finished, delete the ligand temporary folder
        shutil.rmtree(antechamber_temp_folder)        
            
    except Exception as error:
        print(error)
    # # Compress the output file for storage
    output_tar_file = generate_tar_file(output_file)
        
    return output_file, output_tar_file

def gaff_mol2_from_pdb(inchi_key,charge,temp_dir):
    
    pdb_file = f"{temp_dir}/{inchi_key}.pdb"
    output_file = f"{temp_dir}/{inchi_key}_gaff.mol2"
    # Since Antechamber creates temporary files for the calculation, and they should not overwrite each other when executing in parallel, create a temporary folder following the molecule InChIKey, 'cd' into it and run the computation.
    antechamber_temp_folder = f"{temp_dir}/{inchi_key}"
    os.makedirs(antechamber_temp_folder,exist_ok=True) # Create the corresponding temp directory
    # Compute the ligand net charge
    net_charge = compute_charge_from_pdb(pdb_file)
    # Compute mol2 with Sybyl Atom Types - for compatibility with RDKit and Meeko
    antechamber_path = shutil.which('antechamber')
    try: 
        antechamber_command = f'cd {antechamber_temp_folder} && {antechamber_path} -i {pdb_file} -fi pdb -o {output_file} -fo mol2 -c {charge} -nc {net_charge} -at gaff -pf y' # The 'sybyl' atom type convention is used for compatibility with RDKit
        
        subprocess.run(antechamber_command, shell=True, capture_output=True, text=True)

        # Once computation finished, delete the ligand temporary folder
        shutil.rmtree(antechamber_temp_folder)        
            
    except Exception as error:
        print(error)
    
    
    # # Compress the output file for storage
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
    confs = AllChem.EmbedMultipleConfs(mol,nbr_confs,props)
    ps = AllChem.MMFFGetMoleculeProperties(mol)
    AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant='MMFF94s', maxIters=maxIters)
    
    return confs, mol, ps

def get_conformer_rank(confs, mol_hs, ps,conf_rank):
    
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

    # Compute the position of the conformer to be selected based of the percentage provided as 'conf_rank'
    conf_rank_position = int(len(conformers_energies_dict_sorted) * conf_rank / 100)    


    # The following will store the lowest energy conformer as .pdbqt
    selected_conformer = Chem.MolToMolBlock(mol_hs,confId=conformers_energies_dict_sorted[conf_rank_position][0])
    selected_mol = Chem.MolFromMolBlock(selected_conformer, removeHs=False)

    return selected_mol

def save_ligand_pdb_file(selected_mol,dir="/tmp"):
    inchi_key = Chem.MolToInchiKey(selected_mol)
    pdb_file = f'{dir}/{inchi_key}.pdb'
    
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
    
def create_mols_columns(conn, table_name,list_mol_objects_colnames,list_mol_objects_colnames_types):
    # Create all the columns to store mol objects
    try:
        cursor = conn.cursor()
        for index, column in enumerate(list_mol_objects_colnames):
            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column} {list_mol_objects_colnames_types[index]};")
    except Exception as error:
        print(error)
        print("Stopping")
        sys.exit()
    pass

def compute_and_store_pdb(row,db,table_name,temp_dir,conf_rank,timeout):
    # Prepare and store the corresponding .pdb file
    try:
        pdb_file, tar_pdb_file, net_charge = general_functions.timeout_function(pdb_from_smiles,(row["SMILES"],temp_dir,conf_rank),timeout=timeout,on_timeout_args=[row["SMILES"],db,table_name,"Timed out at: 'pdb_from_smiles' computation"])
        store_file_as_blob(db,table_name,'pdb_file',tar_pdb_file,row)
    except Exception as error:
        fail_message = f"Failed at .pdb molecule computation step - Mol id: {row['id']} - deleted from {table_name}"
        general_functions.write_failed_smiles_to_db(row["SMILES"],db,table_name,fail_message)

def compute_and_store_pdb_2(row,db,table_name,temp_dir):
    # Prepare and store the corresponding .pdb file
    try:
        pdb_file, tar_pdb_file, net_charge = pdb_from_smiles(row["SMILES"],temp_dir)
        store_file_as_blob(db,table_name,'pdb_file',tar_pdb_file,row)
    except Exception as error:
        fail_message = f"Failed at .pdb molecule computation step - Mol id: {row['id']}"
        general_functions.write_failed_smiles_to_db(row["SMILES"],db,table_name,fail_message)

def compute_and_store_pdb_safetime(row,db,table_name,temp_dir):
    start = time.time()
    try:
        while True:
            #pdb_file, tar_pdb_file, net_charge = pdb_from_smiles(row["SMILES"],temp_dir)
            pdb_file, tar_pdb_file, net_charge = pdb_from_smiles_safetime(row["SMILES"],temp_dir)
            store_file_as_blob(db,table_name,'pdb_file',tar_pdb_file,row)
            if time.time() - start > 1:
                print("Timed out by internal function")
    
    except Exception as error:
        fail_message = f"Failed at .pdb molecule computation step - Mol id: {row['id']}"
        general_functions.write_failed_smiles_to_db(row["SMILES"],db,table_name,fail_message)

def timeout_compute_and_store_pdb(row,db,table_name,temp_dir):
    return general_functions.timeout_function(compute_and_store_pdb_2,args=(row,db,table_name,temp_dir),timeout=10,on_timeout_args=[row["SMILES"],db,table_name,"Timed out at: 'pdb_from_smiles' computation"])

def compute_and_store_mol2(row,db,table_name,charge,temp_dir,atom_types_dict):
    try: 
        if charge != "bcc-ml":
            # Compute and store in the db the sybyl type mol2
            sybyl_mol2_output_file, sybyl_output_tar_file = sybyl_mol2_from_pdb(row["inchi_key"],charge,temp_dir)
            # Store the Sybyl like file
            store_file_as_blob(db,table_name,'mol2_file_sybyl',sybyl_output_tar_file,row)
            
            # Compute and store in the db the gaff type mol2
            gaff_mol2_output_file, gaff_output_tar_file = gaff_mol2_from_pdb(row["inchi_key"],charge,temp_dir)
            # Store the GAFF like file
            store_file_as_blob(db,table_name,'mol2_file_gaff',gaff_output_tar_file,row)
            # Compute the frcmod file
            output_tar_frcmod = compute_frcmod_file(gaff_mol2_output_file,at="gaff")
            # Store the frcmod file
            store_file_as_blob(db,table_name,'frcmod_file',output_tar_frcmod,row)
            # Store the charge model used
            store_string_in_column(db,table_name,"charge_model",charge,row)
        
        elif charge == "bcc-ml":
            # Compute the bbc-ml charges using the precomputed bbc-ml model
            charge_array = compute_bbc_ml_array(row["SMILES"],row["inchi_key"],temp_dir)
            # Compute the charge using a fake 'gas' model to be replaced computing SYBYL atom types
            sybyl_mol2_output_file, sybyl_output_tar_file = sybyl_mol2_from_pdb(row["inchi_key"],"gas",temp_dir)
            # Replace the charges on the .mol2 files using the precomputed bbc-ml charges
            sybyl_replaced_mol2_file = general_functions.replace_charge_on_mol2_file(sybyl_mol2_output_file,charge_array)
            sybyl_new_tar_file = generate_tar_file(sybyl_replaced_mol2_file)
            # Store the Sybyl like file
            store_file_as_blob(db,table_name,'mol2_file_sybyl',sybyl_new_tar_file,row)
        
            # Compute the charge using a fake 'gas' model to be replaced computing GAFF atom types
            gaff_mol2_output_file, gaff_output_tar_file = gaff_mol2_from_pdb(row["inchi_key"],"gas",temp_dir)                          
            # Replace the charges on the .mol2 files using the precomputed bbc-ml charges
            gaff_replaced_mol2_file = general_functions.replace_charge_on_mol2_file(gaff_mol2_output_file,charge_array)
            gaff_new_tar_file = generate_tar_file(gaff_replaced_mol2_file)
            # Store the GAFF like file
            store_file_as_blob(db,table_name,'mol2_file_gaff',gaff_new_tar_file,row)
            # Compute the frcmod file
            output_tar_frcmod = compute_frcmod_file(gaff_replaced_mol2_file,at="gaff")
            # Store the frcmod file
            store_file_as_blob(db,table_name,'frcmod_file',output_tar_frcmod,row)
            # Store the charge model used
            store_string_in_column(db,table_name,"charge_model",charge,row)
            
            
        else:
            print(f"Charge system: '{charge}' unknown. Stopping...")
            sys.exit()
            
          
    except Exception as error:
        fail_message = f"Failed at .mol2 sybyl atom type computation step - Mol id: {row['id']}"
        general_functions.write_failed_smiles_to_db(row["SMILES"],db,table_name,fail_message)

def compute_charge_from_pdb(pdb_file):
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False) 
    
    molecule_charge = 0
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        molecule_charge = molecule_charge + charge


    return molecule_charge

def rename_sybyl_to_gaff_mol2(row,db,table_name,temp_dir,atom_types_dict):
    sybyl_file = f"{temp_dir}/{row['inchi_key']}_sybyl.mol2"
    gaff_file = f"{temp_dir}/{row['inchi_key']}_gaff.mol2"
    
    with open(atom_types_dict, 'rb') as f:
        atoms_dict = pickle.load(f)
    
    #general_functions.sybyl_to_gaff_mol2(sybyl_file,gaff_file,atom_dict)
    mol2_gaff_file = general_functions.sybyl_to_gaff_mol2_moldf(sybyl_file,gaff_file,atoms_dict)
    
    output_tar_file = generate_tar_file(mol2_gaff_file)
    
    return mol2_gaff_file, output_tar_file
    
def compute_and_store_pdbqt(row,db,table_name,temp_dir):
    try:
        # Prepare and store the corresponding .pdbqt file
        tar_pdbqt_file = pdbqt_from_mol2(f"{temp_dir}/{row['inchi_key']}_sybyl.mol2")
        store_file_as_blob(db,table_name,'pdbqt_file',tar_pdbqt_file,row)
    except Exception as error:
        fail_message = f"Failed at .pdbqt molecule computation step - Mol id: {row['id']}"
        general_functions.write_failed_smiles_to_db(row["SMILES"],db,table_name,fail_message)    
    
def test_dummy_conformer(smiles,nbr_confs,maxIters):
    mol = Chem.MolFromSmiles(smiles)
    mol_hs = Chem.AddHs(mol)
    props = AllChem.ETKDG()
    props.pruneRmsThresh = 0.25
    props.useRandomCoords = True
    props.numThreads = 1
    confs = AllChem.EmbedMultipleConfs(mol_hs,nbr_confs,props)
    
def test_dummy_conformer_on_table(db,table_name):
    conn =  sqlite3.connect(db)
    df = pd.read_sql_query(f"SELECT SMILES FROM {table_name}", conn)
    conn.close()
    
    for idx, row in df.iterrows():
        try:
            run_with_hard_timeout(test_dummy_conformer(row["SMILES"],2,5))
        except Exception as error:
            fail_message = "Failed at dummy conformer generation"
            print(fail_message)
            general_functions.write_failed_smiles_to_db(row["SMILES"],db,table_name,fail_message)

def run_with_hard_timeout(func, args=(), kwargs=None, timeout=5):
    if kwargs is None:
        kwargs = {}
    p = multiprocessing.Process(target=func, args=args, kwargs=kwargs)
    p.start()
    p.join(timeout)
    
    if p.is_alive():
        print("Function timed out and will be killed.")
        p.terminate()
        p.join()
        return None
    
    print("termine")
    #return True  # or return a result via a multiprocessing.Queue
    
def compute_bbc_ml_array(smiles,inchi_key,temp_dir):
    
    molecule = Chem.MolFromSmiles(smiles)
    molecule_hs = Chem.AddHs(molecule)
    # Compute the array of bcc-ml charges using the ESPaloma package
    charge_array = charge(molecule_hs)
    
    try: 
        ## Write the atom, charge pairing to to a file
        with open(f"{temp_dir}/{inchi_key}_charges.txt",'w') as charge_outfile:
            counter = 0 
            for atom in molecule_hs.GetAtoms():
                record = f"{atom.GetSymbol()} : {atom.GetIdx()} {charge_array[counter]} \n"
                charge_outfile.write(record)
                counter += 1
    except Exception as error:
        print(error)
    
    return charge_array

def purge_failed_in_table(db,table_name):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Get values from the first table's column
    cursor.execute("SELECT smiles FROM failed_smiles")
    values = [row[0] for row in cursor.fetchall()]
    
    # Delete rows in the second table where the column matches those values
    for value in values:
        cursor.execute(f"DELETE FROM {table_name} WHERE SMILES = ?", (value,))

    conn.commit()
    conn.close()
    
def drop_mols_columns(db,table_name):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Create a new table without the unwanted columns
    cursor.execute(f"CREATE TABLE new_table AS SELECT id, SMILES, name, flag, stereo_nbr, inchi_key FROM {table_name}")
    
    # Drop the old table
    cursor.execute(f"DROP TABLE {table_name}")
    
    # 3. Rename the new table to the original name
    cursor.execute(f"ALTER TABLE new_table RENAME TO {table_name}")

    conn.commit()
    conn.close()
    
def create_properties_columns(db, table_name, properties_list):
    """
    Create the properties columns if they do not exist in the target table.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    for prop in properties_list:
        # Check if the column already exists
        cursor.execute(f"PRAGMA table_info({table_name})")
        columns = [column[1] for column in cursor.fetchall()]
        
        if prop not in columns:
            # If the column does not exist, create it
            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {prop} REAL")
    
    conn.commit()

def compute_properties(row, db, table_name, properties_list):
    # Create an empty dataframe to store the properties
    prop_df = pd.DataFrame(columns=properties_list)
    # Create the properties columns if they do not exist in the target table
    mol = Chem.MolFromSmiles(row["SMILES"])
    if mol is None:
        props_list = [None] * len(properties_list)
    
    props_list =  [getattr(Descriptors, prop)(mol) for prop in properties_list]
    
    # Store the properties in database 
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Properties columns to update
    set_clause = ', '.join([f"{prop} = ?" for prop in properties_list])
    
    try: 
        cursor.execute(f"UPDATE {table_name} SET {set_clause} WHERE id = ?",(*props_list, row["id"]))
        conn.commit()
    except Exception as error:
        print(f"Error updating properties for SMILES: {row['SMILES']}")
        print(error)    

def write_subset_record_to_db(db,table_name,type,prop_filter,description):
    """
    Write a subset of records to the database based on a property filter.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Create a new table for the subset
    cursor.execute(f"CREATE TABLE IF NOT EXISTS tables_subsets ( id INTEGER PRIMARY KEY AUTOINCREMENT, table_name TEXT, subset_name TEXT, fitering_type TEXT, prop_filter TEXT, description TEXT)")
    
    # Insert records into the subset table - This is the first addition to assign an id
    cursor.execute(f"INSERT INTO tables_subsets (table_name, subset_name, fitering_type, prop_filter, description) VALUES (?, ?, ?, ?, ?)",(table_name, f"{table_name}_subset", type, str(prop_filter), description))
    
    current_id = cursor.lastrowid
    
    # Step 2: Update the row with a string containing the current id
    current_subset = f"{table_name}_subset_{current_id}"
    cursor.execute("UPDATE tables_subsets SET subset_name = ? WHERE id = ?", (current_subset, current_id))
    
    conn.commit()
    conn.close()
    
    return current_subset

def subset_table_by_properties(db,table_name,current_subset,props_filter):
    """
    Subset the table by properties and store the results in a new table.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Create a new table for the subset
    cursor.execute(f"CREATE TABLE IF NOT EXISTS {current_subset} AS SELECT * FROM {table_name} WHERE 1=0")  # Create an empty table with the same structure
        
    # Build the WHERE clause from the filter list
    where_clause = " AND ".join(props_filter)
    
    # Insert records into the subset table
    cursor.execute(f"INSERT INTO {current_subset} SELECT * FROM {table_name} WHERE {where_clause}")
    
    conn.commit()
    conn.close()

def check_smarts_existence(db,smarts_filter):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    try: 
        # Retrieve the SMARTS filter from the database if exists
        cursor.execute("SELECT 1 FROM chem_filters WHERE smarts = ?", (smarts_filter,))
        exists = cursor.fetchone() is not None
        
        # Inform and stop if the SMARTS filter exists
        if exists:
            print(f"The filter: {smarts_filter} already exists in the environment 'chem_filters' table. Stopping...")
            sys.exit()
    
    except SystemExit: # This will be raised if the SMARTS filter already exists in the database
        raise
    
    except Exception as error: # This caught will be raised if the table does not exist yet
        print("SMARTS filter table does not exist yet. Creating it...")
        pass # The failure will pass if the table does not exist yet
    
def check_smarts_reaction_existence(db,smarts_reaction):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    try: 
        # Retrieve the SMARTS filter from the database if exists
        cursor.execute("SELECT 1 FROM smarts_reactions WHERE smarts = ?", (smarts_reaction,))
        exists = cursor.fetchone() is not None
        
        # Inform and stop if the SMARTS filter exists
        if exists:
            print(f"The reaction: {smarts_reaction} already exists in the project 'smarts_reactions' table. Stopping...")
            sys.exit()
    
    except SystemExit: # This will be raised if the SMARTS filter already exists in the database
        raise
    
    except Exception as error: # This caught will be raised if the table does not exist yet
        print("SMARTS reaction table does not exist yet. Creating it...")
        pass # The failure will pass if the table does not exist yet

def insert_smarts_filter_in_table(db,smarts_filter,description):
    """
    Insert a SMARTS filter into the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Create a new table for the SMARTS filters if it does not exist
    cursor.execute(f"CREATE TABLE IF NOT EXISTS chem_filters (id INTEGER PRIMARY KEY AUTOINCREMENT, smarts TEXT, description TEXT)")
    
    # Insert the SMARTS filter into the table
    cursor.execute("INSERT INTO chem_filters (smarts, description) VALUES (?,?)", (smarts_filter,description,))
    
    conn.commit()
    conn.close()
    
def parse_smarts_filters_dict(db_filters,db_workflow,smarts_filters_dict):
    
    # Create a dict of SMARTS filters to be used in the workflow
    conn = sqlite3.connect(db_filters)
    cursor = conn.cursor()
    
    cursor.execute("SELECT * FROM chem_filters")
    rows = cursor.fetchall()
    
    # Build the new dictionary
    filters_instances_dict = {}
    filters_names_list = []
    for row in rows:
        filter_id = row[0]
        if filter_id in smarts_filters_dict:
            filter_name = row[1]
            filter_smarts = row[2]
            filter_instances = smarts_filters_dict[filter_id]
            filters_instances_dict[filter_smarts] = filter_instances
            filters_names_list.append(filter_name)
    
    
    return filters_instances_dict, filters_names_list

def check_workflow_existence(db,filters_instances_dict):
    
    """
    Check if the workflow already exists in the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    pattern_to_check = json.dumps(filters_instances_dict)
    
    try: 
        # Check if the JSON string exists in the column 'dict_column' of 'my_table'
        cursor.execute("SELECT 1 FROM smarts_filters_workflow WHERE filters_instances = ?", (pattern_to_check,))
        exists = cursor.fetchone() is not None
        
        # Inform and stop if the SMARTS filter exists
        if exists:
            print(f"The filter: {filters_instances_dict} already exists as a workflow. Stopping...")
            sys.exit()
    
    except SystemExit: # This will be raised if the SMARTS filter already exists in the database
        raise
    
    except Exception as error: # This caught will be raised if the table does not exist yet
        print("SMARTS workflows table does not exist yet. Creating it...")
        pass # The failure will pass if the table does not exist yet

def check_smarts_reaction_workflow_existence(db,smarts_list):
    
    """
    Check if the workflow already exists in the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    pattern_to_check = json.dumps(smarts_list)
    
    try: 
        # Check if the JSON string exists in the column 'dict_column' of 'my_table'
        cursor.execute("SELECT 1 FROM smarts_reactions_workflow WHERE smarts_ids = ?", (pattern_to_check,))
        exists = cursor.fetchone() is not None
        
        # Inform and stop if the SMARTS filter exists
        if exists:
            print(f"The filter: {smarts_list} already exists as a smart reaction workflow. Stopping...")
            sys.exit()
    
    except SystemExit: # This will be raised if the SMARTS filter already exists in the database
        raise
    
    except Exception as error: # This caught will be raised if the table does not exist yet
        print("SMARTS workflows table does not exist yet. Creating it...")
        pass # The failure will pass if the table does not exist yet

def store_smarts_filters_workflow(db_workflow,filters_instances_dict,filters_names_list,smarts_filters_dict):
    """
    Store the SMARTS filters workflow in the database.
    """
    conn = sqlite3.connect(db_workflow)
    cursor = conn.cursor()
    
    description = input("Plovide a description for the SMARTS filters workflow: ")
    # Create a new table for the SMARTS filters workflow if it does not exist
    cursor.execute(f"CREATE TABLE IF NOT EXISTS smarts_filters_workflow (id INTEGER PRIMARY KEY AUTOINCREMENT, filters_instances TEXT, filters_names TEXT, filters_specification TEXT, description TEXT)")
    
    # Insert the SMARTS filters workflow into the table
    cursor.execute("INSERT INTO smarts_filters_workflow (filters_instances, filters_names, filters_specification, description) VALUES (?,?,?,?)", ((json.dumps(filters_instances_dict)),str(filters_names_list),(json.dumps(smarts_filters_dict)),description,))
    
    conn.commit()
    conn.close()
    
def retrieve_workflow_from_db(db,workflow_id):
    """
    Retrieve the SMARTS filters workflow from the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Retrieve the workflow from the table
    cursor.execute("SELECT filters_instances FROM smarts_filters_workflow WHERE id = ?", (workflow_id,))
    row = cursor.fetchone()
    
    if row is None:
        print(f"No SMARTS filter workflow found with ID: '{workflow_id}'. Stopping...")
        sys.exit()
    
    filters_instances_dict = json.loads(row[0])
    
    return filters_instances_dict

def subset_table_by_smarts_dict(db,table_name,filters_instances_dict):
    
    conn = sqlite3.connect(db)
    df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
    conn.close()    
    
    # Sort by values in decreasing order
    sorted_dict = dict(sorted(filters_instances_dict.items(), key=lambda item: item[1], reverse=True))
    
    for smarts, instances in sorted_dict.items():
        print(f"Filtering by SMARTS: '{smarts}' with '{instances}' instances")
        
        # Use pandarallel for parallel processing and count the instances of each SMARTS pattern
        pandarallel.initialize(progress_bar=True)
        df["smarts_instances"] = df['SMILES'].parallel_apply(lambda smiles: check_smiles_vs_smarts(smiles, smarts))
        # Delete rows where the count of instances is different than the specified number
        df = df[df["smarts_instances"] == instances]
        # Delete the instaces column for the next iteration
        df.drop(columns=["smarts_instances"], inplace=True)
    
    # Reset the index of the DataFrame
    df.reset_index(drop=True, inplace=True)
    
    return df

def check_smiles_vs_smarts(smiles, smarts):
    """
    Check if the SMILES string matches the SMARTS pattern.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    
    pattern = Chem.MolFromSmarts(smarts)
    
    if mol is not None and pattern is not None:
        matches = mol.GetSubstructMatches(pattern)
        matches_instances = len(matches)
        
        return matches_instances
    else:
        return 0
    
def list_available_smarts_filters(db):
    """
    List all available SMARTS filters in the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    cursor.execute("SELECT * FROM chem_filters")
    rows = cursor.fetchall()
    
    if not rows:
        print("No SMARTS filters found in the database.")
        return []
    
    print("Available SMARTS filters:")
    for row in rows:
        print(f"Filter_id: {row[0]}, Filter_Name: {row[1]}, SMARTS: {row[2]}")

def list_available_smarts_reactions(db):
    """
    List all available SMARTS filters in the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    cursor.execute("SELECT * FROM smarts_reactions")
    rows = cursor.fetchall()
    
    if not rows:
        print("No SMARTS reactions found in the database.")
        return []
    
    print("Available SMARTS reactions:")
    for row in rows:
        print(f"Reaction_id: {row[0]}, Reaction_SMARTS: {row[1]}, Description: {row[2]}")

def list_available_smarts_reactions_workflows(db,complete_info):
    """
    List all available SMARTS filters in the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    if complete_info == 0:
        cursor.execute("SELECT id, smarts_descriptions FROM smarts_reactions_workflow")
    else:
        cursor.execute("SELECT * FROM smarts_reactions_workflow")
    
    rows = cursor.fetchall()
    
    if not rows:
        print("No SMARTS reactions workflows found in the database.")
        return []
    
    print("Available SMARTS reactions:")
    for row in rows:
        if complete_info == 0:
            print(f"Reaction_workflow_id: {row[0]}, SMARTS_description: {row[1]}")
        else:
            print(f"Reaction_workflow_id: {row[0]}, reaction_SMARTS_ids: {row[1]}, SMARTS_description: {row[2]}, description: {row[3]}")
    
def insert_smarts_reaction_in_table(db,smarts_reaction,description):
    
    """
    Insert a SMARTS reaction into the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Create a new table for the SMARTS filters if it does not exist
    cursor.execute(f"CREATE TABLE IF NOT EXISTS smarts_reactions (id INTEGER PRIMARY KEY AUTOINCREMENT, smarts TEXT, description TEXT)")
    
    # Insert the SMARTS filter into the table
    cursor.execute("INSERT INTO smarts_reactions (smarts, description) VALUES (?,?)", (smarts_reaction,description,))
    
    conn.commit()
    conn.close()

def parse_smarts_reactions_id_list(db,smarts_reactions_id_list):
    """
    Parse the list of SMARTS reactions IDs and return a dictionary with the SMARTS reactions.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Retrieve the SMARTS reactions from the database
    cursor.execute("SELECT * FROM smarts_reactions")
    rows = cursor.fetchall()
    
    smarts_reactions_list = []
    smarts_descriptions_list = []
    
    for reaction_id in smarts_reactions_id_list:
        for row in rows:
            if row[0] == reaction_id:
                # Append the SMARTS reaction to the list
                smarts_reactions_list.append(row[1])
                smarts_descriptions_list.append(row[2])
    
    return smarts_reactions_list, smarts_descriptions_list

def store_smarts_reactions_workflow(db,smarts_reactions_list,smarts_reactions_id_list,smarts_descriptions_list,description):
    """
    Add a SMARTS reaction workflow to the database.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Create a new table for the SMARTS reactions workflow if it does not exist
    cursor.execute(f"CREATE TABLE IF NOT EXISTS smarts_reactions_workflow (id INTEGER PRIMARY KEY AUTOINCREMENT, smarts_reactions TEXT, smarts_ids TEXT, smarts_descriptions, description TEXT)")
    
    # Insert the SMARTS reactions workflow into the table
    cursor.execute("INSERT INTO smarts_reactions_workflow (smarts_reactions, smarts_ids, smarts_descriptions, description) VALUES (?,?,?,?)", (json.dumps(smarts_reactions_list),json.dumps(smarts_reactions_id_list),json.dumps(smarts_descriptions_list),description,))
    
    conn.commit()
    conn.close()
    
def create_reactions_results_table(db):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Create a new table for the SMARTS reactions workflow if it does not exist
    cursor.execute(f"CREATE TABLE IF NOT EXISTS smarts_reactions_results (id INTEGER PRIMARY KEY AUTOINCREMENT, reaction_workflow_id INTEGER, reactants_lists TEXT, results_table TEXT)")
    
    conn.commit()
    conn.close()
    
def check_reaction_definition(db,reaction_workflow_id,reactants_lists):
    
    try: 
        # Retrieve the reaction workflow from the database
        smarts_reaction_workflow_list = retrieve_reaction_smarts_definition(db,reaction_workflow_id)
        
        # Check if the lists of lists of reactants are constistant with the reaction workflow definition
        check_workflows_vs_reactants_lists(smarts_reaction_workflow_list,reactants_lists)
        
        return smarts_reaction_workflow_list
    
    except:
        print("Writting failed attempt record to the database...")
        write_reaction_attempt_record_to_db(db,smarts_reaction_workflow_list,"Error checking the reaction workflow definition vs the reactants lists. Stopping...")
        sys.exit()
        
def retrieve_reaction_smarts_definition(db,reaction_workflow_id):
    
    # get the corresponding reaction workflow from the database
    conn = sqlite3.connect(db)
    cursor = conn.cursor()  
    cursor.execute("SELECT smarts_reactions FROM smarts_reactions_workflow WHERE id = ?", (reaction_workflow_id,))
    
    row = cursor.fetchone()[0]
    
    # Convert the smarts_reaction_workflow to a list
    smarts_reaction_workflow_list = json.loads(row)
    
    return smarts_reaction_workflow_list
    
def check_workflows_vs_reactants_lists(smarts_reaction_workflow_list,reactants_lists):
    """
    Check if the lists of lists of reactants are consistant with the reaction workflow definition.
    """
    
    # Check if the number of reactants in each list is consistent with the reaction workflow definition
    
    try: 
    
        for index, smarts_reaction_workflow_step in enumerate(smarts_reaction_workflow_list):
            nbr_required_reactants = smarts_reaction_workflow_step.split('>>')[0].count(".") + 1   # Count the number of '.' in the SMARTS string to determine the number of reactants
            
            if nbr_required_reactants != len(reactants_lists[index]):
                print(f"The number of reactants in the list {index} is not consistent with the reaction workflow definition. Stopping...")
                
                sys.exit()
                
        # Inform successful termination of reactino steps checking    
        print("Checking the reaction workflow definition vs the reactants lists: OK")
    
    except:
        print("Error checking the reaction workflow definition vs the reactants lists. Stopping...")
        sys.exit()

def apply_validated_reaction_workflow(db,smarts_reaction_workflow,reactants_lists,remove_mull_inchi_key=True):
    
    # Create an empty dataframe to store the results of the reaction steps
    final_products_df = pd.DataFrame()
        
    for index, reaction in enumerate(smarts_reaction_workflow):
        
        reaction_molecularity = reaction.split('>>')[0].count(".") + 1
        
        if reaction_molecularity == 1:
            # Check if the reactant for the current step is originated in the last reaction step
            if "->" not in reactants_lists[index][0]:
                # Retieve the reactants as a dataframe
                reactants_df = retrieve_smiles_reactants_as_df(db,reactants_lists[index])
                
            else:
                # Parse the reactants_lists to get the previous step reference
                previous_step_reference = abs(int(reactants_lists[index][0].split(":")[-1]))
                # Prepare the column name to retrieve the products from the previous step
                colname = f"SMILES_product_step_{index - previous_step_reference}"
                reactants_serie = final_products_df[colname] # return a pandas series with the SMILES of the products from the previous step
                # Convert the resulting series to a DataFrame and assign "SMILES" as column name
                reactants_df = reactants_serie.to_frame(name='SMILES')
                
            
            # Apply the unimolecular reaction step 
            df_current_step_products = apply_single_unimolecular_reaction_step(db,smarts_reaction_workflow,reaction,reactants_df,index)
        
            # Concatenate the products from the current step to the final products dataframe
            final_products_df = pd.concat([final_products_df, df_current_step_products], ignore_index=True)
        
        
        elif reaction_molecularity == 2:
            # Bimolecular reaction
            
            # Check if the reactant for the current step is originated in the last reaction step
            if reactants_lists[index][0] != "->" and reactants_lists[index][1] != "->":
                
                # Retieve the reactants as a dataframe
                reactants_df1 = retrieve_smiles_reactants_as_df(db,reactants_lists[index],index=0)
                reactants_df2 = retrieve_smiles_reactants_as_df(db,reactants_lists[index],index=1)
                
            else:
                colname = f"SMILES_product_step_{index-1}"
                
                reactants_serie = workflow_df[colname] # return a pandas series with the SMILES of the products from the previous step
                # Convert the resulting series to a DataFrame and assign "SMILES" as column name
                reactants_df = reactants_serie.to_frame(name='SMILES')
            
            
            # Apply the unimolecular reaction step 
            workflow_df = apply_single_bimolecular_reaction_step(db,smarts_reaction_workflow,reaction,reactants_df1,reactants_df2,index)
            
        
        else:
            print(f"Reaction {reaction} is not a unimolecular or bimolecular reaction. Stopping...")
            write_reaction_attempt_record_to_db(db,smarts_reaction_workflow,"Reaction is not a unimolecular or bimolecular reaction. Stopping...")    
            sys.exit()
        
       
    # Once the reaction workflow is applied, get the final products from the last step
    final_products_df = final_products_df.iloc[:,[-1]].rename(columns={final_products_df.columns[-1]: 'SMILES'})
    
    # Assign the name ("synth") to the final products
    final_products_df['name'] = 'synth'
    
    # Assign default flag to the final products
    final_products_df['flag'] = 0
    
    # Compute the inchi_key for the final products
    pandarallel.initialize(progress_bar=False)
    final_products_df["inchi_key"] = final_products_df.parallel_apply(lambda row: compute_inchi_key_refactored(row,db),axis=1)
    
    if remove_mull_inchi_key:
        # If the calculation of the inchi_key failed, remove the corresponding rows
        final_products_df = final_products_df[final_products_df["inchi_key"].notnull()]
    
    return final_products_df
     
def apply_single_unimolecular_reaction_step(db,smarts_reaction_workflow,smarts_reaction,reactants_df,workflow_step_index):
    
    # Apply in a parallelized scheme the corresponding reaction
    pandarallel.initialize(progress_bar=False)
    
    df_current_step_products = pd.DataFrame()
    product_column = f"SMILES_product_step_{workflow_step_index}"
    
    df_current_step_products[product_column] = reactants_df['SMILES'].parallel_apply(lambda smiles: react_molecule_1_component(smiles, smarts_reaction))
    
    # Check if the reaction was successful
    if df_current_step_products[product_column].isnull().all():
        print(f"All reactants in the step {workflow_step_index} failed to react. Stopping...")
        
        write_reaction_attempt_record_to_db(db,smarts_reaction_workflow,"Reaction was not succesfull. Stopping...")
        
        sys.exit()
    
    else:
        # Inform successful termination of the reaction step
        print(f"Reaction step {workflow_step_index} applied successfully. Products stored in column: '{product_column}'")
        
        return df_current_step_products

def apply_single_bimolecular_reaction_step(db,smarts_reaction_workflow,smarts_reaction,reactants_df1,reactants_df2,workflow_step_index):
    pandarallel.initialize(progress_bar=False)
    df_all_products = pd.DataFrame()
    df_current_step_products = pd.DataFrame()
    product_column = f"product_step_{workflow_step_index}"
    
    # In order to combinatorialy prepare the products, the loop on the dataframe should be performed on the one containing the lowest number of compounds, otherwise it will fail. So here the length of both dataframes are compared, and actioned in accordance
            
    if len(reactants_df1) < len(reactants_df2): 
        for index, row in reactants_df1.iterrows():
            current_smiles_1 = row["SMILES"]
            # Apply the reaction on the first reactant
            df_current_step_products[f"SMILES_{product_column}"] = reactants_df2["SMILES"].parallel_apply(lambda smiles_2:react_molecule_2_components(current_smiles_1,smiles_2,smarts_reaction))
            
            # Concatenate the products from the current step to the final products dataframe
            df_all_products = pd.concat([df_all_products, df_current_step_products], ignore_index=True)
    
    else: # The df2 contains less reactants that df1, so loop over it
    
        for index, row in reactants_df2.iterrows():
            current_smiles_2 = row["SMILES"]
            
            # Apply the reaction on the second reactant
            df_current_step_products[f"SMILES_{product_column}"] = reactants_df1["SMILES"].parallel_apply(lambda smiles_1:react_molecule_2_components(smiles_1,current_smiles_2,smarts_reaction))
    
            # Concatenate the products from the current step to the final products dataframe
            df_all_products = pd.concat([df_all_products, df_current_step_products], ignore_index=True)
    
    print("apply_single_bimolecular_reaction_step")
    
    return df_all_products
    
def react_molecule_1_component(smiles,reaction_smarts):
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        reaction_mol = AllChem.ReactionFromSmarts(reaction_smarts)
        product_mol = reaction_mol.RunReactants((mol,))[0][0]
        product_smiles = Chem.MolToSmiles(product_mol)
    
        return product_smiles
    
    except Exception as error:
        pass

def react_molecule_2_components(smiles1, smiles2,reaction_smarts):
    
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        reaction_mol = AllChem.ReactionFromSmarts(reaction_smarts)
        product_mol = reaction_mol.RunReactants((mol1,mol2))[0][0]
        product_smiles = Chem.MolToSmiles(product_mol)
    
        return product_smiles
    
    except Exception as error:
        pass

def retrieve_smiles_reactants_as_df(db,reactants_list,index=0):
    
    conn= sqlite3.connect(db)
    
    df = pd.read_sql_query(f"SELECT SMILES FROM {reactants_list[index]}", conn)
    conn.close()
    
    return df
    
def store_reaction_results(db, final_products_df, reaction_workflow_id, reactants_lists):
    
    table_name = write_reaction_attempt_record_to_db(db,reaction_workflow_id)
        
    save_df_to_db(db,final_products_df,table_name)
    

 # Generate a reaction attempt record in the database
def write_reaction_attempt_record_to_db(db,reaction_workflow_id,message="Reaction attempt successful"):
    
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    # Create a new table for the subset
    cursor.execute(f"CREATE TABLE IF NOT EXISTS reactions_attempts ( id INTEGER PRIMARY KEY AUTOINCREMENT, reaction_workflow_id INTEGER, table_name TEXT, result TEXT)")
    
    # Get the last record of 'my_column'
    cursor.execute("SELECT id FROM reactions_attempts ORDER BY id DESC LIMIT 1")
    result = cursor.fetchone()
    
    if result is None:
        last_id = 0
    else:
        last_id = result[0]     
    
    if message == "Reaction attempt successful":
        table_name = f"reaction_set_{last_id + 1}"
    else:
        table_name = "None"
    
    # Insert a new record with the incremented id
    cursor.execute("INSERT INTO reactions_attempts (id, reaction_workflow_id, table_name, result) VALUES (?, ?, ?, ?)", (last_id + 1, reaction_workflow_id, table_name, message))
    conn.commit()
    conn.close()
    print(f"Reaction attempt record written to the database with id: {last_id + 1}")
    
    return table_name
    
    