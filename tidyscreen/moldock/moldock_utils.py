import tarfile
from tidyscreen import tidyscreen as tidyscreen
import sys
import os
import json
from glob import glob
import subprocess
from ringtail import RingtailCore
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def tar_folder(folder_path,file_prefix):
    # Prepare storing naming
    tar_filename = f"{file_prefix}.tar"
    output_tar = f"/tmp/{tar_filename}"
    
    with tarfile.open(output_tar, "w") as tar:
        tar.add(folder_path, arcname="")  # Keeps original folder name inside tar
    
    # Create a BLOB object with the .tar file 
    with open(output_tar, 'rb') as file:
        blobData = file.read()
    
    file.close()
        
    return tar_filename, blobData

def store_receptor_model(db,tar_filename, receptor_blob):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    
    # Indicate a description for the receptor model
    description = input("Provide a brief description of the receptor model: ")
    
    # Try to create the 'receptors' table if it does not exist
    cursor.execute("CREATE TABLE IF NOT EXISTS receptors (id INTEGER PRIMARY KEY AUTOINCREMENT, filename TEXT, receptor_tar_file BLOB, description TEXT)")
    conn.commit()
    
    # Store the BLOB object in the database
    cursor.execute(f"INSERT INTO receptors (filename, receptor_tar_file, description) VALUES (?,?,?)",(tar_filename,receptor_blob,description))
    conn.commit()
    conn.close()
    
    print("Succesfully stored receptor model")
    
def check_existing_rec_model(db,folder):
    model_to_check = folder.split('/')[-1].replace('-','_') + '.tar'
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    
    try:
        # Query to check if the value exists
        cursor.execute("SELECT EXISTS(SELECT 1 FROM receptors WHERE filename = ?)", (model_to_check,))
        exists = cursor.fetchone()[0]  # 1 if exists, 0 if not
        conn.close()
    
        if exists == 1:
            #print(f"The receptor model: '{model_to_check}' already exists. Stopping.")
            option = input(f"The receptor model: '{model_to_check}' already exists. Continue? (y/n): ")
            if option == 'y':
                print("Continuing...")
            else: # Will raise a valueError to stop the process
                raise ValueError

    except ValueError:
        print("Stopping...")
        sys.exit()    
    
    except Exception as error: # This 'pass' allows the processing for the first receptor, since the first check will due to the table 'receptors' not existing.
        pass

def check_default_docking_conditions(docking_params_db):
    
    try:
        conn = tidyscreen.connect_to_db(docking_params_db)
        cursor = conn.cursor()
    
        default_conditions_description = "default docking conditions"
    
        # Query to check if the string exists
        cursor.execute(f"SELECT description FROM docking_params WHERE description = ?", (default_conditions_description,))
        result = cursor.fetchone()
    
        conn.close()
    
        if result:
            return 1
        else:
            return 0
    
    except:
        return 0 # Will return 0 if the database does not exist

def create_default_params_register(docking_params_db):
    
    general_opts_dict = {'files_ops':{},
                 'conversion_ops':{},
                 'output_ops':{},
                 'setup_ops':{},
                 'search_ops':{},
                 'scoring_ops':{},
                 }

    files_ops_dict = {'--lfile':'LFILE',
                       '--ffile':'FLDFILE',
                       '--flexres':'FLEXFILE',
                       '--filelist':'BATCHFILE',
                       '--import_dpf':"DPFFILE",
                       '--xraylfile': "XRAYLIGFILE",
                       }

    conversion_ops_dict = {'--xml2dlg':'None',
                        }
   
    output_ops_dict = {'--resname':'ligand basename',
                       '--contact_analysis':0,
                       '--xmloutput':1,
                       '--dlgoutput':1,
                       '--dlg2stdout':0,
                       '--rlige':0,
                       '--gfpop':0,
                       '--npdb':0,
                       '--gbest':0,
                       '--clustering':1,
                       '--output-cluster-poses':0,
                       '--hsym':1,
                       '--rmstol':2,
                       }

    setup_ops_dict = {'--devnum':1,
                      '--loadxml':None,
                      '--seed':'time',
                     }
    
    search_ops_dict = {'--heuristics':1,
                       '--heurmax':12000000,
                       '--autostop':1,
                       '--asfreq':5,
                       '--nrun':20,
                       '--nev':2500000,
                       '--ngen':42000,
                       '--lsmet':'ad',
                       '--lsit':300,
                       '--psize':150,
                       '--mrat':2,
                       '--crat':80,
                       '--lsrat':100,
                       '--trat':60,
                       '--dmov':6,
                       '--dang':90,
                       '--rholb':0.01,
                       '--lsmov':2,
                       '--lsang':75,
                       '--cslim':4,
                       '--stopstd':0.15,
                       '--initswgens':0
                       }

    scoring_ops_dict = {'--derivtype':None,
                        '--modpair':None,
                        '--ubmod':0,
                        '--smooth':0.5,
                        '--elecmindist':0.01,
                        '--modqp':0,
                        }
    
    # Configure the general options dictionary with the custom values
    general_opts_dict['files_ops'] = files_ops_dict
    general_opts_dict['conversion_ops'] = conversion_ops_dict
    general_opts_dict['output_ops'] = output_ops_dict
    general_opts_dict['setup_ops'] = setup_ops_dict
    general_opts_dict['search_ops'] = search_ops_dict
    general_opts_dict['scoring_ops'] = scoring_ops_dict
    
    # Params dictionaty description
    description = "default docking conditions"
    
    options = json.dumps(general_opts_dict)
        
    return options, description
    
def compare_conditions_to_default(params_db,conditions_dict):

    default_options = create_default_params_register(params_db)
    default_options_dict = json.loads(default_options[0])

    # This will detect parameters differences between 
    custom_parameter_string = ''
    
    for key in default_options_dict.keys():
        for sub_key in default_options_dict[key].keys():
            if default_options_dict[key][sub_key] != conditions_dict[key][sub_key]:
                custom_parameter_string = custom_parameter_string + f'{sub_key} {conditions_dict[key][sub_key]} '

    return custom_parameter_string

def store_docking_params_register(docking_params_db, options, description):
    conn = tidyscreen.connect_to_db(docking_params_db)
    cursor = conn.cursor()
    
    if description != "default docking conditions":
        description = input("Provide a brief description for the parameter set: ")
        
    # Try to create the registers table in case it does no exists
    cursor.execute("CREATE TABLE IF NOT EXISTS docking_params (id INTEGER PRIMARY KEY AUTOINCREMENT, params_dict TEXT, description TEXT)")
    conn.commit()
    
    # Store the docking parameters register in the database
    cursor.execute(f"INSERT INTO docking_params (params_dict, description) VALUES (?,?)",(options,description,))
    conn.commit()
    
def append_docking_registry(db,table_name,id_receptor_model,id_docking_params):
    # Try to create the 'docking_registries' table if it does not exist
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    description = input("Provide a brief description of the docking assay: ")
    
    # Try to create the registers table in case it does no exists
    cursor.execute("CREATE TABLE IF NOT EXISTS docking_registries (id INTEGER PRIMARY KEY AUTOINCREMENT, table_name TEXT, id_receptor_model INTEGER, id_docking_params INTEGER, description TEXT)")
    conn.commit()
    
    # Store the docking register in the database
    cursor.execute(f"INSERT INTO docking_registries (table_name, id_receptor_model, id_docking_params, description) VALUES (?,?,?,?)",(table_name,id_receptor_model,id_docking_params,description))
    conn.commit()
    
    # Gest the id of the newly create register
    cursor.execute("SELECT id FROM docking_registries ORDER BY id DESC LIMIT 1")
    last_assay_id = cursor.fetchone()[0]
    conn.close()
    
    return last_assay_id
    
def create_assay_folder(assay_folder):
    os.makedirs(f"{assay_folder}/receptor",exist_ok=True)
    os.makedirs(f"{assay_folder}/ligands",exist_ok=True)
    os.makedirs(f"{assay_folder}/dlgs",exist_ok=True)

def retrieve_ligands_for_docking(ligands_db,assay_folder,table_name):
    conn = tidyscreen.connect_to_db(ligands_db)
    cursor = conn.cursor()
    cursor.execute(f"SELECT inchi_key, pdbqt_file FROM {table_name}")
    rows = cursor.fetchall()
    
    # Save each BLOB as a file
    
    counter = 1
    
    for row in rows:
        inchi_key = row[0]
        blob_data = row[1]

        # Define output file path 
        output_file = os.path.join(f"{assay_folder}/ligands", f"{inchi_key}.tar")
        
        # Write the BLOB object
        with open(output_file, "wb") as f:
            f.write(blob_data)
        
        # Decompress the ligand tar file
        with tarfile.open(output_file, "r:*") as tar:
            tar.extractall(path=f"{assay_folder}/ligands")
        
        # Delete the tar file
        os.remove(output_file)
        
    conn.close()
    
def retrieve_receptor(assay_folder,receptor_db,id_receptor_model):
    conn = tidyscreen.connect_to_db(receptor_db)
    cursor = conn.cursor()
    cursor.execute(f"SELECT filename, receptor_tar_file FROM receptors WHERE id = {id_receptor_model}")
    rows = cursor.fetchall()
    
    for row in rows:
        filename = row[0]
        blob_data = row[1]
    
        # Define output file path 
        output_file = os.path.join(f"{assay_folder}/receptor", filename)
        
        # Write the BLOB object
        with open(output_file, "wb") as f:
            f.write(blob_data)
            
        # Decompress the receptor tar file
        with tarfile.open(output_file, "r:*") as tar:
            tar.extractall(path=f"{assay_folder}/receptor")
    
        # Delete the tar file
        os.remove(output_file)

def retrieve_docking_conditions(params_db,id_docking_params):
    conn = tidyscreen.connect_to_db(params_db)
    cursor = conn.cursor()
    cursor.execute(f"SELECT params_dict FROM docking_params WHERE id = {id_docking_params}")
    row = cursor.fetchone()[0]
    
    conditions_dict = json.loads(row)
    
    return conditions_dict

def create_docking_executable(assay_folder,custom_parameter_string, nchunk=1000):
    # Get the list of ligands to dock
    list_of_ligands = glob(f'{assay_folder}/ligands/*.pdbqt')
    # Get the .fld file to configure docking assays
    fld_file = glob(f'{assay_folder}/receptor/*.fld')[0].split('/')[-1]
    
    # Divide the originial lists in sublists 
    divided_list = [list_of_ligands[i * nchunk:(i + 1) * nchunk] for i in range((len(list_of_ligands) + nchunk - 1) // nchunk )]  
    chunk = 1
    
    # Append the execution of docking assays for each ligand
    for list_of_ligands in divided_list:
        counter = 1
        with open(f'{assay_folder}/docking_execution_{chunk}.sh','a') as exec_file:
            for ligand in list_of_ligands:
                # Append instruction at start of file to register initial time
                if counter == 1:
                    exec_file.write("start=$(date +%s)\n")
                    exec_file.write(f"date > timing_{chunk}.txt\n")
                    counter += 1 # This will exit for the loop for this iteration
    
                ligand_file = ligand.split('/')[-1]
                ligand_prefix = ligand_file.replace('.pdbqt','')
                exec_file.write(f'autodock_gpu_128wi --lfile ./ligands/{ligand_file} --ffile ./receptor/{fld_file} {custom_parameter_string}\n')
                exec_file.write(f'mv ./ligands/{ligand_prefix}.dlg ./dlgs\n')
    
            # Append the final timings instruction
            exec_file.write("end=$(date +%s)\n")
            exec_file.write(f"date >> timing_{chunk}.txt\n")
            exec_file.write(f'echo "It took $(($end - $start)) seconds to finish." >> timing_{chunk}.txt\n')
            # Change execution permissions for the generated file
            subprocess.call(['chmod','777',f'{assay_folder}/docking_execution_{chunk}.sh'])
    
        # Update the chunk value
        chunk +=1

    print(f"Succesfully written the docking script/s to: \n \t '{assay_folder}'")
    
def create_docking_conditions_string(conditions_dict):
    """
    Since only non default conditions need to be setted, this function will generate a dict corresponding to the default parameters, compare the provided parameters with the default one, and return a dictionary containing the non-default values to be explicitly indicated in the docking call.
    """    
    
    docking_conditions_string = ""
    
    return  docking_conditions_string
    
def check_docking_assay(registries_db,assay_id):
    # Check if the 'assay_id' exists in the registers database
    conn = tidyscreen.connect_to_db(registries_db)
    cursor = conn.cursor()
    cursor.execute(f"SELECT id FROM docking_registries WHERE id = {assay_id}")
    
    try:
        exists = cursor.fetchone()[0]  # 1 if exists, 0 if not
        pass
    except:
        print(f"Problem retrieving 'assay_id': {assay_id}. Stopping...")
        sys.exit()    
        
def save_df_to_db(db,df,table_name,action="replace"):
    conn = tidyscreen.connect_to_db(db)
    df.to_sql(con=conn, name=table_name,if_exists=action,index=None)
    