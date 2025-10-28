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
import sqlite3
from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run
import io
from contextlib import redirect_stdout
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
import shutil

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

def tag_receptor_model_as_obsolete(db,id_receptor_model):

    conn = sqlite3.connect(db)
    cursor = conn.cursor()  

    # Get the current string 
    cursor.execute(f"SELECT description FROM receptors WHERE id = {id_receptor_model}")
    current_string = cursor.fetchone()[0]

    # Concatenate the "#OBSOLETE#"" tag to the current string
    new_string = current_string + " #OBSOLETE#"

    # Update the description in the database
    cursor.execute(f"UPDATE receptors SET description = ? WHERE id = ?", (new_string, id_receptor_model))
    conn.commit()
    conn.close()    

    print(f"Receptor model with id: {id_receptor_model} has been tagged as obsolete.")  

def check_description_tag(db,id_receptor_model,tag):
    """
    Will check if the receptor model with id: 'id_receptor_model' has the tag: 'tag' in its description.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()  

    # Get the current string 
    cursor.execute(f"SELECT description FROM receptors WHERE id = ?", (id_receptor_model,))
    current_string = cursor.fetchone()[0]

    if tag in current_string:
        print(f"The receptor model with id: {id_receptor_model} has tag: '{tag}' in its description. Exiting...")
        conn.close()
        sys.exit()
    
def process_pdb_with_pdb4amber(pdb_file):

    receptor_output_temp_file = pdb_file.split('/')[-1].replace('.pdb','_temp.pdb')
    receptor_output_file = pdb_file.split('/')[-1].replace('.pdb','_processed.pdb')
    
    receptor_output_temp_file = '/'.join(pdb_file.split('/')[:-1]) + '/' + receptor_output_temp_file
    receptor_output_file = '/'.join(pdb_file.split('/')[:-1]) + '/' + receptor_output_file

    prop = {'remove_tmp': True,
                }

    f = io.StringIO()
    
    with redirect_stdout(f):
        pdb4amber_run(input_pdb_path=pdb_file, output_pdb_path=receptor_output_temp_file, properties=prop)

        output = f.getvalue().splitlines()

    # Delete the temp file
    os.remove(receptor_output_temp_file)

    # Remove log files
    os.remove("log.err")
    os.remove("log.out")

    return output

def get_non_standard_residues(strings):

    for idx, val in enumerate(strings):
            if "Non-standard-resnames" in val and idx + 1 < len(strings):
                non_standard_resids_temp = strings[idx + 1].split(',')
                non_standard_resids = [item.replace(" ", "") for item in non_standard_resids_temp]

    if len(non_standard_resids) == 0:
        return None
    else:
        print("Non-standard residues detected:")
        return non_standard_resids

def reinsert_non_standard_residue(receptor_output_temp_file, receptor_output_file, residue_to_mantain):
    
    # Get matching lines from file1
    with open(receptor_output_temp_file, "r") as f1:
        matching_lines = [line for line in f1 if residue_to_mantain in line]

    print(matching_lines)

    receptor_file = receptor_output_file.split('/')[-1]
    ligand_filename = receptor_output_file.replace(receptor_file,'ligand.pdb')
    
    print(ligand_filename)
    
    with open(ligand_filename, "a") as ligand_file:
        ligand_file.writelines(matching_lines)

    # Append them to file2
    with open(receptor_output_file, "a") as f2:
        f2.write("TER\n")
        f2.writelines(matching_lines)
        f2.write("TER\n")

    print(receptor_output_file)

    return ligand_filename

def add_hydrogens_to_ligand(ligand_filename, residue_to_mantain):
    
    ligand_filename_hs = ligand_filename.replace('ligand.pdb','ligand_hs.pdb')
    ligand_filename_hs_mol2 = ligand_filename.replace('ligand.pdb','ligand_hs.mol2')
    
    mol = Chem.MolFromPDBFile(ligand_filename, removeHs=False)
    mol_with_h = Chem.AddHs(mol)
    
    # Generate 3D coordinates (if not present or to optimize geometry)
    AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
    #AllChem.UFFOptimizeMolecule(mol_with_h)
    
    Chem.MolToPDBFile(mol_with_h, ligand_filename_hs)
    rdmolfiles.MolToMolFile(mol_with_h, ligand_filename_hs_mol2)
    
    
    # Replace 'UNL' label by the residue name to mantain
    with open(ligand_filename_hs, 'r') as file:
        filedata = file.read()
        filedata = filedata.replace('UNL', residue_to_mantain)
    
    with open(ligand_filename_hs, 'w') as file:
        file.write(filedata)    
    
    # Delete CONECT lines
    with open(ligand_filename_hs, 'r') as file:
        lines = file.readlines()
    with open(ligand_filename_hs, 'w') as file:
        for line in lines:
            if not line.startswith('CONECT'):
                file.write(line)    
    
    ligand_filename_hs_fixed = ligand_filename.replace('ligand.pdb','ligand_hs_fixed.pdb')
    new_resnum = 999
    with open(ligand_filename_hs, "r") as infile, open(ligand_filename_hs_fixed, "w") as outfile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                # Replace residue number (columns 23-26, 1-based indexing)
                new_line = line[:22] + f"{new_resnum:4d}" + line[26:]
                outfile.write(new_line)
            else:
                outfile.write(line)
    
    return ligand_filename_hs_fixed

def check_pdb_file_resnumbers(pdb_file):
    """
    Will check if the residue numbers in the pdb file are sequential. If not, will raise a ValueError and stop the process.
    """
    
    output_tempfile = pdb_file.replace('.pdb','_temp.pdb')
    output_file = pdb_file.replace('.pdb','_checked.pdb')


    # Create a file corresponding to the protein
    with open(pdb_file, "r") as infile, open(output_tempfile, "w") as outfile:
        for line in infile:
            if "ATOM" in line:
                    outfile.write(line)

    resnums = []

    with open(output_tempfile, "r") as f1, open(output_file, "w") as f2:
        for line in f1:
            if line.startswith(("ATOM", "HETATM")):
                resnum = int(line[22:26].strip())
                if resnum not in resnums:
                    if len(resnums) > 0:
                        last_resnum = resnums[-1]
                        if last_resnum + 1 != resnum:
                            print(f"gap between: {last_resnum} and {resnum}. 'TER' card inserted")
                            f2.write("TER\n")

                    resnums.append(resnum)
                f2.write(line)

    # Append a TER card to the processed protein
    with open(output_file, "a") as f:
        f.write("TER\n")    

    # Delete the temporary file
    os.remove(output_tempfile)

    return output_file

def prepare_receptor_mol2_only_protein(pdb_file, clean_files):
    
    file_path = pdb_file.rsplit('/', 1)[0]
    
    files_to_remove = []
    
    with open(f"{file_path}/tleap1.in", "w") as f:
        f.write(f"source leaprc.protein.ff14SB\n")
        f.write(f"receptor = loadpdb {pdb_file}\n")
        f.write(f"savepdb receptor {file_path}/receptor.pdb\n")
        f.write(f"saveamberparm receptor {file_path}/receptor.prmtop {file_path}/receptor.inpcrd\n")
        f.write(f"quit\n")
    
    f.close()
            
    tleap_path = shutil.which('tleap')
    
    # Compute mol2 with Sybyl Atom Types - for compatibility with RDKit and Meeko
    tleap_command1 = f'cd {file_path} && {tleap_path} {file_path} -f {file_path}/tleap1.in'
    
    # Execute the command
    subprocess.run(tleap_command1, shell=True, capture_output=True, text=True)
    
    # Check the corresponding output files have been created
    target_files = [f"{file_path}/receptor.prmtop", f"{file_path}/receptor.inpcrd"]
    
    for file in target_files:
        if os.path.getsize(file_path) > 0:
            files_to_remove.append(f"{file_path}/tleap1.in")
            files_to_remove.append(f"{file_path}/receptor.prmtop")
            files_to_remove.append(f"{file_path}/receptor.inpcrd")
            files_to_remove.append(f"{file_path}/leap.log")
        else:
            print(f"File {file} was not created. Stopping.")
            sys.exit()
    
    # Check and inform if alternate residue locations are present
    
    with open(f'{file_path}/leap.log', 'r') as f:
        content = f.read()
        search_string = "Atom names in each residue should be unique."
        if search_string in content:
            print("ATTENTION: check alternate residue conformations. Mantaining the first one.")
        else:
            pass
            
    ### Prepare an Amber-like .mol2 file
    with open(f"{file_path}/cpptraj.in", "w") as f:
        f.write(f"parm {file_path}/receptor.prmtop\n")
        f.write(f"trajin {file_path}/receptor.inpcrd\n")
        f.write(f"trajout {file_path}/receptor_temp.mol2\n")
        f.write(f"go\n")
        f.write(f"quit\n")
    
    f.close()
    
    cpptraj_path = shutil.which('cpptraj')
    
    cpptraj_command = f'cd {file_path} && {cpptraj_path} -i {file_path}/cpptraj.in'
    
    # Check the corresponding output files have been created
    target_files = [f"{file_path}/receptor_temp.mol2"]
    
    for file in target_files:
        if os.path.getsize(file_path) > 0:
            files_to_remove.append(f"{file_path}/cpptraj.in")
            files_to_remove.append(f"{file_path}/receptor_temp.mol2")
        else:
            print(f"File {file} was not created. Stopping.")
            sys.exit()
    
    
    # Execute the command
    subprocess.run(cpptraj_command, shell=True, capture_output=True, text=True)
    
    ### Prepare an autodock-like .mol2 file
    with open(f"{file_path}/tleap2.in", "w") as f:
        f.write(f"source leaprc.protein.ff14SB\n")
        f.write(f"REC = loadmol2 {file_path}/receptor_temp.mol2\n")
        f.write(f"savemol2 REC {file_path}/receptor.mol2 0\n")
        f.write(f"quit\n")
    
    f.close()
    
    tleap_command2 = f'cd {file_path} && {tleap_path} {file_path} -f {file_path}/tleap2.in'
    
    subprocess.run(tleap_command2, shell=True, capture_output=True, text=True)
    
    # Check the corresponding output files have been created
    target_files = [f"{file_path}/receptor.mol2"]
    
    for file in target_files:
        if os.path.getsize(file_path) > 0:
            print(f"File {file} created successfully.")
            files_to_remove.append(f"{file_path}/tleap2.in")
        else:
            print(f"File {file} was not created. Stopping.")
            sys.exit()
    
    # Remove the files added to the 'files_to_remove' list  
    if clean_files:
        files_to_remove_unique = list(set(files_to_remove))
        for file in files_to_remove_unique:
            os.remove(file)
        
    return f"{file_path}/receptor.mol2"

def prepare_pdqbt_file(mol2_file):
    
    file_path = mol2_file.rsplit('/', 1)[0]
    output_file = mol2_file.replace('.mol2','.pdbqt')
    
    prepare_receptor_path = shutil.which('prepare_receptor4.py')
    
    # Compute mol2 with Sybyl Atom Types - for compatibility with RDKit and Meeko
    prepare_receptor_command = f'cd {file_path} && {prepare_receptor_path} -r {mol2_file} -C'
    
    # Execute the command
    subprocess.run(prepare_receptor_command, shell=True, capture_output=True, text=True)
    
    if os.path.getsize(output_file) > 0:
        print(f"File {output_file} created successfully.")
    else:
        print(f"File {output_file} was not created. Stopping.")
        sys.exit()
        
def check_pdb_file_chains(pdb_file):
    """
    Will return a set with the chain identifiers present in the pdb file.
    Args:
        pdb_file (str): Path to the PDB file.
    Returns:
        set: A set containing the unique chain identifiers.
    
    Example:
        chains = check_pdb_file_chains("example.pdb")
        print(chains)  # Output: {'A', 'B', 'C'}    
    
    """
    chains = set()
    
    if not os.path.exists(pdb_file):
        print(f"Error: File {pdb_file} does not exist.")
        return None
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    chain_id = line[21].strip()
                    if chain_id:  # Skip empty chain IDs
                        chains.add(chain_id)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    
    return chains

def process_multichain_pdb_file(pdb_file, chains):
    print(f"Chains detected: {chains}")    
    chain_to_keep = input("Provide the chain identifier to mantain: ('all' to keep all chains') ")
    
    # Check if the user wants to mantain all chains
    if chain_to_keep.lower() == 'all':
        print("All chains will be mantained.")
        return pdb_file
    
    # Check if the provided chain is valid
    if chain_to_keep not in chains:
        print(f"Error: Chain '{chain_to_keep}' not found in the PDB file. Stopping.")
        sys.exit()  
        
    output_file = pdb_file.replace('.pdb',f'_chain_{chain_to_keep}.pdb')
    
    with open(pdb_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith(('ATOM', 'HETATM')):
                current_chain = line[21].strip()
                if current_chain == chain_to_keep:
                    outfile.write(line)
    
    return output_file

def clean_receptor_dir(pdb_file):
    
    pdb_directory = '/'.join(pdb_file.split('/')[:-1])
    
    temp_files = glob(f"{pdb_directory}/*")
    
    files_to_retain = [f"{pdb_directory}/receptor.mol2",
                       f"{pdb_directory}/receptor.pdb",
                       f"{pdb_directory}/receptor.pdbqt",
                       f"{pdb_directory}/reference_ligand.pdb",
                       f"{pdb_directory}/receptor.gpf",
                       pdb_file,
                       ]
    
    # Delete the files that should not be retained
    for file in temp_files:
        if file not in files_to_retain:
            os.remove(file)

def manage_non_standard_residues(pdb_file, non_standard_resids, output_file, clean_files):
    
    mantain_non_standard = input(f"The following non-standard residues were found in the pdb file: {non_standard_resids}. Do you want to mantain ONE of them in the processed receptor pdb file? (y/n): ")
            
    if mantain_non_standard.lower() == 'y':
        residue_to_mantain = input("Type the 3-letter code of the residue you want to mantain: ")
                
        if residue_to_mantain in non_standard_resids:
            ligand_filename = reinsert_non_standard_residue(pdb_file, output_file, residue_to_mantain)
            print(f"The non-standard residue {residue_to_mantain} has been reinserted in the processed pdb file.")
        else:

            print(f"The residue {residue_to_mantain} is not in the list of non-standard residues found. No residues were reinserted.")

    else:
        print("Only the protein will be kept in the processed pdb file.")    
                
        mol2_file = prepare_receptor_mol2_only_protein(output_file, clean_files)
                
        prepare_pdqbt_file(mol2_file)
    
def create_non_standard_ref_file(pdb_file, non_standard_resids, output_file, clean_files):
    
    mantain_non_standard = input(f"The following non-standard residues were found in the pdb file: {non_standard_resids}. Do you want to mantain ONE of them as a REFERENCE pdb file? (y/n): ")
    
    if mantain_non_standard.lower() == 'y':
        residue_to_mantain = input("Type the 3-letter code of the residue you want to save as REFERENCE FILE: ")
        
        if residue_to_mantain in non_standard_resids:
            ligand_filename = save_non_standard_residue_ref_file(pdb_file, residue_to_mantain)
            
def save_non_standard_residue_ref_file(pdb_file, residue_to_mantain):
    
    pdb_directory = '/'.join(pdb_file.split('/')[:-1])
    
    ligand_filename = f'{pdb_directory}/reference_ligand.pdb'
    
    print(ligand_filename)
    
    # Get matching lines from from the pdb file 
    with open(pdb_file, "r") as f1:
        matching_lines = [line for line in f1 if residue_to_mantain in line]

    with open(ligand_filename, "w") as ligand_file:
        ligand_file.writelines(matching_lines)

    print(f"Reference file for non-standard residue {residue_to_mantain} saved as: {ligand_filename}")

def store_receptor_description(description, pdb_file):
    
    pdb_directory = '/'.join(pdb_file.split('/')[:-1])
    
    with open(f"{pdb_directory}/receptor_description.txt", "w") as f:
        f.write(description)    
    f.close()
    
def create_receptor_dlg_file(pdb_file, x_coord, y_coord, z_coord, x_points, y_points, z_points):
    
    pdb_directory = '/'.join(pdb_file.split('/')[:-1])
    gpf_file = f"{pdb_directory}/receptor.gpf"
    
    with open(gpf_file, "w") as f:
        f.write(f"npts {x_points} {y_points} {z_points} # num.grid points in xyz \n")
        f.write("gridfld receptor.maps.fld # grid_data_file \n")
        f.write("spacing 0.375  # spacing(A) \n")
        f.write("receptor_types A C NA OA N SA HD  # receptor atom types \n")
        f.write("ligand_types H A C NA OA N HD S SA F Cl Br # ligand atom types \n")
        f.write("receptor receptor.pdbqt # macromolecule \n")
        f.write(f"gridcenter {x_coord} {y_coord} {z_coord}  # xyz-coordinates or auto \n")
        f.write("smooth 0.5  # store minimum energy w/in rad(A) \n")
        f.write("map receptor.H.map # atom-specific affinity map \n")
        f.write("map receptor.A.map # atom-specific affinity map \n")
        f.write("map receptor.C.map # atom-specific affinity map \n")
        f.write("map receptor.NA.map # atom-specific affinity map \n")
        f.write("map receptor.OA.map # atom-specific affinity map \n")
        f.write("map receptor.N.map # atom-specific affinity map \n")
        f.write("map receptor.HD.map # atom-specific affinity map \n")
        f.write("map receptor.S.map # atom-specific affinity map \n")
        f.write("map receptor.SA.map # atom-specific affinity map \n")
        f.write("map receptor.F.map # atom-specific affinity map \n")
        f.write("map receptor.Cl.map # atom-specific affinity map \n")
        f.write("map receptor.Br.map # atom-specific affinity map \n")
        f.write("elecmap receptor.e.map # electrostatic potential map \n")
        f.write("dsolvmap receptor.d.map  # desolvation potential map \n")
        f.write("dielectric -0.1465  # <0, AD4 distance-dep.diel;>0, constant \n")
    
def check_non_standard_residues_in_pdb(resid):
    
    protein_residues = ['ACE','ALA','ARG','ASH','ASN','ASP','CALA','CARG','CASN','CASP','CCYS','CCYX','CGLN','CGLU','CGLY','CHID','CHIE','CHIP','CHIS','CHYP','CILE','CLEU','CLYS','CMET','CPHE','CPRO','CSER','CTHR','CTRP','CTYR','CVAL','CYM','CYS','CYX','GLH','GLN','GLU','GLY','HID','HIE','HIP','HIS','HYP','ILE','LEU','LYN','LYS','MET','NALA','NARG','NASN','NASP','NCYS','NCYX','NGLN','NGLU','NGLY','NHE','NHID','NHIE','NHIP','NHIS','NILE','NLEU','NLYS','NME','NMET','NPHE','NPRO','NSER','NTHR','NTRP','NTYR','NVAL','PHE','PRO','SER','THR','TRP','TYR','VAL']
    
    if resid in protein_residues:
        return 0
    else:
        return 1
            
            