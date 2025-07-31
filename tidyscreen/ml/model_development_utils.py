from tidyscreen import tidyscreen as tidyscreen
import tarfile
import io
import pandas as pd
import pickle
import sys
import sqlite3
from datetime import datetime
import os
import shutil
import glob
import py3Dmol
import tempfile
import webbrowser
import json

def process_poses_list(assay_id,docking_results_db,training_set_db,pose_id_list,table,flag):
    # Connect to the results database
    conn = tidyscreen.connect_to_db(docking_results_db)
    cursor = conn.cursor()
     
    for pose in pose_id_list:
        # Retrieve the pose information from the docking results database
        ligname, sub_pose, fingerprint_csv_file = retrieve_pose(cursor,pose)
        
        store_fingerprint(training_set_db,assay_id,ligname,sub_pose,fingerprint_csv_file,table,flag,pose)
        
def retrieve_pose(cursor,pose):
    # Retrieve the prolif blob object stored in the results table
    cursor.execute("SELECT LigName, sub_pose, prolif_csv_file FROM prolif_fingerprints WHERE Pose_ID = ?", (pose,))
    result = cursor.fetchall()
    
    if len(result) == 0:
        print(f"No data found for pose ID: {pose}")
        return None, None, None
    
    for ligname, sub_pose, blob in result:
        # Convert BLOB to file-like object
        file_like = io.BytesIO(blob)
        csv_filename = f"/tmp/{sub_pose}.csv"
        
        # Open the tar file
    with tarfile.open(fileobj=file_like, mode='r:*') as tar:
        for member in tar.getmembers():
            # Read the file content
            extracted = tar.extractfile(member)
            if extracted:
                # Save with the name from the database
                with open(f"{csv_filename}", "wb") as out_file:
                    out_file.write(extracted.read())
    
    return ligname, sub_pose, csv_filename
    
def store_fingerprint(training_set_db,assay_id,ligname,sub_pose,fingerprint_csv_file,table,flag,pose):
    # Check if the fingerprint is duplicated in the table
    exists = check_duplicated_fingerprint(training_set_db,assay_id,ligname,sub_pose,table,pose)
    
    if exists == 1:
        print(f"The fingerprint record: '{pose}' corresponding to '{sub_pose}' from assay: '{assay_id}' has already been stored. Skipping storage...")
            
    else:
        # Create a dataframe from the given .csv file
        df = pd.read_csv(fingerprint_csv_file)
        # Add the corresponding flag (0/1 : negative/positive)
        df["binder"] = flag
        df_blob = pickle.dumps(df)
        
        # Connect to the target database and store the fingerprint in the corresponding table
        conn = tidyscreen.connect_to_db(training_set_db)
        cursor = conn.cursor()
        
        cursor.execute(f"INSERT INTO {table} (assay_id, pose_nbr, ligname, sub_pose, fingerprint) VALUES (?,?,?,?,?)", (assay_id,pose,ligname,sub_pose,df_blob,))
        conn.commit()
        conn.close()
    
def check_duplicated_fingerprint(training_set_db,assay_id,ligname,sub_pose,table,pose):
    conn = tidyscreen.connect_to_db(training_set_db)
    cursor = conn.cursor()
    
    # Create the corresponding table if it does not exist
    cursor.execute(f"CREATE TABLE IF NOT EXISTS {table} (assay_id INTEGER, pose_nbr INTEGER, ligname TEXT, sub_pose TEXT, fingerprint BLOB)")
    # Query if the a record matching already exists:
    cursor.execute(f"SELECT * FROM {table} WHERE assay_id = ? AND ligname = ? AND sub_pose = ?", (assay_id,ligname,sub_pose))
    
    # Fetch the result
    result = cursor.fetchone()
    
    # If this is True inform and stop execution
    if result:
        return 1
    else:
        return 0
        
def combine_fingerprints(training_set_db):
    conn = sqlite3.connect(training_set_db)
    cursor = conn.cursor()
    
    # Container for final combined rows
    combined_rows = []
    
    tables = ["positives","negatives"]
    
    # Create a dictionary to store the assay_id and pose_nbr
    assay_pose_dict = {}
    
    for table in tables:
        try: 
        
            # Fetch de data required to construct the dataframe in the POSITIVE table
            cursor.execute(f"SELECT assay_id, pose_nbr, ligname, sub_pose, fingerprint FROM {table}")
            rows = cursor.fetchall()

            # Create a dictionary to store the assay_id and pose_nbr per table
            assay_pose_dict[table] = {}    

            assay_pose_list = []

            for row in rows:
                assay_id, pose_nbr, ligname, sub_pose, fingerprint = row
            
                # Unpickle the blob to get the embedded DataFrame
                df_fingerprint = pickle.load(io.BytesIO(fingerprint))
                
                # Add extra information to the fingerprint dataframe
                df_fingerprint.insert(0, 'sub_pose', sub_pose)
                df_fingerprint.insert(0, 'ligname', ligname)
                df_fingerprint.insert(0, 'assay_id', assay_id)
                
                # Add the df for the current iteration
                combined_rows.append(df_fingerprint)
            
                assay_pose_list.append((assay_id,pose_nbr))

            assay_pose_dict[table] = assay_pose_list
            
            # Concatenate all the rows into one DataFrame
            training_set_df = pd.concat(combined_rows, ignore_index=True)

        except:
            print(f"No table named '{table}' found in the database. Skipping...")
            continue
    
    return training_set_df, assay_pose_dict

def store_training_set(training_set_db, training_set_df,assay_pose_dict):
    
    # Create the table storing datasets if not available
    conn = tidyscreen.connect_to_db(training_set_db) 
    cursor = conn.cursor()
    # Create the projects table only if it does not exists
    cursor.execute('CREATE TABLE IF NOT EXISTS training_datasets (set_id INTEGER PRIMARY KEY AUTOINCREMENT, date DATE, n_positives INTEGER, n_negative INTEGER, members_id BLOB, df_blob BLOB)')
    conn.commit()
    
    # get the time of set preparation
    current_date = datetime.now()
    current_date_formated = current_date.strftime("%Y-%m-%d %H:%M")
    
    # Count positive cases in the dataset
    n_positives = training_set_df[training_set_df["binder"] == 1].shape[0]
    
    # Count negative cases in the dataset
    n_negative = training_set_df[training_set_df["binder"] == 0].shape[0]
    
    # Serialize the dictirionary for storage
    serialized_dict = pickle.dumps(assay_pose_dict)
    
    # Serialize the dataframe for storage
    df_blob = pickle.dumps(training_set_df)
    
    # Create a new entry in the table for the training set
    cursor.execute('INSERT INTO training_datasets (date, n_positives, n_negative, members_id, df_blob) VALUES (?, ?, ?, ?, ?)', (current_date_formated, n_positives, n_negative, serialized_dict, df_blob))
    
    conn.commit()
    conn.close()
    print(f"Training set stored in the database: {training_set_db}")
    
def retrieve_training_set(training_set_db,set_id):
    # Connect to the database
    conn = tidyscreen.connect_to_db(training_set_db)
    cursor = conn.cursor()
    
    # Retrieve the training set from the database
    cursor.execute("SELECT members_id, df_blob FROM training_datasets WHERE set_id = ?", (set_id,))
    result = cursor.fetchone()
    
    if result:
        # Unpickle the blob to get the DataFrame
        members_id_blob = result[0]
        df_blob = result[1]
        members_id = pickle.loads(members_id_blob)
        training_set_df = pickle.loads(df_blob)
        
        return training_set_df, members_id
    else:
        print(f"No training set found with ID: {set_id}")
        return None
    
def save_df_to_file(output_dir,training_set_df, set_id):
    
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    os.makedirs(output_dir)
    
    training_set_df.to_csv(f"{output_dir}/training_set_{set_id}.csv", index=False)
    
def retrieve_pdb_files(docking_assays_path,output_dir,members_id):
    
    
    #for values in members_id.values():
    for key in members_id.keys():
        # Create the target directory to store poses
        target_dir = f"{output_dir}/poses/{key}"
        os.makedirs(target_dir,exist_ok=True)
    
        values = members_id[key]
        for value in values:
            assay_id = value[0]
            pose_id = value[1]
            
            # Connect to the docking assay database
            assay_db = f"{docking_assays_path}/assay_{assay_id}/assay_{assay_id}.db"
            conn = tidyscreen.connect_to_db(assay_db)
            cursor = conn.cursor()
            
            # Retrieve the PDB blob object stored in the results table
            cursor.execute("SELECT sub_pose, docked_pose FROM docked_poses WHERE Pose_ID = ?", (pose_id,))
            result = cursor.fetchall()
            for pose in result:
                ligname, pose_blob = pose
                # Convert BLOB to file-like object
                file_like = io.BytesIO(pose_blob)
                pdb_filename = f"{target_dir}/{ligname}_assay_{assay_id}_Pose_ID_{pose_id}.pdb"
                # Open the tar file
                with tarfile.open(fileobj=file_like, mode='r:*') as tar:
                    for member in tar.getmembers():
                        # Read the file content
                        extracted = tar.extractfile(member)
                        if extracted:
                            # Save with the name from the database
                            with open(f"{pdb_filename}", "wb") as out_file:
                                out_file.write(extracted.read())
                                
def check_fingerprints_results_in_db(db):
    """
    This function will check if the prolif fingerprints table has been computed and stored in the database
    """
    
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    
    try:
        # Select the Pose_IDs from the prolif_fingerprints table
        cursor.execute(f"SELECT Pose_ID FROM prolif_fingerprints ")
        # Count the number of rows
        rows_count = len(cursor.fetchall())
        
        return rows_count
    except Exception as e:
        print(f"Error checking fingerprints in database: {e}. Stopping ...")
        sys.exit(1)

def check_docked_poses(assay_folder):
    """
    This function will check if the docked poses are stored within the folder
    """
    
    folder = f"{assay_folder}/docked_1_per_cluster"
    file_count = len(glob.glob(os.path.join(folder, '*.pdb')))
    
    return file_count

def retrieve_fingerprints_results(db):
    
    conn = tidyscreen.connect_to_db(db)
    # Select the Pose_IDs from the prolif_fingerprints table
    query = "SELECT Pose_ID, sub_pose FROM prolif_fingerprints"
    
    df = pd.read_sql_query(query, conn)
    
    conn.close()
    
    return df

def construct_poses_flag_lists(df, reference_pdb_file, assay_folder):
    
    positive_binders_list = []
    negative_binders_list = []
    
    for index, row in df.iterrows():
        pose_id = row['Pose_ID']
        sub_pose = row['sub_pose']
        
        pose_pdb_file = f"{assay_folder}/docked_1_per_cluster/{sub_pose}.pdb"
        pose_name = f"{sub_pose}.pdb"
        
        
        # Show the 3D visualization and flag the pose
        selected_flag = show_2_molecules_3d_and_flag(reference_pdb_file, pose_pdb_file, pose_name)
  
        # Operate based on the selected flag
        if selected_flag == '+':
            positive_binders_list.append(pose_id)
        elif selected_flag == '-':
            negative_binders_list.append(pose_id)
        elif selected_flag == 's':
            print(f"Skipping pose {pose_id} for sub_pose {sub_pose}.")
            continue
        elif selected_flag == 'f':
            print("Finished flagging poses.")
            break
        
    return positive_binders_list, negative_binders_list
        
def show_2_molecules_3d_and_flag(reference_pdb, molecule_pdb, pose_name):
    
    ref_pdb_str = load_pdb_file(reference_pdb)
    mol_pdb_str = load_pdb_file(molecule_pdb)

    view = py3Dmol.view(width=1000, height=1000)
    view.addModel(ref_pdb_str, 'pdb')
    view.setStyle({'model': 0}, {'stick': {'colorscheme': 'greenCarbon'}})
    view.addModel(mol_pdb_str, 'pdb')
    view.setStyle({'model': 1}, {'stick': {'colorscheme': 'cyanCarbon'}})
    view.addLabel("Your text here", {'position': {'x':0, 'y':0, 'z':0}, 'backgroundColor': 'black', 'fontColor': 'black','fontSize': 50})
    view.zoomTo()
  
    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.html',prefix=f"{pose_name}_") as f:
        html = view._make_html()
        f.write(html)
        temp_html_path = f.name
        
        
        
    webbrowser.open('file://' + os.path.abspath(temp_html_path))
    print(f"Flaging the pose for {pose_name} file {f.name}")
    
    selected_flag = request_flag_from_user()
    
    os.remove(temp_html_path)
    
    return selected_flag

def request_flag_from_user():
    valid_inputs = ['+', '-', 's', 'f']
    while True:
        user_input = input("Flag the pose as '+' or '-' (positive or negative). 's' to skip, 'f' to finish flagging: ")
        if user_input in valid_inputs:
            return user_input
        else:
            print(f"Invalid input. Please enter one of {valid_inputs}.")
    
def load_pdb_file(pdb_path):
    with open(pdb_path, 'r') as f:
        return f.read()
    
def register_fingerprints_addition(db, assays_id_list,poses_id,flag):
    """
    This function will register the addition of a fingerprint to the training set
    """
    
    # Connect to the database
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    
    assays_id_list_str = json.dumps(assays_id_list)
    poses_id_list_str = json.dumps(poses_id)
    
    # Create the fingerprints_additions table if it does not exist
    cursor.execute('CREATE TABLE IF NOT EXISTS fingerprints_additions (register_action INTEGER PRIMARY KEY AUTOINCREMENT, assays_id_list TEXT, poses_id_lists TEXT, flag TEXT, executed_at DATETIME DEFAULT (CURRENT_TIMESTAMP))')
    
    cursor.execute("INSERT INTO fingerprints_additions (assays_id_list, poses_id_lists, flag) VALUES (?,?,?)", (assays_id_list_str,poses_id_list_str,flag,))
    conn.commit()
    conn.close()