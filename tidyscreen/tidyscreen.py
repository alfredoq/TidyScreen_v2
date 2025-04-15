import sys
import os
import sqlite3
from pathlib import Path
import shutil
import pickle

# This module is dedicated to interact with the general actions on the TidyScreen installed locally as an environment 
main_db = sys.modules['tidyscreen'].__file__.replace('__init__.py','projects_db/projects_database.db')

def check_db(print_output=1):
    # Will check the existence of the database that should be installed when TidyScreen is installed. 

    if os.path.exists(main_db):
        if print_output == 1:
            print("TidyScreen main database found! Continuing...")
    else:
        print("TidyScreen main database does NOT exist. Creating it...")
        create_projects_table()

def create_projects_table():
    conn = connect_to_db(main_db) # Connect to the main database
    cur = conn.cursor()
    # Create the projects table only if it does not exists
    cur.execute('CREATE TABLE IF NOT EXISTS projects (name VARCHAR, path VARCHAR, description VARCHAR)')
    conn.commit()

def connect_to_db(db_file):
    conn = sqlite3.connect(f'{db_file}')
    return conn

def projects(print_output=1):
    # Will list available projects in the database.
    check_db(print_output) # Verify the existence of the database
    create_projects_table() # Create the projects table in case it is not present in the .db
    conn = connect_to_db(main_db) # Connect to the main database
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM projects")
    rows = cursor.fetchall()

    if print_output == 1: # Check if the intention is to print to screen (print_output == 1) or returnt 'rows' variable (print_output == 0)
        if len(rows) == 0:
            print("No projects available in database")
        else:
            for row in rows:
                print(f'Project: {row[0]} \n \t located at {row[1]}')
    else:
        return rows

def check_project_exists(rows,proj_name):
    # Will check the existence of a given project
    exists = any(tup[0] == proj_name for tup in rows)

    if exists:
        #print(f'Project {proj_name} exists.')
        return 1
    else:
        #print(f'Project {proj_name} does NOT exists.')
        return 0

def get_project_local_path(proj_name):
    # Get all available projects in the main database
    rows = projects(print_output=0)
    # Iterate the rows in whithin the projects list
    for row in rows:
        if row[0] == proj_name:
            return row[1] # Return the project path

def create_project(path, proj_name):
    ### Will create a project and associated file structure
    rows = projects(print_output = 0)
    # Check the existence of the project to delete
    exists = check_project_exists(rows,proj_name)
    if exists == 0:
        # Create the corresponding project
        description = input("Input project short description: ")
        proj_path = f"{path}/{proj_name}"
        all_paths_dict = create_project_folder(proj_path)
        store_project_in_main_database(proj_name,proj_path,all_paths_dict,description)
        print(f"Project '{proj_name}' created at: '{proj_path}'.")
    else:
        print(f"Project '{proj_name}' already exists. Cannot be created. Stopping...")
        sys.exit()

def del_project(proj_name):
    ### Will delete a single project and all associated files
    rows = projects(print_output = 0)
    # Check the existence of the project to delete
    exists = check_project_exists(rows,proj_name)
    if exists == 1:
        # Delete the corresponding project and associated files
        delete_project(proj_name)
    else:
        print(f"The project '{proj_name}' does NOT exist. Stopping...")

def create_project_folder(path):
    # This will create the project structure at the specified path
    # Define the folder structure of a new project
    project_folders_structure = {'chemspace':["raw_data","processed_data","misc"],
                                'docking':["docking_assays","docking_registers","params",'raw_data','receptors'],
                                'dynamics':["md_assays","md_registers","md_params"],
                                '.project_vars':["paths"]
                                }

    # Initialize a dictionary to contain all project paths
    all_paths_dict = {}

    for base_folder in project_folders_structure.keys():
        current_level_dict = {}
        for folder in project_folders_structure[base_folder]:
            new_folder = f"{path}/{base_folder}/{folder}"
            # Create the corresponding folder
            Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
            # Create a dictionary with folder for the current level of folder
            current_level_dict[folder]=new_folder
            # Create a dictionary for the base level of folders
            all_paths_dict[base_folder]=current_level_dict

    # Save the project paths to a pickle file
    project_paths_dir = f"{path}/.project_vars/paths"
    with open(f"{project_paths_dir}/paths.pkl","wb") as paths_file:
        pickle.dump(all_paths_dict, paths_file)

def store_project_in_main_database(name,path,all_paths_dict,description):
    conn = connect_to_db(main_db) # Connect to the main database
    cur = conn.cursor()
    # Create the projects table only if it does not exists
    cur.execute('CREATE TABLE IF NOT EXISTS projects (name VARCHAR, path VARCHAR, description VARCHAR)')
    # Store the project register
    cur.execute(f'INSERT INTO projects (name, path, description) values (?,?,?)', (name, path,description))
    conn.commit()

def delete_project(name):
    # Provide a confimartion prior to deletion
    confirm = input("Do you confirm '{name}' deletion? (y/n): ")
    
    if confirm == 'y':
    # Check if the project exists
        conn = connect_to_db(main_db) # Connect to the main database
        delete_project_folder(conn,name)
        delete_db_registry(conn,name)
    else:
        print(f"Aborting deleteion of project: '{name}'")
    
def delete_db_registry(conn,name):
    # Get the ID of the project to delete
    cur = conn.cursor()
    cur.execute(f"DELETE FROM projects WHERE name = '{name}'")
    conn.commit()
    print(f"Project '{name}' deleted from main database.")

def delete_project_folder(conn,name):
    # Get the path of the project to delete
    cur = conn.cursor()
    cur.execute(f"SELECT path FROM projects WHERE name = '{name}'")
    project_path = cur.fetchall()[0][0]
    shutil.rmtree(project_path)
    print(f"Project folder '{project_path}' deleted filesystem.")

class ActivateProject:

    def __init__(self, projname):
        # Check if the called project exists
        rows = projects(print_output = 0)
        # Check the existence of the project to delete
        exists = check_project_exists(rows,projname)
        
        if exists == 1:
            self.projname = projname
            self.base_proj_path = get_project_local_path(projname)
            self.maindbpath = main_db
            self.proj_folders_path = ActivateProject.get_proj_paths(self)

    def get_proj_paths(self):
        with open(f"{self.base_proj_path}/.project_vars/paths/paths.pkl", "rb") as paths_file:
            my_dict = pickle.load(paths_file)
            return my_dict




class Engine:
    def __init__(self, horsepower):
        self.horsepower = horsepower

    def start(self):
        return "Engine started!"
