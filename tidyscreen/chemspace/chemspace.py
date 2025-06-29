import warnings
warnings.filterwarnings("ignore")
from tidyscreen import tidyscreen as tidyscreen
from tidyscreen.chemspace import cs_utils as cs_utils
import os
import pandas as pd
from tidyscreen.GeneralFunctions import general_functions as general_functions
from rdkit import RDLogger
from pandarallel import pandarallel
import shutil
import sqlite3
# Suppress RDKit warnings and errors
RDLogger.DisableLog('rdApp.*')

class ChemSpace:
    def __init__(self, project):
        self.env_path = '/'.join(tidyscreen.__file__.split('/')[:-1])
        self.project = project
        self.projects_db = f"{self.env_path}/projects_db/projects_database.db"
        self.cs_db_path = self.project.proj_folders_path["chemspace"]["processed_data"]
        self.cs_database_file = f"{self.cs_db_path}/chemspace.db"

    def check_cs_database(self):
        if not os.path.exists(f"{self.cs_db_path}/chemspace.db"):
            print("Database does not exists.")

    def input_csv(self, file, stereo_enum=0):
        """
        Will read a .csv file and store it into de corresponding database
        """
        # Read the csv file and check if the first element is a valid SMILES        
        target_table_name, df, first_element = general_functions.csv_reader(file) 
        # Check if the first element is a valid SMILES - Will stop the process if not
        general_functions.check_smiles(first_element)
        # Make all the processing on the generated df (sanitization, enumeration, inchi key calculation)
        df = cs_utils.process_input_df(df,self.cs_database_file,file,stereo_enum)
        # Store the final df into de database
        general_functions.save_df_to_db(self.cs_database_file,df,target_table_name)
        
    def list_ligand_tables(self):
        """
        Will list all ligand tables available in the project
        """
        cs_utils.list_ligands_tables(f"{self.cs_db_path}/chemspace.db")

    def delete_table(self,table_name):
        cs_utils.delete_ligands_table(f"{self.cs_db_path}/chemspace.db",table_name)

    def depict_ligand_table(self,table_name):
        db = f"{self.cs_db_path}/chemspace.db"
        #print(self.project.proj_folders_path["chemspace"]["raw_data"])
        output_path = f"{self.project.proj_folders_path['chemspace']['misc']}/{table_name}_depict"
        # Check if the folder is already present
        cs_utils.check_folder_presence(output_path,create=1)
        # call the depiction function
        cs_utils.depict_ligands_table(db,table_name,output_path)
        
        print(f"Successfully depicted ligands in table: '{output_path}'")

    def generate_mols_in_table(self,table_name,charge="bcc-ml",pdb=1,mol2=1,pdbqt=1,conf_rank=0,timeout=10,delete_temp_dir=1):
        """
        Will process all SMILES present in a given table an generate molecules stored in different formats
        """
        # Define the target database (i.e. chemspace)
        db = f"{self.cs_db_path}/chemspace.db"
        # Create the list of columns to create
        list_mol_objects_colnames, list_mol_objects_colnames_types = cs_utils.create_mols_columns_from_selection(pdb,mol2,pdbqt)
        # Connect to the chemspace database
        conn = tidyscreen.connect_to_db(db)
        # Check the existence of the columns in target table. Inform and exist if exists
        cs_utils.check_columns_existence_in_table(conn, table_name, list_mol_objects_colnames)
        # Create the corresponding mols columns in the target table
        cs_utils.create_mols_columns(conn, table_name,list_mol_objects_colnames,list_mol_objects_colnames_types)
        # Read target table into dataframe
        df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
        # Define a temp_dir were all intermeadiate files will be stored
        temp_dir = f"{self.cs_db_path}/temp_dir"
        os.makedirs(temp_dir,exist_ok=True) # Create the corresponding temp directory
        if pdb == 1:
            print("Computing .pdb files for ligands")
            # Compute and store the .pdb files using pandarallel
            pandarallel.initialize(progress_bar=True) # Activate Progress Bar
            df.parallel_apply(lambda row: cs_utils.compute_and_store_pdb(row,db,table_name,temp_dir,conf_rank,timeout), axis=1)
            # Delete all rows in the target table in which .mol2 computation may have failed (errors were registered accondingly)
            general_functions.delete_nulls_table(db,table_name,"pdb_file")
        if mol2 == 1:
            print("Computing .mol2 files for ligands")
            # Get the atom types dictionary for sybyl to gaff2 conversion
            atom_types_dict = f"{self.env_path}/chemspace/sybyl_to_gaff_at_dict.pkl"
            # Compute and store the sybyl mol2 files using pandarallel
            pandarallel.initialize(progress_bar=True) # Activate Progress Bar
            df.parallel_apply(lambda row: cs_utils.compute_and_store_mol2(row,db,table_name,charge,temp_dir,atom_types_dict), axis=1)
            # Purge the rows in which the .mol2 computation may have failed
            general_functions.delete_nulls_table(db,table_name,"mol2_file_sybyl")
            general_functions.delete_nulls_table(db,table_name,"mol2_file_gaff")
        if pdbqt == 1:
            print("Computing .pdbqt files for ligands")
            df.parallel_apply(lambda row: cs_utils.compute_and_store_pdbqt(row,db,table_name,temp_dir), axis=1)
            # Purge the rows in which the .pdbqt computation may have failed
            general_functions.delete_nulls_table(db,table_name,"mol2_file_sybyl")
        
        if delete_temp_dir == 1:
            # Delete the temp directory
            shutil.rmtree(temp_dir, ignore_errors=True)
        
    def retrieve_mols_in_table(self,table_name,outpath=None,ligname=None,pdb=1,mol2_sybyl=1,mol2_gaff2=1,frcmod=1,pdbqt=1,inform=1):
        database_folder = self.project.proj_folders_path["chemspace"]["processed_data"]
        db = f"{database_folder}/chemspace.db"
        
        if outpath == None:
            outpath = f"{self.project.proj_folders_path['chemspace']['misc']}/{table_name}_lig_files"
            os.makedirs(outpath,exist_ok=True)
        
        ## Retrieve each file type
        # pdb files
        if pdb == 1:
            cs_utils.retrieve_blob_ligfiles(db,table_name,outpath,ligname,blob_colname="pdb_file")
        if mol2_sybyl == 1:
            cs_utils.retrieve_blob_ligfiles(db,table_name,outpath,ligname,blob_colname="mol2_file_sybyl")
        if mol2_gaff2 == 1:
            cs_utils.retrieve_blob_ligfiles(db,table_name,outpath,ligname,blob_colname="mol2_file_gaff")
        if frcmod == 1:
            cs_utils.retrieve_blob_ligfiles(db,table_name,outpath,ligname,blob_colname="frcmod_file")
        if pdbqt == 1:
            cs_utils.retrieve_blob_ligfiles(db,table_name,outpath,ligname,blob_colname="pdbqt_file")
        
        if inform == 1:
            print(f"Ligands extracted to: \n \t '{outpath}")

    def subset_chemspace_table(self,source_table,dest_table,colname,filter):
        db = f"{self.project.proj_folders_path['chemspace']['processed_data']}/chemspace.db"
        try:
            general_functions.subset_table(db,source_table,dest_table,colname,filter)
            print(f"Succesfully subseted table: {source_table}")
        except Exception as error:
            print(f"Error substing table {source_table} \n {error}")
            
    def purge_failed_in_table(self,table_name):
        """
        Will purge the table from all rows with null values in any of the mols columns
        """
        db = f"{self.cs_db_path}/chemspace.db"
        # Purge the original table from failed rows
        cs_utils.purge_failed_in_table(db,table_name)
        # Delete the mols columns from the table in order to recompute
        cs_utils.drop_mols_columns(db,table_name)
        
        print(f"Purged table '{table_name}' from failed rows")
        
    def compute_properties(self,table_name,properties_list=["MolWt","MolLogP","NumHDonors","NumHAcceptors","NumRotatableBonds","TPSA"]):
        """
        Will compute the properties for the ligands in the table
        """
        db = f"{self.cs_db_path}/chemspace.db"
        # Read the table into a dataframe
        df = pd.read_sql_query(f"SELECT * FROM {table_name}", tidyscreen.connect_to_db(db))
        # Create the columns in the table to store the properties if they do not exists
        cs_utils.create_properties_columns(db, table_name, properties_list)
        # # Compute the properties using pandarallel
        pandarallel.initialize(progress_bar=True)
        df.parallel_apply(lambda row: cs_utils.compute_properties(row, db, table_name,properties_list), axis=1)
        
    def subset_table_by_properties(self,table_name,props_filter,):
        """
        Will subset the source table by a given property and store the result in the destination table
        """
        db = f"{self.cs_db_path}/chemspace.db"
        try:
            description = input("Enter a description for the subset: ")
            current_subset = cs_utils.write_subset_record_to_db(db,table_name,"by_props",props_filter,description)
            cs_utils.subset_table_by_properties(db,table_name,current_subset,props_filter)
            print(f"Succesfully subseted table: '{table_name}' by properties: '{props_filter}'")
        except Exception as error:
            print(f"Error subseting table {table_name} by property {props_filter} \n {error}")
    
    def add_smarts_filter(self,smarts_filter,description=None):
        try: 
            db = self.projects_db
            cs_utils.check_smarts_existence(db,smarts_filter)
            cs_utils.insert_smarts_filter_in_table(db,smarts_filter,description)
            print(f"Succesfully added SMARTS filter: '{smarts_filter}'")
        except Exception as error:
            print(f"Error inserting SMARTS filter: '{smarts_filter}' \n {error}")
    
    def create_smarts_filters_workflow(self,smarts_filters_dict):
        """
        Will create a workflow to subset a table by SMARTS filters
        """
        db_workflow = f"{self.cs_db_path}/chemspace.db"
        db_filters = self.projects_db
        
        filters_instances_dict, filters_names_list = cs_utils.parse_smarts_filters_dict(db_filters,db_workflow,smarts_filters_dict)
        cs_utils.check_workflow_existence(db_workflow,filters_instances_dict)
        cs_utils.store_smarts_filters_workflow(db_workflow,filters_instances_dict,filters_names_list,smarts_filters_dict)
        print(f"Succesfully created SMARTS filters workflow with filters: '{filters_instances_dict}'")
    
    def subset_table_by_smarts_workflow(self,table_name,workflow_id):
        try:
            db = f"{self.cs_db_path}/chemspace.db"
            # Retrieve the filters instances dict from the database using the workflow_id
            filters_instances_dict = cs_utils.retrieve_workflow_from_db(db,workflow_id)
            # Generate the target table name
            description = input("Enter a description for the subset: ")
            current_subset = cs_utils.write_subset_record_to_db(db,table_name,"by_smarts",filters_instances_dict,description)
            # Filter the table by the SMARTS filters
            filtered_df = cs_utils.subset_table_by_smarts_dict(db,table_name,filters_instances_dict)
            # Store the final df into de database
            general_functions.save_df_to_db(db,filtered_df,current_subset)
            # Inform the user
            print(f"Succesfully subseted table: '{table_name}' by SMARTS filters workflow with ID: '{workflow_id}'")
        except:
            print("An error occurred while trying to subset the table by SMARTS filters workflow. Please check the inputs and try again.")
    
    def list_available_smarts_filters(self):
        """
        Will list all available SMARTS filters in the project
        """
        db = self.projects_db
        cs_utils.list_available_smarts_filters(db)
        
    def copy_table_to_new_name(self, old_table_name, new_table_name):
        """
        Copy the content of an existing table to a new table with a different name.
        """
        db = f"{self.cs_db_path}/chemspace.db"
        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        
        # Create the new table with the same structure as the old table
        cursor.execute(f"CREATE TABLE IF NOT EXISTS {new_table_name} AS SELECT * FROM {old_table_name} WHERE 1=0")
        
        # Copy the content from the old table to the new table
        cursor.execute(f"INSERT INTO {new_table_name} SELECT * FROM {old_table_name}")
        
        conn.commit()
        conn.close()
        
        print(f"Successfully copied table '{old_table_name}' to '{new_table_name}'")
        
    def add_smarts_reaction(self,smarts_reaction,description=None):
        try: 
            db = f"{self.cs_db_path}/chemspace.db"
            cs_utils.check_smarts_reaction_existence(db,smarts_reaction)
            cs_utils.insert_smarts_reaction_in_table(db,smarts_reaction,description)
            print(f"Succesfully added SMARTS reaction: '{smarts_reaction}'")
        except Exception as error:
            print(f"Error inserting SMARTS filter: '{smarts_reaction}' \n {error}")
    
    def add_smarts_reaction_workflow(self,smarts_reactions_id_list):
        """
        Will create a workflow to subset a table by SMARTS reactions
        """
        try: 
            db = f"{self.cs_db_path}/chemspace.db"
            cs_utils.check_smarts_reaction_workflow_existence(db,smarts_reactions_id_list)
            description = input("Enter a description for the SMARTS reactions workflow: ") 
            smarts_reactions_list, smarts_descriptions_list = cs_utils.parse_smarts_reactions_id_list(db,smarts_reactions_id_list)
            cs_utils.store_smarts_reactions_workflow(db,smarts_reactions_list,smarts_reactions_id_list,smarts_descriptions_list,description)
            print(f"Succesfully created SMARTS reactions workflow with reactions: '{smarts_reactions_id_list}'")
            
        except Exception as error:
            print("An error occurred while trying to create the SMARTS reactions workflow. Please check the inputs and try again.")
            print(error)