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
# Suppress RDKit warnings and errors
RDLogger.DisableLog('rdApp.*')


class ChemSpace:
    def __init__(self, project):
        self.env_path = '/'.join(tidyscreen.__file__.split('/')[:-1])
        self.project = project
        self.cs_db_path = self.project.proj_folders_path["chemspace"]["processed_data"]
        self.cs_database_file = f"{self.cs_db_path}/chemspace.db"

    def check_cs_database(self):
        if not os.path.exists(f"{self.cs_db_path}/chemspace.db"):
            print("Database does not exists.")

    def input_csv(self, file):
        """
        Will read a .csv file and store it into de corresponding database
        """
        # Read the csv file and check if the first element is a valid SMILES        
        target_table_name, df, first_element = general_functions.csv_reader(file) 
        print(df.columns)
        # Check if the first element is a valid SMILES - Will stop the process if not
        general_functions.check_smiles(first_element)
        # Make all the processing on the generated df (sanitization, enumeration, inchi key calculation)
        df = cs_utils.process_input_df(df,self.cs_database_file,file)
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
        
        print("finished")

    def generate_mols_in_table(self,table_name,charge="bcc",pdb=1,mol2=1,pdbqt=1):
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
                # Evaluate if molecules will fail upon conformer generation - Use timeout function in serial computation.
                print("Evaluating molecules conformational sanity")
                #TODO
                print("Computing .pdb files for ligands")
                # Compute and store the .pdb files using pandarallel
                pandarallel.initialize(progress_bar=True) # Activate Progress Bar
                df.parallel_apply(lambda row: cs_utils.compute_and_store_pdb(row,db,table_name,temp_dir), axis=1)
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
            
