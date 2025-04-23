from tidyscreen import tidyscreen as tidyscreen
from tidyscreen.chemspace import cs_utils as cs_utils
import os
import pandas as pd

class ChemSpace:
    def __init__(self, project):
        self.project = project
        self.cs_db_path = self.project.proj_folders_path["chemspace"]["processed_data"]

    def check_cs_database(self):
        if not os.path.exists(f"{self.cs_db_path}/chemspace.db"):
            print("no existe")

    def input_csv(self, file):
        """
        Will read a .csv file and store it into de corresponding database
        """
        target_table_name = file.split("/")[-1].replace(".csv","").replace(".smi","").replace("-", "_") # The last replace will avoid SQL selection actions conflicts
        df = pd.read_csv(file,header=None,index_col=False)
        print(df)
        #df = df.reset_index() # Will add the index as a column to generate an 'id' column
        #print(df)
        first_element = df.iloc[0, 0]  # First row, second column (i.e. the first SMILES)
        # Check if the first element is a valid SMILES
        cs_utils.check_smiles(first_element) # Will stop execution if 'first_element' not a valid SMILES
        
        df = cs_utils.process_input_df(df)
        cs_utils.save_df_to_db(f"{self.cs_db_path}/chemspace.db",df,target_table_name)

        print(f"Table '{target_table_name}' created in: '{self.cs_db_path}/chemspace.db'")

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

    def generate_mols_in_table(self,table_name,charge="bcc",pdb=1,mol2_sybyl=1,mol2_gaff2=1,pdbqt=1):
        """
        Will process all SMILES present in a given table an generate molecules stored in different formats
        """
        db = f"{self.cs_db_path}/chemspace.db"
        cs_utils.process_all_mols_in_table(db,table_name,charge,pdb,mol2_sybyl,mol2_gaff2,pdbqt) # Will generate and store the corresponding .mol2 files in the given table
        # Clean the /tmp directory
        #cs_utils.clean_dir("/tmp")
        
    def retrieve_mols_in_table(self,table_name,outpath=None,ligname=None,pdb=1,mol2_sybyl=1,mol2_gaff2=1,frcmod=1,pdbqt=1):
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
        
        print(f"Ligands extracted to: \n \t '{outpath}")
