from tidyscreen.moldock import moldock_utils as moldock_utils
from tidyscreen import tidyscreen as tidyscreen
import sys
from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run
import io
from contextlib import redirect_stdout
from tidyscreen.chemspace import cs_utils as cs_utils
import os


class MolDock:
    
    def __init__(self, project):
        self.project = project
        self.docking_assays_path = self.project.proj_folders_path["docking"]["docking_assays"]    
        self.docking_registers_path = self.project.proj_folders_path["docking"]["docking_registers"]    
        self.docking_params_path = self.project.proj_folders_path["docking"]["params"]    
        self.receptor_models_path = self.project.proj_folders_path["docking"]["receptors"]
        self.ligands_db = self.project.proj_folders_path["chemspace"]['processed_data'] + "/chemspace.db"
        print("MolDock dimension activated.")
    
    def process_raw_pdb(self, pdb_file, create_dlg_file=True, clean_files=True, x_coord='X', y_coord='Y', z_coord='Z', x_points="X", y_points="Y", z_points="Z"):  
        """
        Will process a raw pdb file using pdb4amber from biobb_amber package. The processed pdb file will be used to generate the receptor mol2 and pdqbt files required for docking with MolDock.
        Args:
            pdb_file (str): Path to the raw pdb file to process
            clean_files (bool): If True, intermediate files generated during the process will be removed. Defaults to True.
        Returns:
            None        
        
        Example:
            >>> la_workshop_moldock.process_raw_pdb("/path/to/raw_receptor.pdb", clean_files=True)  
        
        """
        if os.path.isfile(pdb_file):
            original_pdb_file = pdb_file
            pass
        else:
            print(f"File {pdb_file} does not exist.")
            sys.exit()

        # Create a description file to store with the receptor model
        description = input("Please provide a brief description to store with the receptor model: ") # Will be stored when cleaning the receptor directory

        moldock_utils.store_receptor_description(description, pdb_file)

        # Will check if multiple chainsmoldock_utils. are present in the pdb file
        
        chains = moldock_utils.check_pdb_file_chains(pdb_file)

        # If there are multiple chains, inform the user and ask which chain to keep
        if len(chains) > 1:
            pdb_file = moldock_utils.process_multichain_pdb_file(pdb_file, chains)

        output_file = moldock_utils.check_pdb_file_resnumbers(pdb_file)

        output_info = moldock_utils.process_pdb_with_pdb4amber(pdb_file) 

        non_standard_resids = moldock_utils.get_non_standard_residues(output_info)
        
        if len(non_standard_resids) > 0:
            
            # Query if a non-standard residue is to be kept as part of the receptor
            moldock_utils.manage_non_standard_residues(pdb_file, non_standard_resids, output_file, clean_files)
            
            # Query if any non-standard residue is to be kept as a reference file
            moldock_utils.create_non_standard_ref_file(pdb_file, non_standard_resids, output_file, clean_files)
            
        else:
            
            mol2_file = moldock_utils.prepare_receptor_mol2_only_protein(output_file, clean_files)
            
            moldock_utils.prepare_pdqbt_file(mol2_file)

        if create_dlg_file:
            moldock_utils.create_receptor_dlg_file(pdb_file, x_coord, y_coord, z_coord, x_points, y_points, z_points)

        if clean_files:
            moldock_utils.clean_receptor_dir(original_pdb_file)


    def input_receptor(self,folder):
        """
        Will read a folder in which all files corresponding to a receptor model (including grids) are present. A '.tar' file containing all subfiles will be stored in the corresponding database
        """    
        # Define the database file to store receptor models
        db = f"{self.receptor_models_path}/receptors.db"
        file_prefix = folder.split('/')[-1].replace('-','_') # The last replacement is to avoid conflicts with SQL
        # Will check if the receptor model is already stored by analyzing the folder name 
        ### Function
        moldock_utils.check_existing_rec_model(db,folder)
        # Will generate a .tar file containing all the receptor information
        tar_filename, receptor_blob = moldock_utils.tar_folder(folder,file_prefix)
        # Store the receptor model in the corresponding 'receptors' db
        moldock_utils.store_receptor_model(db,tar_filename, receptor_blob)
    
        print("Successfully stored the receptor model in the database.")

    def tag_obsolete_receptor(self, id_receptor_model):
        """
        Will tag a receptor model as obsolete in the corresponding database
        """
        db = f"{self.receptor_models_path}/receptors.db"
        
        response = input(f"Are you sure you want to tag the receptor model with id: ' {id_receptor_model}' as obsolete? (y/n): ")
        
        moldock_utils.check_description_tag(db,id_receptor_model,tag="#OBSOLETE#")


        if response.lower() == 'y':
            moldock_utils.tag_receptor_model_as_obsolete(db,id_receptor_model)
        else:
            print("Operation cancelled. The receptor model has not been tagged as obsolete.")
            sys.exit()

    def create_docking_params_set(self):
        docking_params_db = f"{self.docking_params_path}/docking_params.db"
        # Will check if the default docking params set has already been created. If exists, returns without further action, otherwise will create the default parameters continue.
        exists = moldock_utils.check_default_docking_conditions(docking_params_db)
        # Create a default params register if it does not exist
        if exists == 0:
            options, description = moldock_utils.create_default_params_register(docking_params_db)
            moldock_utils.store_docking_params_register(docking_params_db,options,description)
            print("Default parameters did not existed, so I created a default set.")
            
        elif exists == 1:
            print("A Default docking condition set already exists. \n The recommended way to create a new docking parameter set is to edit the Default set using SQLiteBrowser. \n Stopping.") 
            sys.exit()   

    def dock_table(self,table_name,id_receptor_model,id_docking_params):

        # Check if the receptor model is valid
        moldock_utils.check_description_tag(f"{self.receptor_models_path}/receptors.db",id_receptor_model,tag="#OBSOLETE#")

        # Will append a docking registry to the corresponding 'db' and return the assay_id
        registries_db = f"{self.docking_registers_path}/docking_registries.db"
        assay_id = moldock_utils.append_docking_registry(registries_db,table_name,id_receptor_model,id_docking_params)
        
        # Create the corresponding docking assay folder containing all required files 
        assay_folder = self.docking_assays_path + f'/assay_{assay_id}'
        moldock_utils.create_assay_folder(assay_folder)
    
        ## Retrieve all ligands from the table to dock
        moldock_utils.retrieve_ligands_for_docking(self.ligands_db,assay_folder,table_name)
    
        ## Retrieve the receptor model to the assay folder
        receptor_db = self.receptor_models_path+'/receptors.db'
        moldock_utils.retrieve_receptor(assay_folder,receptor_db,id_receptor_model)
    
        ## Retrieve docking condictions as dict, compare to default conditions and return non default ones
        params_db = self.docking_params_path+'/docking_params.db'
        conditions_dict = moldock_utils.retrieve_docking_conditions(params_db,id_docking_params)    
    
        custom_parameter_string = moldock_utils.compare_conditions_to_default(params_db,conditions_dict)

        # Create the docking script
        
        moldock_utils.create_docking_executable(assay_folder,custom_parameter_string)