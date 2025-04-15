from tidyscreen.moldock import moldock_utils as moldock_utils
from tidyscreen import tidyscreen as tidyscreen
import sys

class MolDock:
    
    def __init__(self, project):
        self.project = project
        self.docking_assays_path = self.project.proj_folders_path["docking"]["docking_assays"]    
        self.docking_registers_path = self.project.proj_folders_path["docking"]["docking_registers"]    
        self.docking_params_path = self.project.proj_folders_path["docking"]["params"]    
        self.receptor_models_path = self.project.proj_folders_path["docking"]["receptors"]
        self.ligands_db = self.project.proj_folders_path["chemspace"]['processed_data'] + "/chemspace.db"
    
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
    
    def create_docking_params_set(self):
        docking_params_db = f"{self.docking_params_path}/docking_params.db"
        # Will check if the default docking params set has already been created. If exists, returns without further action, otherwise will create the default parameters continue.
        exists = moldock_utils.check_default_docking_conditions(docking_params_db)
        # Create a default params register if it does not exist
        if exists == 0:
            options, description = moldock_utils.create_default_params_register(docking_params_db)
            moldock_utils.store_docking_params_register(docking_params_db,options,description)
            
        elif exists == 1:
            print("A Default docking condition set already exists. \n The recommended way to create a new docking parameter set is to edit the Default set using SQLiteBrowser. \n Stopping.") 
            sys.exit()   

    def dock_table(self,table_name,id_receptor_model,id_docking_params):
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
        
