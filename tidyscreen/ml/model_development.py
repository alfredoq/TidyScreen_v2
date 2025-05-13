from tidyscreen import tidyscreen as tidyscreen
from tidyscreen.ml import model_development_utils as mdevel_utils
import os
import shutil

class ModelDevelopment:
    def __init__(self, project):
        self.project = project
        self.docking_assays_path = self.project.proj_folders_path["docking"]["docking_assays"]
        self.training_set_path = self.project.proj_folders_path["ml"]["training_sets"]
        self.training_set_db = f"{self.training_set_path}/training_sets.db"
    
    def flag_pose_as_positive(self,assay_id_list,pose_id_list_of_lists):
        
        # Loop over the list of assay_ids
        for index, assay_id in enumerate(assay_id_list):
            # Define the results_db and target_db as a variables
            training_set_db = self.training_set_db
            docking_results_db = f"{self.docking_assays_path}/assay_{assay_id}/assay_{assay_id}.db"        

            # Get the list of poses corresponding to the assay_id
            poses_id_list = pose_id_list_of_lists[index]
    
            ### Process the list of poses and store them into the target db
            mdevel_utils.process_poses_list(assay_id,docking_results_db,training_set_db,poses_id_list,"positives",1)
        
    def flag_pose_as_negative(self,assay_id_list,pose_id_list_of_lists):
        
        # Loop over the list of assay_ids
        for index, assay_id in enumerate(assay_id_list):
            # Define the results_db and target_db as a variables
            training_set_db = self.training_set_db
            docking_results_db = f"{self.docking_assays_path}/assay_{assay_id}/assay_{assay_id}.db"        
        
            # Get the list of poses corresponding to the assay_id
            poses_id_list = pose_id_list_of_lists[index]
        
            ### Process the list of poses and store them into the target db
            mdevel_utils.process_poses_list(assay_id,docking_results_db,training_set_db,poses_id_list,"negatives",0)
    
    def construct_taining_set(self):
        # Define the fingerprints db
        training_set_db = self.training_set_db
        training_set_df, assay_pose_dict = mdevel_utils.combine_fingerprints(training_set_db)
        mdevel_utils.store_training_set(training_set_db, training_set_df,assay_pose_dict)
        
    
    def retrieve_training_set(self,set_id,get_poses=1):
        # Define the fingerprints db
        training_set_db = self.training_set_db
        # Retrieve the training set from the database   
        training_set_df, members_id = mdevel_utils.retrieve_training_set(training_set_db,set_id)
        
        # Save the training set to a CSV file
        output_dir = f"{self.training_set_path}/set_{set_id}"
        mdevel_utils.save_df_to_file(output_dir,training_set_df, set_id)
        
        if get_poses == 1:
            # Retrieve the pdb files from the database
            docking_assays_path = self.docking_assays_path
            mdevel_utils.retrieve_pdb_files(docking_assays_path,output_dir,members_id)
        