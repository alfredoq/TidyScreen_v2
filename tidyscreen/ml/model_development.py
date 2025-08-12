from tidyscreen import tidyscreen as tidyscreen
from tidyscreen.ml import model_development_utils as mdevel_utils
import os
import shutil
import sys
import glob
from tidyscreen.docking_analysis import docking_analysis_utils as docking_analysis_utils

class ModelDevelopment:
    def __init__(self, project):
        self.project = project
        self.docking_assays_path = self.project.proj_folders_path["docking"]["docking_assays"]
        self.training_set_path = self.project.proj_folders_path["ml"]["training_sets"]
        self.training_set_db = f"{self.training_set_path}/training_sets.db"
    
    def flag_poses_as_positive(self,assay_id_list,pose_id_list_of_lists,comment="manually flagged as positive"):
        
        # Loop over the list of assay_ids
        for index, assay_id in enumerate(assay_id_list):
            # Define the results_db and target_db as a variables
            training_set_db = self.training_set_db
            docking_results_db = f"{self.docking_assays_path}/assay_{assay_id}/assay_{assay_id}.db"        

            # Get the list of poses corresponding to the assay_id
            poses_id_list = pose_id_list_of_lists[index]
    
            # ### Process the list of poses and store them into the target db
            mdevel_utils.process_poses_list(assay_id,docking_results_db,training_set_db,poses_id_list,"positives",1,comment)
        
            ### Store the fingerprint registration
            mdevel_utils.register_fingerprints_addition(training_set_db, assay_id_list,pose_id_list_of_lists,"positives")
        
    def flag_poses_as_negative(self,assay_id_list,pose_id_list_of_lists,comment="manually flagged as negative"):
        
        # Loop over the list of assay_ids
        for index, assay_id in enumerate(assay_id_list):
            # Define the results_db and target_db as a variables
            training_set_db = self.training_set_db
            docking_results_db = f"{self.docking_assays_path}/assay_{assay_id}/assay_{assay_id}.db"        
        
            # Get the list of poses corresponding to the assay_id
            poses_id_list = pose_id_list_of_lists[index]
        
            ### Process the list of poses and store them into the target db
            mdevel_utils.process_poses_list(assay_id,docking_results_db,training_set_db,poses_id_list,"negatives",0,comment)
            
            ### Store the fingerprint registration
            mdevel_utils.register_fingerprints_addition(training_set_db, assay_id_list,pose_id_list_of_lists,"negatives")
    
    def construct_training_set(self, filter_by="all", assay_id=None, date_start=None, date_end=None, comment=None,fp_start=None,fp_end=None,prolif_params_set=None):
        # Define the fingerprints db
        training_set_db = self.training_set_db
        training_set_df, assay_pose_dict = mdevel_utils.combine_fingerprints(training_set_db, filter_by, assay_id, date_start, date_end, comment, fp_start, fp_end, prolif_params_set)
        
        mdevel_utils.store_training_set(training_set_db, training_set_df,assay_pose_dict, filter_by)
    
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
        
        print(f"Training set retrieved and saved successfully to: {output_dir}.")
    
    def visualize_and_flag_docked_poses(self, assay_id, reference_pdb_file, comment="flagged by serial visualization"):
        """
        This function will generate a 3D visualization of the docked pose and allow the user to flag it as positive or negative binding pose.
        """
        # Setup the assay folder path variable
        assay_folder = f"{self.docking_assays_path}/assay_{assay_id}"
        assay_db = self.docking_assays_path + f'/assay_{assay_id}/assay_{assay_id}.db'
        
        # Check if the prolif fingerprints table exists - > return number of registers
        fingerprints_in_db = mdevel_utils.check_fingerprints_results_in_db(assay_db)    
        
        # Check if the docked poses are stored within the folder - > return the number of docked poses
        number_of_docked_poses = mdevel_utils.check_docked_poses(assay_folder)
        
        
        # Check if fingerprints registers match the number of docked poses by a multiplicity 
        if fingerprints_in_db % number_of_docked_poses != 0:
            print(f"Error: The number of fingerprints ({fingerprints_in_db}) does not match the number of docked poses ({number_of_docked_poses}).")
            sys.exit(1)
            return
        
        # Create dataframe containing the subposes and their corresponding Pose_IDs - > return the dataframe
        df = mdevel_utils.retrieve_fingerprints_results(assay_db)
        df_unique_poses = df.drop_duplicates(subset=["sub_pose"])
        
        # Process the dataframe to generate a 3D visualization of the docked pose and create flag lists - > returns the positive and negative lists
        positive_binders_list, negative_binders_list = mdevel_utils.construct_poses_flag_lists(df_unique_poses, reference_pdb_file, assay_folder)
        
        # Store the positive and negative lists in the database
        if len(positive_binders_list) > 0:
            ModelDevelopment.flag_poses_as_positive(self, [assay_id], [positive_binders_list], comment)
        if len(negative_binders_list) > 0:
            ModelDevelopment.flag_poses_as_negative(self, [assay_id], [negative_binders_list], comment)
        
        print("Flagging completed successfully.")
        
    def reprocess_fingerprints_set(self, reference_pdb_file, set_id, vmd_path=None, comment=None):
        
        # In case the comment is not provided, set a default comment
        if comment is None:
            comment = f"reprocessed from set_id {set_id}"
        
        positive_pdb_poses_folder = f"{self.training_set_path}/set_{set_id}/poses/positives"
        negative_pdb_poses_folder = f"{self.training_set_path}/set_{set_id}/poses/negatives"
        
        ### Get the list of pdb files in the folder
        positive_pdb_files = glob.glob(os.path.join(positive_pdb_poses_folder, '*.pdb'))
        negative_pdb_files = glob.glob(os.path.join(negative_pdb_poses_folder, '*.pdb'))
        
        ### Reprocess the stored poses using VMD only the files have not been processed before
        ## Process positive pdb files
        counter = 1
        processed = 0 # Assume that files have not been processed before with VMD for visualization
        for file in positive_pdb_files:
            if counter == 1:
                processed = mdevel_utils.check_pdb_file(file,flag="CRYST")
            
            if processed == 0:
                output_pdb_file = file.replace('.pdb','_vmd.pdb')
                if vmd_path is None:
                    vmd_path = input("Enter the path to VMD executable: ")
                # Report action
                print("Processing .pdb files with VMD for visualization...")
                docking_analysis_utils.process_with_vmd(file,output_pdb_file,vmd_path)
                # Rename the output file to the original one
                os.rename(output_pdb_file, file)

            counter += 1
        
        ## Process negative pdb files
        counter = 1
        processed = 0 # Assume that files have not been processed before with VMD for visualization
        for file in negative_pdb_files:
            if counter == 1:
                processed = mdevel_utils.check_pdb_file(file,flag="CRYST")
            
            if processed == 0:
                output_pdb_file = file.replace('.pdb','_vmd.pdb')
                if vmd_path is None:
                    vmd_path = input("Enter the path to VMD executable: ")
                # Report action
                print("Processing .pdb files with VMD for visualization...")
                docking_analysis_utils.process_with_vmd(file,output_pdb_file,vmd_path)
                # Rename the output file to the original one
                os.rename(output_pdb_file, file)

            counter += 1
        
        ### Visualize and construct the flag lists       
        positive_binders_list, negative_binders_list, possitive_assay_id_list, negative_assay_id_list = mdevel_utils.construct_poses_flag_lists_from_pdb_reprocess_set(reference_pdb_file, positive_pdb_files)
        
        # Store the positive and negative lists in the database
        comment = f"Fingerprints reprocessed from set_id: {set_id}"
        
        for index, pose in enumerate(positive_binders_list):
            assay_id = possitive_assay_id_list[index]
            ModelDevelopment.flag_poses_as_positive(self, assay_id, [[pose]], comment)
            
        for index, pose in enumerate(negative_binders_list):
            assay_id = negative_assay_id_list[index]
            ModelDevelopment.flag_poses_as_negative(self, assay_id, [[pose]], comment)
        
        print("Flagging completed successfully.")