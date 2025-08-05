import warnings
warnings.filterwarnings("ignore")
from tidyscreen.docking_analysis import docking_analysis_utils as docking_analysis_utils
from tidyscreen.moldyn import moldyn_utils as md_utils
from tidyscreen.GeneralFunctions import general_functions as general_functions
import time
import shutil
import sys
import json
import py3Dmol
import ast

class DockingAnalysis:
    
    def __init__(self, project,amberhome=None):
        self.project = project
        self.docking_assays_path = self.project.proj_folders_path["docking"]["docking_assays"]
        self.docking_registers_path = self.project.proj_folders_path["docking"]["docking_registers"]    
        self.docking_params_path = self.project.proj_folders_path["docking"]["params"]    
        self.receptor_models_path = self.project.proj_folders_path["docking"]["receptors"]
        self.ligands_db = self.project.proj_folders_path["chemspace"]['processed_data'] + "/chemspace.db"

        if amberhome is None:
            self.amberhome = input("Please, input the AMBERHOME path: ")
            if self.amberhome is None:
                print("Error: AMBERHOME environment variable is not set.")
                sys.exit()
        else:
            self.amberhome = amberhome
    
    def process_docking_assay(self, assay_id, max_poses=10, vmd_path=None, extract_poses=0):
        registries_db = f"{self.docking_registers_path}/docking_registries.db"
        # Check if the 'assay_id' existis in docking registries
        docking_analysis_utils.check_docking_assay(registries_db,assay_id)
        # Process the corresponding 'dlg' files
        assay_folder = self.docking_assays_path + f'/assay_{assay_id}'
        results_db_file = docking_analysis_utils.process_dlg_files(assay_folder,assay_id,max_poses)
        # Add docking subposes in database
        docking_analysis_utils.add_docking_subposes_nbr(results_db_file)
        if extract_poses == 1:
            # Extract 1 pdb pose per cluster
            print(f"\n Extracting '{max_poses}' PDB poses per cluster \n")
            docking_analysis_utils.extract_1_pdb_per_cluster(assay_folder,results_db_file,max_poses,vmd_path)
    
    def compute_fingerprints_for_docked_pose(self, assay_id, results_pose_id, mmgbsa=1, prolif=1, clean_files=1, clean_folder=1, solvent="implicit", min_steps=5000, store_docked_poses=1, prolif_parameters_set=1, iteration=1,  ligresname="UNL"):
    
    ### Start to log the time
        start_time = time.time()
    ### Create a custom folder for the analysis and copy/generate relevant files
        
        # Setup the assay folder path variable
        assay_folder = self.docking_assays_path + f'/assay_{assay_id}'
        
        # Set the main fingerprints folder
        main_fingerprints_folder = f"{assay_folder}/fingerprints_analyses"
        
        # Create the fingerprints analysis folder
        complex_pdb_file, output_path, receptor_filename, ligname, sub_pose, pose_pdb_file  = docking_analysis_utils.create_fingerprints_analysis_folder(self,assay_folder,assay_id, results_pose_id)
        
        if iteration == 1:
            # Create a tleap vs cristal reference dictionary
            tleap_vs_cristal_reference_dict = docking_analysis_utils.compute_tleap_vs_cristal_reference_dict(receptor_filename, f"{assay_folder}/fingerprints_analyses")
            
        else:
            # Load the existing tleap vs cristal reference dictionary
            with open(f"{main_fingerprints_folder}/tleap_vs_cristal_resnames_dict.csv", 'r') as f:
                tleap_vs_cristal_reference_dict = json.load(f)
        
        # Inform ligand under processing:
        print(f"\n Processing ligand: {ligname} - Pose: {int(results_pose_id):.0f} \n")
    
        if mmgbsa == 1:
        
        ### Compute fingerprints using MMPBSA on target folder
            
            # Set the AMBERHOME environment variable
            amberhome = self.amberhome
            
            # Compute the MMPBSA based per-residue interaction fingerprint
            
            try:
                ## Compute the MMGBSA fingerprints
                prmtop_file, crd_file, tleap_vs_cristal_reference_dict,mmpbsa_decomp_csv_output = docking_analysis_utils.compute_fingerprints(output_path,main_fingerprints_folder,complex_pdb_file,receptor_filename,solvent,min_steps,iteration,ligresname,amberhome)
        
                # Store the MMGBSA results in the database
                docking_analysis_utils.store_mmbgsa_fingerprints_results_in_db(assay_folder,assay_id,results_pose_id,ligname,sub_pose,complex_pdb_file,mmpbsa_decomp_csv_output)
        
            except Exception as e:
                print(f"Error during MMGBSA fingerprint computation: {e}")
                sys.exit()
        
        
        if prolif == 1:
            
        ### Compute ProLIFfingerprints for minimized pose ###
            interactions_list = ['Anionic','CationPi','Cationic','EdgeToFace','FaceToFace','HBAcceptor','HBDonor','Hydrophobic','MetalAcceptor','MetalDonor','PiCation','PiStacking','VdWContact','XBAcceptor','XBDonor']
            # First create empty df of interactions mapping to the whole receptor matching crystallographic residue numbering
            all_residues_plf_df = general_functions.create_prolif_reference_df(tleap_vs_cristal_reference_dict,interactions_list)
            
            # In case no MMGBSA computation was done, prepare the input prmtop and crd files for the ProLIF computation
            if mmgbsa == 1: 
                # If MMGSA computation was already done, use the prmtop and crd files generated by the MMPBSA computation which will include minimization for the corresponding pose
                pass
            else:
                # Compute the input prmtop and crd files for the ProLIF computation
                prmtop_file, crd_file = docking_analysis_utils.prepare_prolif_input_coordinates(output_path,complex_pdb_file,solvent,min_steps)
            
            # Retrieve the specific parameters for ProLIF computation
            prolif_parameters_db = self.docking_params_path + "/prolif_parameters.db"
            
            parameters_dict = docking_analysis_utils.retrieve_prolif_parameters_set(prolif_parameters_db,prolif_parameters_set)
            
            # Compute the fingerprint for the docked pose
            figerprints_df = docking_analysis_utils.compute_prolif_fps_for_docked_pose(prmtop_file,crd_file,interactions_list,tleap_vs_cristal_reference_dict, parameters_dict)
                
            # Map the fingerprint df to the crystallographic numbering
            df_mapped_to_cristal = general_functions.map_prolif_fingerprints_df_to_crystal_sequence(figerprints_df,tleap_vs_cristal_reference_dict)
            
            # Append the mapped row to the reference the whole protein interactions dataframe
            merged_df = general_functions.merge_calculated_and_reference_fingerprints_df(df_mapped_to_cristal,all_residues_plf_df)
    
            # Save the interactions df to the corresponding assay folder
            prolif_output_csv = f"{output_path}/prolif_fingerprints_renum.csv"
            merged_df.to_csv(prolif_output_csv,index=False)
    
            # Store the ProLIF results in the database
            docking_analysis_utils.store_prolif_fingerprints_results_in_db(assay_folder,assay_id,results_pose_id,ligname,sub_pose,complex_pdb_file,prolif_output_csv,prolif_parameters_set)
            
        
    ### Store the docked pose in the results database if required
        if store_docked_poses == 1:
            docking_analysis_utils.store_docked_pose_in_db(assay_folder,assay_id,results_pose_id,ligname,sub_pose,pose_pdb_file)
        
    ### Delete intermediate files after all calculations if required ###
        
        if clean_files == 1:
            md_utils.clean_MMPBSA_files(output_path)
        
    ### Delete the whole folder of computations if required
        if clean_folder == 1:
            shutil.rmtree(main_fingerprints_folder)
        
    ### Log the time at finish and calculate elapsed
        end_time = time.time()
        elapsed_time = f"{end_time - start_time:.2f}"
        
        print(f"Finished computing the fingerprint for the docked pose. - {elapsed_time} seconds")
        
    def compute_fingerprints_for_whole_assay(self, assay_id, mmgbsa=1, prolif=1, clean_files=0, clean_folder=0, solvent="implicit", min_steps=5000, stored_docked_poses=1, clean_assay_folder=1, prolif_parameters_set=1):
        assay_folder = self.docking_assays_path + f'/assay_{assay_id}'
        assay_results_db = f"{assay_folder}/assay_{assay_id}.db"
        
        docked_poses_list = docking_analysis_utils.retrieve_docked_poses_id(assay_results_db)
        
              
        ## Compute the fingerprints using a for loop:
        iteration = 1
        for pose in docked_poses_list:
            # Execute the fingerprint computation for each pose
            DockingAnalysis.compute_fingerprints_for_docked_pose(self, assay_id, pose, mmgbsa, prolif, clean_files, clean_folder, solvent, min_steps, stored_docked_poses, prolif_parameters_set, iteration=iteration)
            
            # Add to iteration counter
            iteration += 1

        if mmgbsa == 1:
            ## Sort the 'mmgbsa_fingerprints' table based on Pose_ID
            general_functions.sort_table(assay_folder,assay_id,"mmgbsa_fingerprints","Pose_ID")
            
        if prolif == 1:
            ## Sort the 'prolif_fingerprints' table based on Pose_ID
            general_functions.sort_table(assay_folder,assay_id,"prolif_fingerprints","Pose_ID")

        ## Delete the general fingerprint folders if required
        if clean_assay_folder == 1:
            shutil.rmtree(f"{assay_folder}/fingerprints_analyses")
    
    def set_prolif_custom_parameters(self,comment=None):
        """ This function allows the user to set custom parameters for ProLIF computation."""
        
        prolif_parameters_db = self.docking_params_path + "/prolif_parameters.db"    
        
        values_dict = {}
        
        while True:
            interaction = input("Input the interaction to parameterize ('f' to finish): ")
            if interaction and interaction.lower() != 'f':
                parameter = input("Input the parameter to set: ")
                value = input("Input the value for the parameter: ")
                # Evaluate the value to handle different types (int, float, str)
                processed_value = ast.literal_eval(value)
                values_dict[interaction] = {parameter: processed_value}
            elif interaction.lower() == 'f':
                if len(values_dict) == 0:
                    print("No parameters set. Storing an empty set.")
                    parameter = ''
                break
        
        docking_analysis_utils.save_prolif_parameters_set(prolif_parameters_db, parameter, values_dict,comment)
        
        
    
    