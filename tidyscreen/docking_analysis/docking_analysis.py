import warnings
warnings.filterwarnings("ignore")
from tidyscreen.docking_analysis import docking_analysis_utils as docking_analysis_utils
from tidyscreen.moldyn import moldyn_utils as md_utils
from tidyscreen.GeneralFunctions import general_functions as general_functions
import time
import shutil


class DockingAnalysis:
    
    def __init__(self, project):
        self.project = project
        self.docking_assays_path = self.project.proj_folders_path["docking"]["docking_assays"]
        self.docking_registers_path = self.project.proj_folders_path["docking"]["docking_registers"]    
        self.docking_params_path = self.project.proj_folders_path["docking"]["params"]    
        self.receptor_models_path = self.project.proj_folders_path["docking"]["receptors"]
        self.ligands_db = self.project.proj_folders_path["chemspace"]['processed_data'] + "/chemspace.db"
    
    def process_docking_assay(self,assay_id):
        registries_db = f"{self.docking_registers_path}/docking_registries.db"
        # Check if the 'assay_id' existis in docking registries
        docking_analysis_utils.check_docking_assay(registries_db,assay_id)
        # Process the corresponding 'dlg' files
        assay_folder = self.docking_assays_path + f'/assay_{assay_id}'
        results_db_file = docking_analysis_utils.process_dlg_files(assay_folder,assay_id)
        # Add docking subposes in database
        docking_analysis_utils.add_docking_subposes_nbr(results_db_file)
        # Extract 1 pdb pose per cluster
        docking_analysis_utils.extract_1_pdb_per_cluster(assay_folder,results_db_file)
    
    def compute_fingerprints_for_docked_pose(self,assay_id,results_pose_id,clean_files=1,clean_folder=1,solvent="implicit",min_steps=5000,store_docked_poses=1, iteration=1, ligname="UNL"):
    ### Start to log the time
        start_time = time.time()
    ### Create a custom folder for the analysis and copy/generate relevant files
        assay_folder = self.docking_assays_path + f'/assay_{assay_id}'
        complex_pdb_file, output_path, receptor_filename, ligname, sub_pose, pose_pdb_file  = docking_analysis_utils.create_fingerprints_analysis_folder(self,assay_folder,assay_id, results_pose_id)
        
        # Inform ligand under processing:
        print(f"\n Processing ligand: {ligname} \n")
    
    
    ### Compute fingerprints using MMPBSA on target folder
        # Compute the MMPBSA based per-residue interaction fingerprint
        main_fingerprints_folder = f"{assay_folder}/fingerprints_analyses"
        prmtop_file, crd_file, tleap_vs_cristal_reference_dict,mmpbsa_decomp_csv_output = docking_analysis_utils.compute_fingerprints(output_path,main_fingerprints_folder,complex_pdb_file,receptor_filename,solvent,min_steps,iteration,ligname)
        
    ### Compute ProLIFfingerprints for minimized pose ###
        interactions_list = ['Anionic','CationPi','Cationic','EdgeToFace','FaceToFace','HBAcceptor','HBDonor','Hydrophobic','MetalAcceptor','MetalDonor','PiCation','PiStacking','VdWContact','XBAcceptor','XBDonor']
        # First create empty df of interactions mapping to the whole receptor matching crystallographic residue numbering
        all_residues_plf_df = general_functions.create_prolif_reference_df(tleap_vs_cristal_reference_dict,interactions_list)
        
        # Compute the fingerprint for the docked pose
        figerprints_df = docking_analysis_utils.compute_prolif_fps_for_docked_pose(prmtop_file,crd_file,interactions_list,tleap_vs_cristal_reference_dict)
        
        # Map the fingerprint df to the crystallographic numbering
        df_mapped_to_cristal = general_functions.map_prolif_fingerprints_df_to_crystal_sequence(figerprints_df,tleap_vs_cristal_reference_dict)
        
        # Append the mapped row to the reference the whole protein interactions dataframe
        merged_df = general_functions.merge_calculated_and_reference_fingerprints_df(df_mapped_to_cristal,all_residues_plf_df)
        # Save the interactions df to the corresponding assay folder
        prolif_output_csv = f"{output_path}/prolif_fingerprints_renum.csv"
        merged_df.to_csv(prolif_output_csv,index=False)
        
    ### Store relevant computed files in the results database
        docking_analysis_utils.store_fingerprints_results_in_db(assay_folder,assay_id,results_pose_id,ligname,sub_pose,complex_pdb_file,mmpbsa_decomp_csv_output,prolif_output_csv)
    
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
        
    def compute_fingerprints_for_whole_assay(self,assay_id,clean_files=0,clean_folder=0,solvent="implicit",min_steps=5000,stored_docked_poses=1,clean_assay_folder=1):
        assay_folder = self.docking_assays_path + f'/assay_{assay_id}'
        assay_results_db = f"{assay_folder}/assay_{assay_id}.db"
        
        docked_poses_list = docking_analysis_utils.retrieve_docked_poses_id(assay_results_db)
              
        ## Compute the fingerprints using a for loop:
        iteration = 1
        for pose in docked_poses_list:
            # Execute the fingerprint computation for each pose
            DockingAnalysis.compute_fingerprints_for_docked_pose(self,assay_id,pose,clean_files,clean_folder,solvent,min_steps,stored_docked_poses,iteration=iteration)
            # Add to iteration counter
            iteration += 1
                        

        # Sort the 'fingerprints' table based on Pose_ID
        general_functions.sort_table(assay_folder,assay_id,"fingerprints","Pose_ID")

        ## Delete the general fingerprint folders if required
        if clean_assay_folder == 1:
            shutil.rmtree(f"{assay_folder}/fingerprints_analyses")
                
    