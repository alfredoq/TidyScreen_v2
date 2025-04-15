from tidyscreen.docking_analysis import docking_analysis_utils as docking_analysis_utils

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
    
    def compute_fingerprints_for_docked_pose(self,assay_id,results_pose_id,clean_files):
        # Create a custom folder for the analysis and copy/generate relevant files
        assay_folder = self.docking_assays_path + f'/assay_{assay_id}'
        complex_pdb_file, output_path, receptor_filename  = docking_analysis_utils.create_fingerprints_analysis_folder(self,assay_folder,assay_id, results_pose_id)
        # Compute fingerprints on target folder
        docking_analysis_utils.compute_fingerprints(output_path,complex_pdb_file,receptor_filename,clean_files)
        print("Finished computing the fingerprint for the docked pose.")