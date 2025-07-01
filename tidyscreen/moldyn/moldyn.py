import warnings
warnings.filterwarnings("ignore")
from tidyscreen.moldyn import moldyn_utils as moldyn_utils
from tidyscreen.docking_analysis import docking_analysis_utils as docking_analysis_utils
import sys
import os


class MolDyn:
    
    def __init__(self, project):
        self.project = project
        self.mdyn_path = self.project.proj_folders_path['dynamics']
        self.docking_path = self.project.proj_folders_path["docking"]
        self.docking_assays_path = self.project.proj_folders_path["docking"]["docking_assays"]
    
    def create_md_assay(self,docking_assay_id,docking_pose_id):
        
        try:
        
            # # Create the registry in the corresponding assays db
            registries_db = f"{self.mdyn_path['md_registers']}/md_registries.db"
            assay_id = moldyn_utils.create_md_assay_registry(registries_db,docking_assay_id,docking_pose_id)
            
            # # Create the corresponding assay folder and retrieve relevant data
            assay_folder = f"{self.mdyn_path['md_assays']}/assay_{assay_id}"
            moldyn_utils.create_md_assay_folder(assay_folder)
            
            ## Retrieve the .dlg to extract the selected binding pose
            md_folder = self.mdyn_path
            docking_folder = f"{self.docking_assays_path}/assay_{docking_assay_id}"
            ligname, sub_pose, dlg_file, run_number = docking_analysis_utils.retrieve_dlg_file(docking_folder,docking_assay_id,docking_pose_id)
            
            # Extract the .pdb file by parsing the 'run_number' in the 'dlg_file'
            pose_pdb_file = docking_analysis_utils.parse_dlg_by_run_number(ligname,dlg_file,run_number,assay_folder)
            
            # Get the receptor filename from the dlg file
            receptor_path = f"{docking_folder}/receptor"
            receptor_filename = docking_analysis_utils.get_receptor_name_from_dlg(dlg_file,receptor_path,assay_folder)
            
            # Generate the complexed file into de corresponding folder
            complex_pdb_file = docking_analysis_utils.generate_complex_from_pose(ligname,pose_pdb_file, assay_folder, receptor_filename,docking_pose_id)
            
            # Retrieve the table name matching docking assay
            table_name = docking_analysis_utils.retrieve_table_name_from_assay_id(self,docking_assay_id)
            
            # Retrieve .mol2 and .frcmod files
            docking_analysis_utils.retrieve_tleap_ligand_param_files(self,table_name,assay_folder,ligname,pdb=0,mol2_sybyl=0,mol2_gaff2=1,frcmod=1,pdbqt=0)
            
            # From the complex filename prepare the corresponding files naming:
            ligand_prefix = complex_pdb_file.split('/')[-1].split('_')[1]
            ligand_mol2_ref_file = f'{ligand_prefix}_gaff.mol2'
            ligand_frcmod_ref_file = f'{ligand_prefix}.frcmod'
            
            # Prepare the input files to compute prepare the .prmtop and .inpcrd
            moldyn_utils.prepare_md_initial_files(assay_folder,complex_pdb_file,ligand_mol2_ref_file,ligand_frcmod_ref_file, solvent="explicit",min_steps=5000,dynamics=1)
            
            # Execute tLeap initialization
            moldyn_utils.run_tleap_input(assay_folder,input_file='tleap.in')
        
            # Inform successful creation of the MD assay
            print(f"MD assay {assay_id} created successfully in {assay_folder}")
        
        except Exception as error:
            print("Error creating MD assays. Stopping...")
            print(error)
            sys.exit()
            
    def analyze_trajectory(self,assay_id):
        """
        Analyze the trajectory of a given MD assay.
        """
        try:
            # Retrieve the assay folder
            assay_folder = f"{self.mdyn_path['md_assays']}/assay_{assay_id}"
            
            
            # Check if the folder exists
            if os.path.isdir(assay_folder):
                pass
            else:
                print(f"Assay folder {assay_folder} does not exist. Stopping...")
                sys.exit()
            
            # Analyze the trajectory
            moldyn_utils.prepare_MD_folder_for_MMGBSA(assay_folder)
            
        except Exception as error:
            print("Error analyzing trajectory. Stopping...")
            print(error)
            sys.exit()