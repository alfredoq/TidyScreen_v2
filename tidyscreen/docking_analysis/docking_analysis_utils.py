import sys
from tidyscreen import tidyscreen as tidyscreen
from ringtail import RingtailCore
import pandas as pd
from pathlib import Path
import os
from tidyscreen.chemspace import chemspace as chemspace
from tidyscreen.moldyn import moldyn_utils as md_utils
import shutil

def check_docking_assay(registries_db,assay_id):
    # Check if the 'assay_id' exists in the registers database
    conn = tidyscreen.connect_to_db(registries_db)
    cursor = conn.cursor()
    cursor.execute(f"SELECT id FROM docking_registries WHERE id = {assay_id}")
    
    try:
        exists = cursor.fetchone()[0]  # 1 if exists, 0 if not
        pass
    except:
        print(f"Problem retrieving 'assay_id': {assay_id}. Stopping...")
        sys.exit()
        
def process_dlg_files(assay_folder,assay_id,max_poses=10):
    
    try:
        results_db_file = f"{assay_folder}/assay_{assay_id}.db"
        rtc = RingtailCore(results_db_file)
        rtc.add_results_from_files(file_path = f"{assay_folder}/dlgs/", recursive = True, save_receptor = False, max_poses = max_poses)
            
        return results_db_file
    
    except Exception as error:
        print(error)
        print("Error processing with Ringtail")        
        
def add_docking_subposes_nbr(results_db_file):
    conn = tidyscreen.connect_to_db(results_db_file)
    
    try:
        sql = """SELECT *
               FROM Results;"""
        
        df = pd.read_sql(sql,conn)

        # Add the '_#' value to the corresponding name

        unique_names_list = []
        for index, row in df.iterrows():
            ligname = row['LigName']
            if ligname not in unique_names_list:
                counter=1
                df.at[index,'sub_pose'] =  ligname+'_'+str(counter)
                unique_names_list.append(ligname)
            else:
                counter+=1 
                df.at[index,'sub_pose'] =  ligname+'_'+str(counter)

        # Reorder df columns and store to the database replacing the original one
        df = df[['Pose_ID', 'LigName', 'sub_pose','receptor', 'pose_rank', 'run_number','docking_score', 'leff', 'deltas', 'cluster_rmsd', 'cluster_size','reference_rmsd', 'energies_inter', 'energies_vdw', 'energies_electro','energies_flexLig', 'energies_flexLR', 'energies_intra','energies_torsional', 'unbound_energy', 'nr_interactions', 'num_hb','about_x', 'about_y', 'about_z', 'trans_x', 'trans_y', 'trans_z','axisangle_x', 'axisangle_y', 'axisangle_z', 'axisangle_w', 'dihedrals','ligand_coordinates', 'flexible_res_coordinates']]
        
        save_df_to_db(results_db_file,df,"Results")

    except Exception as error:
        print(error)
        
def save_df_to_db(db,df,table_name,action="replace"):
    conn = tidyscreen.connect_to_db(db)
    df.to_sql(con=conn, name=table_name,if_exists=action,index=None)

def extract_1_pdb_per_cluster(assay_folder,results_db_file):
    """
    This function will parse a .dlg file and will extract 1 .pdb per identified cluster into a folder named 'cluster_pdb_files' located within the docking assay folder
    """
    conn = tidyscreen.connect_to_db(results_db_file)
    
    # Configure the output path to store corresponding .sdf files    
    output_path = f'{assay_folder}/docked_1_per_cluster'
    Path(output_path).mkdir(parents=True, exist_ok=True)
    
    # Retrieve the names of ligand in order to get the .dlg to process and that contains the ligands
    sql = f"""SELECT LigName FROM Ligands;"""
    ligands_name_df = pd.read_sql_query(sql,conn)
    
    activation_keywords =  ['CLUSTERING', 'HISTOGRAM'] # This is the opening line of the Clustering Histogram
    shutdown_keywords = ['RMSD', 'TABLE']
    number_of_poses_to_extract = 999
    
    for index, row in ligands_name_df.iterrows():
        ligname = row['LigName']
        dlg_file_to_parse = f'{assay_folder}/dlgs/{ligname}.dlg'
        
        trigger = 0 
        
        # This loop will construct the list of docked poses to store based on the clustering histogram
        counter = 1
        pose_rank_list = []
        pose_id_list = []

        for line in open(dlg_file_to_parse):
            new_line = line.rsplit()
            if new_line == activation_keywords:
                trigger = 1
            elif new_line == shutdown_keywords:
                trigger = 0
    
            if trigger == 1 and len(new_line) == 10 and counter <= number_of_poses_to_extract and "#" in new_line[-1]:
                pose_rank_list.append(int(new_line[0]))
                pose_id_list.append(new_line[4])
                counter+=1        

        # This will extract all the pre-identified pdb files and store them into the indicated folder
        for pose_rank in pose_rank_list:
            output_file = f'{output_path}/{ligname}_{pose_rank}.pdb'
            pose_id = pose_id_list[pose_rank-1]
            trigger_2 = 0
            activation_keywords_2 =  ['DOCKED:', 'USER', 'Run', '=', pose_id] # This is the opening line of the .pdbqt run to extract
            shutdown_keywords_2 = ['DOCKED:', 'ENDMDL']
            with open(output_file,'w') as pdb_output:
                for line in open(dlg_file_to_parse):
                    new_line = line.rsplit()
                    if new_line == activation_keywords_2:
                        trigger_2 = 1
                    elif new_line == shutdown_keywords_2:
                        trigger_2 = 0
        
                    if trigger_2 == 1 and new_line[1] == 'ATOM': # This extract the lines constituting the .pdb file
                        # This will output the .pdb file with the corresponding format
                        #pdb_output.write(f"{new_line[1]}{new_line[2]:>7}{'':<2}{new_line[3]:<4}{new_line[4]:<8}{new_line[5]:<7}{new_line[6]:<8}{new_line[7]:<8}{new_line[8]:<8}\n")
                        atom_field = "HETATM"
                        pdb_output.write(f"{atom_field:<6}{new_line[2]:>5}{'':<1}{new_line[3]:>4}{'':<1}{new_line[4]:>3}{'':<2}{new_line[5]:>4}{'':<4}{new_line[6]:>8}{new_line[7]:>8}{new_line[8]:>8}\n")
        
            pdb_output.close()
        
def create_fingerprints_analysis_folder(self,assay_folder,assay_id,results_pose_id):
    # Creat the for to store all files
    output_path = f'{assay_folder}/fingerprints_analyses/pose_{results_pose_id}'
    Path(output_path).mkdir(parents=True, exist_ok=True)
    # Retrieve the dlg corresponding to the ligand
    ligname, dlg_file, run_number = retrieve_dlg_file(assay_folder,assay_id,results_pose_id)
    # Extract the .pdb file by parsing the 'run_number' in the 'dlg_file'
    pose_pdb_file = parse_dlg_by_run_number(ligname,dlg_file,run_number,output_path)
    # Get the receptor filename from the dlg file
    receptor_path = f'{assay_folder}/receptor'
    receptor_filename = get_receptor_name_from_dlg(dlg_file,receptor_path,output_path)
    print(receptor_filename)
    # Generate the complexed file into de corresponding folder
    complex_pdb_file = generate_complex_from_pose(ligname,pose_pdb_file, output_path, receptor_filename,results_pose_id)
    # Retrieve table name matching docking assay
    table_name = retrieve_table_name_from_assay_id(self,assay_id)
    # Retrieve .mol2 and .frcmod files
    retrieve_tleap_ligand_param_files(self,table_name,output_path,ligname,pdb=1,mol2_sybyl=1,mol2_gaff2=1,frcmod=1,pdbqt=1)
    
    return complex_pdb_file, output_path, receptor_filename
    
def retrieve_dlg_file(assay_folder,assay_id,results_pose_id):
    """
    This function will return the name of the .dlg file corresponding the provided 'results_pose_id', as well as the 'run_number' corresponding the 'pose_id'
    """
    
    results_db_file = f"{assay_folder}/assay_{assay_id}.db"
    print(results_db_file)
    conn = tidyscreen.connect_to_db(results_db_file)
    cursor = conn.cursor()
    cursor.execute(f"SELECT LigName, run_number FROM Results WHERE Pose_ID = {results_pose_id}")
    
    try:
        
        data = cursor.fetchall()
        ligname = data[0][0]
        dlg_file = assay_folder + '/dlgs/' + ligname + '.dlg'  
        run_number = data[0][1]
        return ligname, dlg_file, run_number
        
    except Exception as error:
        print(error)
        print(f"Error retrieving pose '{results_pose_id}'. Stopping...")
        sys.exit()
    
def parse_dlg_by_run_number(ligname,dlg_file,run_number,output_path):
    
    output_filename = f'{output_path}/{ligname}_{run_number}.pdb'
    trigger = 0
    with open(dlg_file,'r') as dlg_file, open(output_filename,'w') as output_file:
        for line in dlg_file:
            new_line = line.rsplit()
            if len(new_line) == 4 and new_line[0] == 'Run:' and int(new_line[1]) == run_number:
                trigger = 1
    
            if len(new_line) == 2 and new_line[1] == 'ENDMDL' and trigger == 1:
                trigger = 0

            if trigger == 1 and len(new_line) > 5 and new_line[1] == 'ATOM':
                atom_field = "HETATM"
                pdb_string = f"{atom_field:<6}{new_line[2]:>5}{'':<1}{new_line[3]:>4}{'':<1}{new_line[4]:>3}{'':<2}{new_line[5]:>4}{'':<4}{new_line[6]:>8}{new_line[7]:>8}{new_line[8]:>8}"
                output_file.write(f"{pdb_string}\n")

    return output_filename

def get_receptor_name_from_dlg(dlg_file,receptor_path,output_path):
        
    with open(dlg_file,'r') as input_file:
        for line in input_file:
            line_split = line.rsplit()

            if len(line_split) > 2 and line_split[0] == "Receptor" and line_split[1] == "name:":
                receptor_pdb_filename = line_split[2]+'.pdb'
    
    receptor_filename = f"{receptor_path}/{receptor_pdb_filename}"
    
    # Check if the 'receptor_filename' retrieved exists:
    if os.path.exists(receptor_filename):
        pass
    else:
        print(f"The receptor {receptor_pdb_filename} does not exist. Stopping...")
    
    # Copy the receptor .pdb file to the fingerprint folder
    #shutil.copy(receptor_filename, output_path)
    
    return receptor_filename

def generate_complex_from_pose(ligname,pose_pdb_file, output_path, receptor_filename,results_pose_id):
    
    complex_output_file = f"{output_path}/complex_{ligname}_pose_{results_pose_id}.pdb"
    
    with open(receptor_filename,'r') as receptor_file, open(complex_output_file, "w") as outfile:
       lines =receptor_file.readlines()
       for line in lines:
           if line.strip() != "END":
               outfile.write(line)
    outfile.close()
    
    # Check of the last line at this point is 'TER'. If not, add it
    with open(complex_output_file, "r") as outfile:
        lines = outfile.readlines()
        if lines and lines[-1].strip() == "TER":
            add_ter = 0
        else:
            add_ter = 1
    outfile.close()
    
    if add_ter == 1:
        with open(complex_output_file, "a") as outfile:
            outfile.write("TER\n")
    
    # Read and write the contents of the pdb_pose file
    with open(pose_pdb_file, "r") as f2, open(complex_output_file,"a") as outfile:
        outfile.write(f2.read())
    outfile.close()
    
    # Check of the last line at this point is 'TER'. If not, add it
    with open(complex_output_file, "r") as outfile:
        lines = outfile.readlines()
        if lines and lines[-1].strip() == "TER":
            add_ter = 0
        else:
            add_ter = 1
    outfile.close()
    
    if add_ter == 1:
        with open(complex_output_file, "a") as outfile:
            outfile.write("TER\n")
            outfile.write("END\n")

    return complex_output_file
            
def retrieve_tleap_ligand_param_files(self,table_name,output_path,ligname,pdb,mol2_sybyl,mol2_gaff2,frcmod,pdbqt):
    chemspace.ChemSpace.retrieve_mols_in_table(self,table_name,output_path,ligname,pdb,mol2_sybyl,mol2_gaff2,frcmod,pdbqt)

def retrieve_table_name_from_assay_id(self,assay_id):
    registries_db = f"{self.project.proj_folders_path['docking']['docking_registers']}/docking_registries.db"
    conn = tidyscreen.connect_to_db(registries_db)
    
    cursor = conn.cursor()
    cursor.execute(f"SELECT table_name FROM docking_registries WHERE id = {assay_id}")
    table_name = cursor.fetchone()[0]
    
    return table_name
    
def compute_fingerprints(assay_folder,complex_pdb_file,receptor_filename,clean_files,solvent):
    # From the complex filename prepare the corresponding files naming:
    ligand_prefix = complex_pdb_file.split('/')[-1].split('_')[1]
    ligand_mol2_ref_file = f'{ligand_prefix}_gaff.mol2'
    ligand_frcmod_ref_file = f'{ligand_prefix}.frcmod'
    # Prepare the input files to compute prepare the .prmtop and .inpcrd
    md_utils.prepare_md_initial_files(assay_folder,complex_pdb_file,ligand_mol2_ref_file,ligand_frcmod_ref_file, solvent, dynamics=0)
    # Execute tLeap initialization
    md_utils.run_tleap_input(assay_folder,input_file='tleap.in')
    # Perform minimizations
    md_utils.perform_minimization(assay_folder,ligand_prefix,solvent)
    # Create the MMGBSA input file
    md_utils.write_mmgbsa_input(assay_folder)
    # Use ante ante-MMPBSA.py to create the files required for computation
    md_utils.apply_ante_MMPBSA(assay_folder)
    
    # Strip solven (if explicit model was used) and/or compute the MMPBSA 
    if solvent == "explicit":
        # Strip waters from min2.crd to
        md_utils.strip_waters(assay_folder,"min2.crd","complex.prmtop")
        # perform the MMPBSA analysis
        assay_folder, decomp_file = md_utils.compute_MMPBSA(assay_folder,"min2_strip.crd")
    
    if solvent == "implicit":
        # perform the MMPBSA analysis
        assay_folder, decomp_file = md_utils.compute_MMPBSA(assay_folder,"min1.crd")
    
    # perform a cleaning of the target directory if required
    if clean_files == 1:
        md_utils.clean_MMPBSA_files(assay_folder)

    # renumber the output of MMGBSA according to the original receptor.
    md_utils.renumber_mmgbsa_output(assay_folder,decomp_file,receptor_filename)