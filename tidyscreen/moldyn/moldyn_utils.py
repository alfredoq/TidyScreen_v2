import warnings
warnings.filterwarnings("ignore")
import shutil
import subprocess
import glob
import os
from tidyscreen.GeneralFunctions import general_functions as general_functions
import sys
import pandas as pd
from tidyscreen import tidyscreen
import json


def prepare_md_initial_files(output_path,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm,solvent,min_steps,dynamics=1):
    ## Prepare the corresponding 'tleap' input file
    # These preparations is common for fingerprints and MD assays
    
    if solvent == "explicit":
        prepare_tleap_input(output_path,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm)
        prepare_min1_input(output_path,min_steps)
        prepare_min2_input(output_path,min_steps)
    
    if solvent == "implicit":
        prepare_tleap_input_implicit_solvent(output_path,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm)
        prepare_min1_input_implicit_solvent(output_path,min_steps)

    if dynamics == 1: # This will only be executed for simulations involving MD preparations
        prepare_heat1_input(output_path)
        prepare_heat2_input(output_path)
        prepare_equil_input(output_path)
        prepare_prod_input(output_path)
        prepare_md_execution_script(output_path)
    
def prepare_tleap_input(output_path,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm):
    with open(f"{output_path}/tleap.in","w") as tleap_input:
        tleap_input.write("source leaprc.protein.ff14SB \n")
        tleap_input.write("source leaprc.gaff \n")
        tleap_input.write("source leaprc.water.tip3p \n")
        tleap_input.write(f"UNL = loadmol2 {output_path}/{mol2_lig_parm} \n")
        tleap_input.write(f"loadamberparams {output_path}/{frcmod_lig_parm} \n")
        tleap_input.write(f"COM = loadpdb {complex_pdb_file} \n")
        tleap_input.write(f"addions COM Na+ 0 \n")
        tleap_input.write(f"addions COM Cl- 0 \n")
        tleap_input.write("solvateoct COM TIP3PBOX 10.0 \n")
        tleap_input.write(f"saveamberparm COM {output_path}/complex.prmtop {output_path}/complex.inpcrd \n")
        tleap_input.write("quit")
    
    tleap_input.close()

def prepare_tleap_input_implicit_solvent(output_path,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm):

    with open(f"{output_path}/tleap.in","w") as tleap_input:
        tleap_input.write("source leaprc.protein.ff14SB \n")
        tleap_input.write("source leaprc.gaff \n")
        tleap_input.write("set default PBRadii bondi \n")
        tleap_input.write(f"UNL = loadmol2 {output_path}/{mol2_lig_parm} \n")
        tleap_input.write(f"loadamberparams {output_path}/{frcmod_lig_parm} \n")
        tleap_input.write(f"COM = loadpdb {complex_pdb_file} \n")
        tleap_input.write(f"saveamberparm COM {output_path}/complex.prmtop {output_path}/complex.inpcrd \n")
        tleap_input.write("quit")
    
    tleap_input.close()

def prepare_min1_input(output_path,min_steps):
    with open(f"{output_path}/min1.in","w") as min1_input:    
        min1_input.write("Initial minimization of water molecules\n")
        min1_input.write("&cntrl\n")
        min1_input.write("imin = 1,\n")
        min1_input.write(f"maxcyc = {min_steps},\n")
        min1_input.write(f"ncyc = {int(min_steps/2)},\n")
        min1_input.write("ntxo = 1,\n")
        min1_input.write("ntb = 1,\n")
        min1_input.write("ntr = 1,\n")
        min1_input.write("cut = 8.0,\n")
        min1_input.write("restraintmask = '(:* & !:WAT)',\n")
        min1_input.write("restraint_wt = 500,\n")
        min1_input.write("&end\n")
        min1_input.write("END\n")
    
    min1_input.close()

def prepare_min1_input_implicit_solvent(output_path,min_steps):

    with open(f"{output_path}/min1.in","w") as min1_input:    
        min1_input.write("Initial minimization of water molecules\n")
        min1_input.write("&cntrl\n")
        min1_input.write("imin = 1,\n")
        min1_input.write(f"maxcyc = {min_steps},\n")
        min1_input.write(f"ncyc = {int(min_steps/2)},\n")
        min1_input.write("ntx = 1,\n")
        min1_input.write("ntxo = 1,\n")
        min1_input.write("igb = 8,\n")
        min1_input.write("gbsa = 3,\n")
        min1_input.write("ntb = 0,\n")
        min1_input.write("ntr = 0,\n")
        min1_input.write("cut = 999,\n")
        min1_input.write("surften = 0.07,\n")
        min1_input.write("saltcon = 0.0,\n")
        min1_input.write("&end\n")
        min1_input.write("END\n")
    
    min1_input.close()

def prepare_min2_input(output_path,min_steps):
    with open(f"{output_path}/min2.in","w") as min2_input:
        min2_input.write("Initial minimization of water molecules\n")
        min2_input.write("&cntrl\n")
        min2_input.write("imin = 1,\n")
        min2_input.write(f"maxcyc = {min_steps},\n")
        min2_input.write(f"ncyc = {int(min_steps)},\n")
        min2_input.write("ntxo = 1,\n")
        min2_input.write("ntb = 1,\n")
        min2_input.write("ntr = 0,\n")
        min2_input.write("cut = 8.0,\n")
        min2_input.write("&end\n")
        min2_input.write("END\n")
    
    min2_input.close()

def prepare_heat1_input(output_path):
    
    with open(f"{output_path}/heat1.in","w") as heat1_input:
        heat1_input.write(f"First heating stage of the solvent, with restraints in the solute\n")
        heat1_input.write("&cntrl\n")
        heat1_input.write("imin = 0,\n")
        heat1_input.write("irest = 0,\n")
        heat1_input.write("ntx = 1,\n")
        heat1_input.write("ntb = 2,\n")
        heat1_input.write("ntp = 1,\n")
        heat1_input.write("cut = 10.0,\n")
        heat1_input.write("ntr = 1,\n")
        heat1_input.write("ntc = 2,\n")
        heat1_input.write("ntf = 2,\n")
        heat1_input.write("tempi = 0.0,\n") 
        heat1_input.write("temp0 = 310.0,\n")
        heat1_input.write("ntt = 3,\n")
        heat1_input.write("gamma_ln = 1.0,\n")
        heat1_input.write("nstlim = 200000,\n")
        heat1_input.write("dt = 0.002,\n")
        heat1_input.write("ntpr = 100,\n")
        heat1_input.write("ntwx = 100,\n")
        heat1_input.write("ntwr = 1000,\n")
        heat1_input.write("restraintmask = '(:* & !:WAT & !:Na+ & !:Cl-)',\n")
        heat1_input.write("restraint_wt = 50.0,\n")
        heat1_input.write("nmropt = 1,\n")
        heat1_input.write("&end\n")
        heat1_input.write("&wt\n")
        heat1_input.write("TYPE = 'TEMP0',\n")
        heat1_input.write("ISTEP1 = 0,\n")
        heat1_input.write("ISTEP2 = 100000,\n")
        heat1_input.write("VALUE1 = 0.0,\n")
        heat1_input.write("VALUE2 = 310.0,\n")
        heat1_input.write("/\n")
        heat1_input.write("&wt TYPE = 'END'\n")
        heat1_input.write("/\n")
        heat1_input.write("END\n")
    
    heat1_input.close()
    
def prepare_heat2_input(output_path):
   
    with open(f"{output_path}/heat2.in","w") as heat2_input:
        heat2_input.write(f"Second heating stage of the solvent and solute, with restraints in the protein backbone\n")
        heat2_input.write("&cntrl\n")
        heat2_input.write("imin = 0,\n")
        heat2_input.write("irest = 0,\n")
        heat2_input.write("ntx = 1,\n")
        heat2_input.write("ntb = 2,\n")
        heat2_input.write("ntp = 1,\n")
        heat2_input.write("cut = 10.0,\n")
        heat2_input.write("ntr = 1,\n")
        heat2_input.write("ntc = 2,\n")
        heat2_input.write("ntf = 2,\n")
        heat2_input.write("tempi = 310.0,\n") 
        heat2_input.write("temp0 = 310.0,\n")
        heat2_input.write("ntt = 3,\n")
        heat2_input.write("gamma_ln = 1.0,\n")
        heat2_input.write("nstlim = 200000,\n")
        heat2_input.write("dt = 0.002,\n")
        heat2_input.write("ntpr = 100,\n")
        heat2_input.write("ntwx = 100,\n")
        heat2_input.write("ntwr = 1000,\n")
        heat2_input.write("restraintmask = '(@C,CA,N,O & !:WAT)',\n")
        heat2_input.write("restraint_wt = 5.0,\n")
        heat2_input.write("nmropt = 0,\n")
        heat2_input.write("&end\n")
        heat2_input.write("/\n")
        heat2_input.write("END\n")
    
    heat2_input.close()
   
def prepare_equil_input(output_path):
   
    with open(f"{output_path}/equi.in","w") as equi_input:
        equi_input.write(f"Equilibration stage in which only the protein backbone is restrained\n")
        equi_input.write(f"&cntrl\n")
        equi_input.write(f"imin = 0,\n") 
        equi_input.write(f"irest = 1,\n")
        equi_input.write(f"ntx = 5,\n")
        equi_input.write(f"ntb = 1,\n")
        equi_input.write(f"ntp = 0,\n")
        equi_input.write(f"cut = 10.0,\n")
        equi_input.write(f"ntr = 1,\n")
        equi_input.write(f"ntc = 2,\n")
        equi_input.write(f"ntf = 2,\n")
        equi_input.write(f"tempi = 310.0,\n")
        equi_input.write(f"temp0 = 310.0,\n")
        equi_input.write(f"ntt = 3,\n")
        equi_input.write(f"gamma_ln = 1.0,\n")
        equi_input.write(f"nstlim = 500000,\n")
        equi_input.write(f"dt = 0.002,\n")
        equi_input.write(f"ntpr = 100,\n")
        equi_input.write(f"ntwx = 100,\n")
        equi_input.write(f"ntwr = 1000,\n")
        equi_input.write(f"restraintmask = '(@C,CA,N,O & !:WAT)',\n")
        equi_input.write(f"restraint_wt = 0.5,\n")
        equi_input.write(f"&end\n")
        equi_input.write(f"END\n")
    
    equi_input.close()
    
def prepare_prod_input(output_path):
    
    with open(f"{output_path}/prod.in","w") as prod_input:
        prod_input.write(f"Production stage in which only the protein backbone is restrained\n")
        prod_input.write(f"&cntrl\n")
        prod_input.write(f"imin = 0,\n") 
        prod_input.write(f"irest = 1,\n")
        prod_input.write(f"ntx = 5,\n")
        prod_input.write(f"ntb = 1,\n")
        prod_input.write(f"ntp = 0,\n")
        prod_input.write(f"cut = 10.0,\n")
        prod_input.write(f"ntr = 1,\n")
        prod_input.write(f"ntc = 2,\n")
        prod_input.write(f"ntf = 2,\n")
        prod_input.write(f"tempi = 310.0,\n")
        prod_input.write(f"temp0 = 310.0,\n")
        prod_input.write(f"ntt = 3,\n")
        prod_input.write(f"gamma_ln = 1.0,\n")
        prod_input.write(f"nstlim = 5000000,\n")
        prod_input.write(f"dt = 0.002,\n")
        prod_input.write(f"ntpr = 100,\n")
        prod_input.write(f"ntwx = 100,\n")
        prod_input.write(f"ntwr = 1000,\n")
        prod_input.write(f"restraintmask = '(@C,CA,N,O & !:WAT)',\n")
        prod_input.write(f"restraint_wt = 0.5,\n")
        prod_input.write(f"&end\n")
        prod_input.write(f"END\n")
        
def prepare_md_execution_script(output_path):
    pmemd_cuda_path = shutil.which('pmemd.cuda')
    with open(f'{output_path}/md_execution.sh','a') as exec_file:
        exec_file.write("echo 'Executing Min1'\n")
        exec_file.write(f"{pmemd_cuda_path} -O -i min1.in -o min1.out -p complex.prmtop -c complex.inpcrd -r min1.crd -ref complex.inpcrd \n")
        exec_file.write("echo 'Executing Min2'\n")
        exec_file.write(f"{pmemd_cuda_path} -O -i min2.in -o min2.out -p complex.prmtop -c min1.crd -r min2.crd -ref min1.crd \n")
        exec_file.write("echo 'Executing Heat1'\n")
        exec_file.write(f"{pmemd_cuda_path} -O -i heat1.in -o heat1.out -p complex.prmtop -c min2.crd -r heat1.crd -ref min2.crd \n")
        exec_file.write("echo 'Executing Heat2'\n")
        exec_file.write(f"{pmemd_cuda_path} -O -i heat2.in -o heat2.out -p complex.prmtop -c heat1.crd -r heat2.crd -ref heat1.crd \n")
        exec_file.write("echo 'Executing Equilibration'\n")
        exec_file.write(f"{pmemd_cuda_path} -O -i equi.in -o equi.out -p complex.prmtop -c heat2.crd -r equi.crd -ref heat2.crd \n")
        exec_file.write("echo 'Executing Production'\n")
        exec_file.write(f"{pmemd_cuda_path}  -O -i prod.in -o prod.out -p complex.prmtop -c equi.crd -r prod.crd -x prod.nc -ref equi.crd \n")

    exec_file.close()

def prepare_md_execution_script_mendieta(output_path,mendieta_assays_path="NONE"):
    with open(f'{output_path}/md_execution_mendieta.sh','a') as exec_file:
        exec_file.write("#!/bin/bash \n")
        exec_file.write("### Nombre de la tarea\n")
        exec_file.write("#SBATCH --job-name=nombre\n")
        exec_file.write("### Cantidad de nodos a usar\n")
        exec_file.write("#SBATCH --nodes=1\n")
        exec_file.write("### GPUs por nodo (<= 2)\n")
        exec_file.write("### OJO: Todos los procesos en un nodo ven ambas GPU con gpu:2!\n")
        exec_file.write("#SBATCH --gres=gpu:1\n")
        exec_file.write("### Procesos por nodo\n")
        exec_file.write("#SBATCH --ntasks-per-node=1\n")
        exec_file.write("### Cores por proceso (OpenMP/Pthreads/etc)\n")
        exec_file.write("### Recordar exportar OMP_NUM_THREADS/MKL_NUM_THREADS/etc con el mismo valor\n")
        exec_file.write("#SBATCH --cpus-per-task=1\n")
        exec_file.write("export OMP_NUM_THREADS=1\n")
        exec_file.write("### Tiempo de ejecucion. Formato dias-horas:minutos. Maximo una semana.\n")
        exec_file.write("#SBATCH --time 7-0:00\n")
        exec_file.write("### Environment setup\n")
        exec_file.write(". /etc/profile\n")
        exec_file.write("### Environment modules\n")
        exec_file.write("module load cuda/5.0\n")
        exec_file.write("### Ejecutar la tarea\n")
        exec_file.write("srun algun_programa_gpu\n")

    exec_file.close()

def run_tleap_input(output_path,input_file):
    
    tleap_path = shutil.which('tleap')
    command = f'{tleap_path} -f {output_path}/{input_file}' 
    subprocess.run(command, shell=True, capture_output=True, text=True)
    
def perform_minimization(output_path,ligand_prefix,solvent,inform=1):
    # Compute mol2 with Sybyl Atom Types - for compatibility with RDKit and Meeko
    pmemd_path = shutil.which('pmemd.cuda')
    # Run first minimization procedure
    try:
        command = f'{pmemd_path} -O -i {output_path}/min1.in -o {output_path}/min1.out -p {output_path}/complex.prmtop -c {output_path}/complex.inpcrd -r {output_path}/min1.crd -ref {output_path}/complex.inpcrd' 
        if inform == 1:
            print(f"Runnning 'min1' for ligand: '{ligand_prefix}'")
        subprocess.run(command, shell=True, capture_output=True, text=True)
    except Exception as error:
        print(error)
    
    if solvent == "explicit":
        # Run second minimization procedure
        command = f'{pmemd_path} -O -i {assay_folder}/min2.in -o {assay_folder}/min2.out -p {assay_folder}/complex.prmtop -c {assay_folder}/min1.crd -r {assay_folder}/min2.crd' 
        if inform == 1:
            print(f"Runnning 'min2' for ligand: '{ligand_prefix}'")
        subprocess.run(command, shell=True, capture_output=True, text=True)
    
def write_mmgbsa_input(output_path,interval=100):
    with open(f'{output_path}/mmgbsa.in','w') as mmgbsa_file:
        mmgbsa_file.write("Per-residue GB and PB decomposition \n")
        mmgbsa_file.write("&general \n")
        mmgbsa_file.write(f"startframe = 1, interval = {interval}, verbose=1, \n")
        mmgbsa_file.write("/ \n")
        mmgbsa_file.write("&gb \n")
        mmgbsa_file.write("igb=5, saltcon=0.100, \n")
        mmgbsa_file.write("/ \n")
        mmgbsa_file.write("&decomp \n")
        mmgbsa_file.write("idecomp=2, csv_format=0, \n")
        mmgbsa_file.write("/ \n")
        
    mmgbsa_file.close()

def strip_waters(output_path,traj_file,prmtop_file):
    traj_file_prefix = traj_file.split('.')[0]
    traj_file_extension = traj_file.split('.')[1]
    with open(f'{output_path}/cpptraj_strip_wat.in','w') as cpptraj_file:
        cpptraj_file.write(f"parm {output_path}/{prmtop_file}\n")
        cpptraj_file.write(f"trajin {output_path}/{traj_file}\n")
        cpptraj_file.write(f"strip :WAT,Na+,Cl- \n")
        cpptraj_file.write(f"trajout {output_path}/{traj_file_prefix}_strip.{traj_file_extension} trajectory  \n")
        cpptraj_file.write(f"go \n")
        cpptraj_file.write(f"quit \n")
    
    cpptraj_file.close()
    
    # Determine the ante_MMPBSA path
    cpptraj_path = shutil.which('cpptraj')
    # Run cpptraj computation to strip water molecules
    command = f'{cpptraj_path} -i {output_path}/cpptraj_strip_wat.in' 
    subprocess.run(command, shell=True, capture_output=True, text=True)
    
def apply_ante_MMPBSA(output_path,ligresname,amberhome):
    # Determine the ante_MMPBSA path
    anteMMPBSA_path = shutil.which('ante-MMPBSA.py')
    # Set $AMBERHOME environment variable
    os.environ["AMBERHOME"] = amberhome
    # Run ante_MMPBSA computation
    command = f'{anteMMPBSA_path} -p {output_path}/complex.prmtop -s :WAT,Na+,Cl- -c {output_path}/complex_MMGBSA.prmtop -r {output_path}/receptor_MMGBSA.prmtop -l {output_path}/ligand_MMGBSA.prmtop -n :{ligresname}' 
    
    print(f"Computing ante_MMPBSA")
    
    subprocess.run(command, shell=True, capture_output=True, text=True)

def compute_MMPBSA(output_path,traj_input):
    
    # Check if the MMPBSA.py script is available in the conda environment, if so delete it so as to use the system one
    conda_prefix = os.environ.get('CONDA_PREFIX')
    # Path to the file in the bin directory (replace 'my_executable' with your filename)
    file_path = os.path.join(conda_prefix, 'bin', 'MMPBSA.py')
    if os.path.isfile(file_path):
        os.remove(file_path)
        print(f"Removed MMPBSA.py from conda environment. Using system version instead.")
    
    try: 
        # Get the MMPBSA.py path
        MMPBSA_path = shutil.which('MMPBSA.py')
        # Run ante_MMPBSA computation
        command = f'{MMPBSA_path} -i {output_path}/mmgbsa.in -o {output_path}/mmgbsa.out -do {output_path}/mmgbsa_DECOMP.out -cp {output_path}/complex_MMGBSA.prmtop -rp {output_path}/receptor_MMGBSA.prmtop -lp {output_path}/ligand_MMGBSA.prmtop -y {output_path}/{traj_input}'
        print(f"Computing MMPBSA")
        subprocess.run(command, shell=True, capture_output=True, text=True)
        
        # Clean the temporal files
        command = f'{MMPBSA_path} -i {output_path}/mmgbsa.in -o {output_path}/mmgbsa.out -do {output_path}/mmgbsa_DECOMP.out -cp {output_path}/complex_MMGBSA.prmtop -rp {output_path}/receptor_MMGBSA.prmtop -lp {output_path}/ligand_MMGBSA.prmtop -y {output_path}/{traj_input} --clean'
        
        subprocess.run(command, shell=True, capture_output=True, text=True)
        
        return output_path, f"{output_path}/mmgbsa_DECOMP.out"
    
    except Exception as error:
        print(f"Error during MMPBSA computation: {error}")
        return None, None
    
def clean_MMPBSA_files(target_dir):
    patterns_to_retain = ['renum.csv','_pose_']
    
    # Construct a list of files to retain
    files_to_retain = []
    for pattern in patterns_to_retain:
        for filename in os.listdir(target_dir):
            if pattern in filename:
                files_to_retain.append(filename)

    # Delete file if not marked to be retained
    for filename in os.listdir(target_dir):
        if filename not in files_to_retain:
            file_path = os.path.join(target_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
    
    
    # for pattern in patterns_to_clean:
    #     for filename in os.listdir(target_dir):
    #         if pattern in filename:
    #             for pattern_to_exclude in patterns_to_retain:
    #                 if pattern_to_exclude not in filename
    #                     file_path = os.path.join(target_dir, filename)
    #                 if os.path.isfile(file_path):
    #                 os.remove(file_path)

def renumber_mmgbsa_output(output_path,main_fingerprints_folder,decomp_file,receptor_filename,iteration):
    
    # # Parse iformation regarding the receptor file to manage renumbering
    # if iteration == 1:
    #     resname_field,resnumber_field,chain_name_field = parse_receptor_fields(receptor_filename, main_fingerprints_folder)
    
    # else:
    #     # This will load the parameters from the 'params.json' file
    #     with open(f"{main_fingerprints_folder}/receptor_fields.json",'r') as file:
    #         loaded_params = json.load(file)

    #     resname_field = loaded_params['resname_field']
    #     resnumber_field = loaded_params['resnumber_field']
    #     chain_name_field = loaded_params['chain_name_field']
    
    # This will load the parameters from the 'params.json' file
    with open(f"{main_fingerprints_folder}/receptor_fields.json",'r') as file:
            loaded_params = json.load(file)

    resname_field = loaded_params['resname_field']
    resnumber_field = loaded_params['resnumber_field']
    chain_name_field = loaded_params['chain_name_field']
    
    # Get a dictionary (resname_resnum) for the original receptor
    receptor_sequence_dict = general_functions.get_pdb_sequence_dict(receptor_filename,resname_field,resnumber_field,chain_name_field)
        
    # Get a dictionary (resname_resnum_vales) for the data in the decomp file
    decomp_file_dict = parse_mmgbsa_output(decomp_file)
    
    # Check if both dictionaries match. If NOT, stop execution
    if len(receptor_sequence_dict) != len(decomp_file_dict):
        print("Dictionaries for renumbering do NOT match. Stopping...")
        sys.exit()

    ## Process both dictionaries and output the final .csv file
    # This will out a .csv file with components and renumbered residues values
    decomp_csv_file_renum = prepare_mmgbsa_renumbered_output(receptor_sequence_dict,decomp_file_dict,output_path)
    # This will out a .csv file with components and amber residues values
    decomp_csv_file = prepare_mmgbsa_original_output(receptor_sequence_dict,decomp_file_dict,output_path)

    return decomp_csv_file,decomp_csv_file_renum

def renumber_mmgbsa_output2(output_path,decomp_file,tleap_vs_cristal_resnames_dict,iteration):
    
    # Get a dictionary (resname_resnum_vales) for the data in the decomp file
    decomp_file_dict = parse_mmgbsa_output(decomp_file)
    
    # Check if both dictionaries match. If NOT, stop execution
    if len(tleap_vs_cristal_resnames_dict) != len(decomp_file_dict):
        print("Dictionaries for renumbering do NOT match. Stopping...")
        sys.exit()

    
    print(tleap_vs_cristal_resnames_dict)
    print("#"*10)
    print("#"*10)
    print(decomp_file_dict)
    

    

def parse_mmgbsa_output(filename):
    with open(filename,'r') as input_file:
        file_resname_resnum_dict = {}
        for line in input_file:
            line_split = line.rsplit()
            if len(line_split) == 30 and line_split[3] == "R":
                resname = line_split[0]
                resnumber = line_split[1]
                vdw_comp = line_split[11]
                ele_comp = line_split[15]
                gas_comp = str(round(float(vdw_comp) + float(ele_comp),3))
                pol_solv_comp = line_split[19]
                npol_solv_comp = line_split[23]
                total_comp = line_split[27]
                # Create a dictionary registry value containing all relevant information
                file_resname_resnum_dict[line_split[1]] = resname + '_' + resnumber + '_' + vdw_comp + '_' + ele_comp + '_' + gas_comp + '_' + pol_solv_comp + '_' + npol_solv_comp + '_' + total_comp
                
    input_file.close()
    
    return file_resname_resnum_dict

def prepare_mmgbsa_renumbered_output(dict1,dict2,assay_folder):
    # This will write a .csv file with the components with renumbered residue numbers
    renumbered_file = f'{assay_folder}/mmgbsa_DECOMP_renum.csv'
    with open(renumbered_file,'w') as output_file_renum:
        counter = 1
        for key, value in dict2.items():
            resname,chain,resnum = dict1[int(key)].split('_')
            vdw,ele,gas,pol_solv,np_solv,total = value.split('_')[2:]
            # Write a record into the output .csv file
            record = f"{resname},{chain},{resnum},{vdw},{ele},{gas},{pol_solv},{np_solv},{total}\n"
            # Write the header for the first record
            if counter == 1:
                output_file_renum.write("resname,chain,resnumber,vdw,ele,gas,pol_solv,np_solv,total\n")
                counter += 1

            output_file_renum.write(record)

    return renumbered_file

def prepare_mmgbsa_original_output(dict1,dict2,output_path):
    # This will write a .csv file with the components with renumbered residue numbers
    original_file = f'{output_path}/mmgbsa_DECOMP.csv'
    with open(original_file,'w') as output_file:
        counter = 1
        for key, value in dict2.items():
            resname,chain,resnum = dict1[int(key)].split('_')
            vdw,ele,gas,pol_solv,np_solv,total = value.split('_')[2:]
            # Write a record into the output .csv file
            record = f"{resname},{chain},{key},{vdw},{ele},{gas},{pol_solv},{np_solv},{total}\n"
            # Write the header for the first record
            if counter == 1:
                output_file.write("resname,chain,resnumber,vdw,ele,gas,pol_solv,np_solv,total\n")
                counter += 1

            output_file.write(record)

    return original_file

def create_md_assay_registry(db,docking_assay_id,docking_pose_id):
    # Try to create the 'dynamics_registries' table if it does not exist
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    description = input("Provide a brief description of the MD assay: ")
    
    # Try to create the registers table in case it does no exists
    cursor.execute("CREATE TABLE IF NOT EXISTS dynamics_registries (id INTEGER PRIMARY KEY AUTOINCREMENT, docking_assay_id INTEGER, docking_pose_id INTEGER, description TEXT)")
    conn.commit()
    
    # Store the MD register in the database
    cursor.execute(f"INSERT INTO dynamics_registries (docking_assay_id, docking_pose_id, description) VALUES (?,?,?)",(docking_assay_id,docking_pose_id,description))
    conn.commit()
    
    # Gest the id of the newly create register
    cursor.execute("SELECT id FROM dynamics_registries ORDER BY id DESC LIMIT 1")
    last_assay_id = cursor.fetchone()[0]
    conn.close()
    
    return last_assay_id

def create_md_assay_folder(assay_folder):
    os.makedirs(f"{assay_folder}",exist_ok=True)
    
def parse_receptor_fields(receptor_filename, main_fingerprints_folder):
    with open(receptor_filename,'r') as file:
        residue = 1
        for line in file:
            # In the first line require parsing parameters:
            line_split = line.rsplit()
            if residue == 1 and len(line_split) > 5 and line_split[0] == 'ATOM':
                print("Showing a receptor sample line to input required parameters: \n")
                print(line)
                resname_field = int(input("indicate the 0-based field number of the residue name: "))
                resnumber_field = int(input("indicate the 0-based field number of the residue number: "))
                chain_name_field = int(input("indicate the 0-based field number of the Chain Name (-1 if not present): "))
                
                residue += 1
    
    ## Write the parameters to a file in the assay folder
    params = {"resname_field": resname_field, "resnumber_field": resnumber_field, "chain_name_field": chain_name_field}
    
    with open(f"{main_fingerprints_folder}/receptor_fields.json",'w') as file:
        json.dump(params, file, indent=4)
    
    file.close()
    
    return resname_field,resnumber_field,chain_name_field

def prepare_MD_folder_for_MMGBSA(assay_folder,ligname,amberhome):
    # Prepare the corresponding .prmtop and .inpcrd files for MMGBSA
    apply_ante_MMPBSA(assay_folder,ligname,amberhome)
    write_mmgbsa_input(assay_folder)
    write_MMGBSA_computation_script(assay_folder)
    
def write_MMGBSA_computation_script(assay_folder,traj_input='prod_strip.nc'):
    
    # Get the MMPBSA.py path
    MMPBSA_path = shutil.which('MMPBSA.py')
    
    # Generate the command to run MMPBSA
    command1 = f'{MMPBSA_path} -i {assay_folder}/mmgbsa.in -o {assay_folder}/mmgbsa.out -do {assay_folder}/mmgbsa_DECOMP.out -cp {assay_folder}/complex_MMGBSA.prmtop -rp {assay_folder}/receptor_MMGBSA.prmtop -lp {assay_folder}/ligand_MMGBSA.prmtop -y {assay_folder}/{traj_input}'
    
    # Generate the command to clean MMPBSA
    command2 = f'{MMPBSA_path} -i {assay_folder}/mmgbsa.in -o {assay_folder}/mmgbsa.out -do {assay_folder}/mmgbsa_DECOMP.out -cp {assay_folder}/complex_MMGBSA.prmtop -rp {assay_folder}/receptor_MMGBSA.prmtop -lp {assay_folder}/ligand_MMGBSA.prmtop -y {assay_folder}/{traj_input} --clean'
    
    with open(f'{assay_folder}/mmgbsa_execution.sh','w') as exec_file:
        exec_file.write("#!/bin/bash \n")
        exec_file.write(f"{command1}\n")
        exec_file.write(f"{command2}\n")
    
    exec_file.close()
    
    print(f"Finished preparing MMGBSA execution script in: \n \t {assay_folder}")
    
def print_mmgbsa_info(assay_folder):
    
    print("Assay folder: ", assay_folder)
    
    with open(f"{assay_folder}/mmgbsa.out",'r') as mmgbsa_file:
        activate = 0
        for line in mmgbsa_file:
            if "Differences (Complex - Receptor - Ligand):" in line:
                activate = 1
    
            elif activate == 1:
                line_split = line.split()
                if "VDWAALS" in line:
                    print(f"VDWAALS: {line_split[1]}")
                elif "EEL" in line:
                    print(f"EEL: {line_split[1]}")
                elif "EGB" in line:
                    print(f"EGB: {line_split[1]}")
                elif "ESURF" in line:
                    print(f"ESURF: {line_split[1]}")
                elif "DELTA G gas" in line:
                    print(f"DELTA G gas: {line_split[3]}")
                elif "DELTA G solv" in line:
                    print(f"DELTA G solv: {line_split[3]}")
                elif "DELTA TOTAL" in line:
                    print(f"DELTA G total: {line_split[2]}")
                
                