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


def prepare_md_initial_files(assay_folder,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm,solvent,min_steps,dynamics=1):
    ## Prepare the corresponding 'tleap' input file
    # These preparations is common for fingerprints and MD assays
    
    if solvent == "explicit":
        prepare_tleap_input(assay_folder,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm)
        prepare_min1_input(assay_folder,min_steps)
        prepare_min2_input(assay_folder,min_steps)
    
    if solvent == "implicit":
        prepare_tleap_input_implicit_solvent(assay_folder,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm)
        prepare_min1_input_implicit_solvent(assay_folder,min_steps)

    if dynamics == 1: # This will only be executed for simulations involving MD preparations
        prepare_heat1_input(assay_folder)
        prepare_heat2_input(assay_folder)
        prepare_equil_input(assay_folder)
        prepare_prod_input(assay_folder)
        prepare_md_execution_script(assay_folder)
    
def prepare_tleap_input(assay_folder,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm):
    with open(f"{assay_folder}/tleap.in","w") as tleap_input:
        tleap_input.write("source leaprc.protein.ff14SB \n")
        tleap_input.write("source leaprc.gaff \n")
        tleap_input.write("source leaprc.water.tip3p \n")
        tleap_input.write(f"UNL = loadmol2 {assay_folder}/{mol2_lig_parm} \n")
        tleap_input.write(f"loadamberparams {assay_folder}/{frcmod_lig_parm} \n")
        tleap_input.write(f"COM = loadpdb {complex_pdb_file} \n")
        tleap_input.write(f"addions COM Na+ 0 \n")
        tleap_input.write(f"addions COM Cl- 0 \n")
        tleap_input.write("solvateoct COM TIP3PBOX 10.0 \n")
        tleap_input.write(f"saveamberparm COM {assay_folder}/complex.prmtop {assay_folder}/complex.inpcrd \n")
        tleap_input.write("quit")
    
    tleap_input.close()

def prepare_tleap_input_implicit_solvent(assay_folder,complex_pdb_file,mol2_lig_parm,frcmod_lig_parm):

    with open(f"{assay_folder}/tleap.in","w") as tleap_input:
        tleap_input.write("source leaprc.protein.ff14SB \n")
        tleap_input.write("source leaprc.gaff \n")
        tleap_input.write("set default PBRadii bondi \n")
        tleap_input.write(f"UNL = loadmol2 {assay_folder}/{mol2_lig_parm} \n")
        tleap_input.write(f"loadamberparams {assay_folder}/{frcmod_lig_parm} \n")
        tleap_input.write(f"COM = loadpdb {complex_pdb_file} \n")
        tleap_input.write(f"saveamberparm COM {assay_folder}/complex.prmtop {assay_folder}/complex.inpcrd \n")
        tleap_input.write("quit")
    
    tleap_input.close()

def prepare_min1_input(assay_folder,min_steps):
    with open(f"{assay_folder}/min1.in","w") as min1_input:    
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

def prepare_min1_input_implicit_solvent(assay_folder,min_steps):

    with open(f"{assay_folder}/min1.in","w") as min1_input:    
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

def prepare_min2_input(assay_folder,min_steps):
    with open(f"{assay_folder}/min2.in","w") as min2_input:
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

def prepare_heat1_input(assay_folder):
    
    with open(f"{assay_folder}/heat1.in","w") as heat1_input:
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
    
def prepare_heat2_input(assay_folder):
   
    with open(f"{assay_folder}/heat2.in","w") as heat2_input:
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
   
def prepare_equil_input(assay_folder):
   
    with open(f"{assay_folder}/equi.in","w") as equi_input:
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
    
def prepare_prod_input(assay_folder):
    
    with open(f"{assay_folder}/prod.in","w") as prod_input:
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
        
def prepare_md_execution_script(assay_folder):
    pmemd_cuda_path = shutil.which('pmemd.cuda')
    with open(f'{assay_folder}/md_execution.sh','a') as exec_file:
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

def prepare_md_execution_script_mendieta(assay_folder,mendieta_assays_path="NONE"):
    with open(f'{assay_folder}/md_execution_mendieta.sh','a') as exec_file:
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

def run_tleap_input(assay_folder,input_file):
    
    tleap_path = shutil.which('tleap')
    command = f'{tleap_path} -f {assay_folder}/{input_file}' 
    subprocess.run(command, shell=True, capture_output=True, text=True)
    
def perform_minimization(assay_folder,ligand_prefix,solvent,inform=1):
    # Compute mol2 with Sybyl Atom Types - for compatibility with RDKit and Meeko
    pmemd_path = shutil.which('pmemd.cuda')
    # Run first minimization procedure
    try:
        command = f'{pmemd_path} -O -i {assay_folder}/min1.in -o {assay_folder}/min1.out -p {assay_folder}/complex.prmtop -c {assay_folder}/complex.inpcrd -r {assay_folder}/min1.crd -ref {assay_folder}/complex.inpcrd' 
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
    
def write_mmgbsa_input(assay_folder):
    with open(f'{assay_folder}/mmgbsa.in','w') as mmgbsa_file:
        mmgbsa_file.write("Per-residue GB and PB decomposition \n")
        mmgbsa_file.write("&general \n")
        mmgbsa_file.write("startframe = 1, interval = 100, verbose=1, \n")
        mmgbsa_file.write("/ \n")
        mmgbsa_file.write("&gb \n")
        mmgbsa_file.write("igb=5, saltcon=0.100, \n")
        mmgbsa_file.write("/ \n")
        mmgbsa_file.write("&decomp \n")
        mmgbsa_file.write("idecomp=2, csv_format=0, \n")
        mmgbsa_file.write("/ \n")
        
    mmgbsa_file.close()

def strip_waters(assay_folder,traj_file,prmtop_file):
    traj_file_prefix = traj_file.split('.')[0]
    traj_file_extension = traj_file.split('.')[1]
    with open(f'{assay_folder}/cpptraj_strip_wat.in','w') as cpptraj_file:
        cpptraj_file.write(f"parm {assay_folder}/{prmtop_file}\n")
        cpptraj_file.write(f"trajin {assay_folder}/{traj_file}\n")
        cpptraj_file.write(f"strip :WAT,Na+,Cl- \n")
        cpptraj_file.write(f"trajout {assay_folder}/{traj_file_prefix}_strip.{traj_file_extension} trajectory  \n")
        cpptraj_file.write(f"go \n")
        cpptraj_file.write(f"quit \n")
    
    cpptraj_file.close()
    
    # Determine the ante_MMPBSA path
    cpptraj_path = shutil.which('cpptraj')
    # Run cpptraj computation to strip water molecules
    command = f'{cpptraj_path} -i {assay_folder}/cpptraj_strip_wat.in' 
    subprocess.run(command, shell=True, capture_output=True, text=True)
    
def apply_ante_MMPBSA(assay_folder):
    # Determine the ante_MMPBSA path
    anteMMPBSA_path = shutil.which('ante-MMPBSA.py')
    # Run ante_MMPBSA computation
    command = f'{anteMMPBSA_path} -p {assay_folder}/complex.prmtop -s :WAT,Na+,Cl- -c {assay_folder}/complex_MMGBSA.prmtop -r {assay_folder}/receptor_MMGBSA.prmtop -l {assay_folder}/ligand_MMGBSA.prmtop -n :UNL' 
    print(f"Computing ante_MMPBSA")
    subprocess.run(command, shell=True, capture_output=True, text=True)

def compute_MMPBSA(assay_folder,traj_input):
    MMPBSA_path = shutil.which('MMPBSA.py')
    
    # Run ante_MMPBSA computation
    command = f'{MMPBSA_path} -i {assay_folder}/mmgbsa.in -o {assay_folder}/mmgbsa.out -do {assay_folder}/mmgbsa_DECOMP.out -cp {assay_folder}/complex_MMGBSA.prmtop -rp {assay_folder}/receptor_MMGBSA.prmtop -lp {assay_folder}/ligand_MMGBSA.prmtop -y {assay_folder}/{traj_input}'
    print(f"Computing MMPBSA")
    subprocess.run(command, shell=True, capture_output=True, text=True)
    
    # Clean the temporal files
    command = f'{MMPBSA_path} -i {assay_folder}/mmgbsa.in -o {assay_folder}/mmgbsa.out -do {assay_folder}/mmgbsa_DECOMP.out -cp {assay_folder}/complex_MMGBSA.prmtop -rp {assay_folder}/receptor_MMGBSA.prmtop -lp {assay_folder}/ligand_MMGBSA.prmtop -y {assay_folder}/{traj_input} --clean'
    
    subprocess.run(command, shell=True, capture_output=True, text=True)
    
    return assay_folder, f"{assay_folder}/mmgbsa_DECOMP.out"
    
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

def renumber_mmgbsa_output(assay_folder,decomp_file,receptor_filename):
    # Get a dictionary (resname_resnum) for the original receptor
    receptor_sequence_dict = general_functions.get_pdb_sequence_dict(receptor_filename,resname_field=3,resnumber_field=4)
    
    # Get a dictionary (resname_resnum_vales) for the data in the decomp file
    decomp_file_dict = parse_mmgbsa_output(decomp_file)
    
    # Check if both dictionaries match. If NOT, stop execution
    if len(receptor_sequence_dict) != len(decomp_file_dict):
        print("Dictionaries for renumbering do NOT match. Stopping...")
        sys.exit()

    ## Process both dictionaries and output the final .csv file
    # This will out a .csv file with components and renumbered residues values
    decomp_csv_file_renum = prepare_mmgbsa_renumbered_output(receptor_sequence_dict,decomp_file_dict,assay_folder)
    # This will out a .csv file with components and amber residues values
    decomp_csv_file = prepare_mmgbsa_original_output(receptor_sequence_dict,decomp_file_dict,assay_folder)

    return decomp_csv_file,decomp_csv_file_renum

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
            resname,resnum = dict1[int(key)].split('_')
            vdw,ele,gas,pol_solv,np_solv,total = value.split('_')[2:]
            # Write a record into the output .csv file
            record = f"{resname},{resnum},{vdw},{ele},{gas},{pol_solv},{np_solv},{total}\n"
            # Write the header for the first record
            if counter == 1:
                output_file_renum.write("resname,resnumber,vdw,ele,gas,pol_solv,np_solv,total\n")
                counter += 1

            output_file_renum.write(record)

    return renumbered_file

def prepare_mmgbsa_original_output(dict1,dict2,assay_folder):
    # This will write a .csv file with the components with renumbered residue numbers
    original_file = f'{assay_folder}/mmgbsa_DECOMP.csv'
    with open(original_file,'w') as output_file:
        counter = 1
        for key, value in dict2.items():
            resname,resnum = dict1[int(key)].split('_')
            vdw,ele,gas,pol_solv,np_solv,total = value.split('_')[2:]
            # Write a record into the output .csv file
            record = f"{resname},{key},{vdw},{ele},{gas},{pol_solv},{np_solv},{total}\n"
            # Write the header for the first record
            if counter == 1:
                output_file.write("resname,resnumber,vdw,ele,gas,pol_solv,np_solv,total\n")
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