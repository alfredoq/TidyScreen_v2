import warnings
warnings.filterwarnings("ignore")
import sys
import sqlite3
from tidyscreen import tidyscreen as tidyscreen
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
from itertools import islice
import tarfile
import io
from rdkit import Chem
import concurrent.futures
import moldf
import copy

def renumber_pdb_file_using_crystal(crystal_file,target_file,renumbered_file,resname_field=3,resnumber_field=5):
    crystal_dict = get_pdb_sequence_dict(crystal_file,resname_field=3,resnumber_field=5)
    amber_dict = get_pdb_sequence_dict(target_file,resname_field=3,resnumber_field=4)
    combined_dict = combine_dictionaries(amber_dict,crystal_dict)
    renumber_pdb_file(target_file,renumbered_file,combined_dict,resname_field=3,resnumber_field=4)
    
def get_pdb_sequence_dict(pdb_file,resname_field,resnumber_field):
    """
    This function will create a dictionary from parsing the crystallographic .pdb file, returning:

        - key: the 'sequential' numbering of the residue
        - value: a string in the form: "X1_X2" in which X1 is the residue name and X2 is the corresponding numbering
    """
    with open(pdb_file,'r') as file:
        sequence_dict = {}
        residue = 1
        for line in file:
            line_split = line.rsplit()
            if len(line_split) > 5 and line_split[0] == 'ATOM':
                residue_value = f"{line_split[resname_field]}_{line_split[resnumber_field]}"
                if residue_value not in sequence_dict.values():
                    sequence_dict[residue] = residue_value
                    residue += 1
    
    file.close()
    
    return sequence_dict

def combine_dictionaries(dict1,dict2):
    """
    Will return a dictionry in which:
    
        - keys: 'values' from dict1
        - values: 'values' from dict2
    """
    # Check if both dictionaries are of the same length. If not, inform and exit.
    if len(dict1) == len(dict2):
        combined_dict = {dict1[k]: dict2[k] for k in dict1 if k in dict2}
        return combined_dict
    else:
        print("The lenght of the dictionaries corresponding to residue names are not equal. Stopping...")
        sys.exit()
    
def renumber_pdb_file(target_file,renumbered_file,combined_dict,resname_field,resnumber_field):
    with open(target_file,'r') as input_file, open(renumbered_file,'w') as output_file:
        for line in input_file:
            line_split = line.rsplit()
            if len(line_split) > 5 and line_split[0] == 'ATOM':
                reference_key = f"{line_split[resname_field]}_{line_split[resnumber_field]}"
                destination_key_value = combined_dict[reference_key]
                destination_resnumber = destination_key_value.split('_')[1]
                
                if len(line_split[2]) <= 3:
                    # Parse a new line be written to the renamed .pdbqt file
                    column_widths = [7, 4, 4, 6, 3, 12, 8, 8]
                    # Format the line in accordance to .pdbqt 
                    formated_line = line_split[0].ljust(column_widths[0]) + line_split[1].rjust(column_widths[1]) + '  ' + line_split[2].ljust(column_widths[2]) + line_split[3].ljust(column_widths[3]) + destination_resnumber.rjust(column_widths[4]) + line_split[5].rjust(column_widths[5]) + line_split[6].rjust(column_widths[6]) + line_split[7].rjust(column_widths[7])
                    # Write the formatted line to the output file
                    output_file.write(f"{formated_line} \n")
                
                if len(line_split[2]) == 4:
                    # Parse a new line be written to the renamed .pdbqt file
                    column_widths = [7, 4, 5, 6, 3, 12, 8, 8]
                    # Format the line in accordance to .pdbqt 
                    formated_line = line_split[0].ljust(column_widths[0]) + line_split[1].rjust(column_widths[1]) + ' ' + line_split[2].ljust(column_widths[2]) + line_split[3].ljust(column_widths[3]) + destination_resnumber.rjust(column_widths[4]) + line_split[5].rjust(column_widths[5]) + line_split[6].rjust(column_widths[6]) + line_split[7].rjust(column_widths[7])
                    # Write the formatted line to the output file
                    output_file.write(f"{formated_line} \n")
                    
                if line_split[0] == "TER":
                    output_file.write(line)

def subset_table(db,source_table,dest_table,colname,filter):
    conn = tidyscreen.connect_to_db(db) # Connect to the main database
    cursor = conn.cursor()
    
    # Build a WHERE clause dynamically to use all columns
    where_clause = ' AND '.join([f"{col} = ?" for col in colname])

    # Construct the full query
    query = f"""CREATE TABLE IF NOT EXISTS {dest_table}_temp AS
                SELECT * FROM {source_table} WHERE {where_clause};"""

    # Execute with parameters
    cursor.execute(query, filter)

    # Reset the index of temp table
    reset_id_in_column(conn,f"{dest_table}_temp",dest_table)

def reset_id_in_column(conn,source_table,target_table):

    # Get all column names except 'id'
    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info({target_table})")
    columns = [row[1] for row in cursor.fetchall() if row[1] != 'id']
    col_str = ', '.join(columns)

    # Drop the target table if it exists
    cursor.execute(f"DROP TABLE IF EXISTS {target_table}")

    # Recreate the target table with same schema, resetting the id
    # This assumes the schema is known or you can adapt it as needed
    # Here's an example assuming same schema except for id
    # Use PRAGMA or manual schema copy logic here if needed
    cursor.execute(f"""
        CREATE TABLE {target_table} (
            'id' INTEGER PRIMARY KEY AUTOINCREMENT,
            {', '.join(f"{col} TEXT" for col in columns)}  -- adjust types if needed
        )
    """)

    # Insert rows without id
    cursor.execute(f"""
        INSERT INTO {target_table} ({col_str})
        SELECT {col_str} FROM {source_table}
    """)

    # Drop the temporary table after reindexing
    cursor.execute(f"DROP TABLE IF EXISTS {source_table}")


    conn.commit()
    conn.close()

def create_prolif_reference_df(tleap_vs_cristal_reference_dict,interactions_list):
    # Generate the column names by combining each dict value with each list item (the interactions to be computed)
    colnames = [f"{val}_{item}" for val in tleap_vs_cristal_reference_dict.values() for item in interactions_list]
    all_residues_plf_df = pd.DataFrame(columns=colnames)
    
    return all_residues_plf_df

def map_prolif_fingerprints_df_to_crystal_sequence(df,tleap_vs_cristal_reference_dict):
            
    # Drop the first level of the 'df' (corresponds to the Ligand name)
    df.columns = df.columns.droplevel('ligand')  # or use the level number
    
    # Rename upper level of the multindex column ProLIF fingerprints df
    df.columns = pd.MultiIndex.from_tuples([(tleap_vs_cristal_reference_dict.get(upper, upper), lower) for upper, lower in df.columns]
)
    # Combine names in the two levels of the fingerprints dataframe
    df.columns = ['{}_{}'.format(upper, lower) for upper, lower in df.columns]
    
    return df

def merge_calculated_and_reference_fingerprints_df(calc_df,reference_df):
    
    merged_df = pd.concat([reference_df, calc_df], ignore_index=True, sort=False).fillna(0)
    merged_df = merged_df.replace({'True':1}).astype(int)
    
    return merged_df

def generate_tar_file(file):
    """Compress the pdb file into a TAR archive and return its binary content."""
    filename = file.split("/")[-1]
    tar_buffer = io.BytesIO()  # Create an in-memory buffer
    with tarfile.open(fileobj=tar_buffer, mode="w") as tar:
        tar.add(file, arcname=filename)  # Store only filename in TAR
    
    return tar_buffer.getvalue()  # Return TAR file as binary

def sort_table(assay_folder,assay_id,table,column):
    results_db = f"{assay_folder}/assay_{assay_id}.db"
    conn = tidyscreen.connect_to_db(results_db)
    cursor = conn.cursor()
    
    # Replace 'your_table' and 'sort_column' with your table and column names
    cursor.execute(f"""
        CREATE TABLE sorted_table AS
        SELECT * FROM {table}
        ORDER BY {column};
    """)
    
    # Drop the original table
    cursor.execute(f"DROP TABLE {table};")
    
    # Rename the sorted table to original name
    cursor.execute(f"ALTER TABLE sorted_table RENAME TO {table};")
    
    conn.commit()
    conn.close()

def csv_reader(file):
    """
    Will read a .csv file and return a pandas dataframe
    """
    # Generate the target table name while also assuring compatibility SQL standards in the name
    
    target_table_name = file.split("/")[-1].replace(".csv","").replace(".smi","").replace("-", "_") 
    
    df = pd.read_csv(file,header=None,index_col=False)
    #df = pd.read_csv(file,header=None)
    #df = df.reset_index()
    
    # First row, second column (i.e. the first SMILES)
    first_element = df.iloc[0, 0]  
    
    return target_table_name, df, first_element

def check_table_presence(conn,table_name):
    "Will return 1 if exists, otherwise returns 0"
    cursor = conn.cursor()
    cursor.execute("""SELECT name FROM sqlite_master WHERE type='table' AND name=?;""", (table_name,))

    table_exists = cursor.fetchone()

    if table_exists:
        print(f"Table {table_name} exists.")
        return 1
    else:
        #print(f"Table '{table_name}' DOES NOT exists.")
        return 0

def save_df_to_db(db,df,table_name):
    conn = tidyscreen.connect_to_db(db)
    exists = check_table_presence(conn,table_name)

    if exists == 1:
        replace_table_action = input(f"The table named '{table_name}' already exists in the database. I will replace it, are you ok with that? (y/n): ")

        if replace_table_action == 'y':
            print(f"I will overwrite table '{table_name}' as indicated.'")
        else:
            print(f"Quiting to safely retain table {table_name}. Stopping...")
            sys.exit()

    try:
        df.to_sql(con=conn, name=table_name,if_exists="replace",index=None)
        # Inform processing
        print(f"Table '{table_name}' created in: '{db}'")
    
    except Exception as error:
        print(error)
        
def check_smiles(smiles):
    
    try:     
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            print(f"Problem reading SMILES columns - Example {smiles} \n Stopping...")
            sys.exit()
        else:
            print("SMILES column valid...")
            
    except:
        print(f"Problem reading SMILES columns - Example {smiles} \n Stopping...")
        sys.exit()
        
def write_failed_smiles_to_db(smiles,db,file,cause):
    conn = tidyscreen.connect_to_db(db)
    cursor = conn.cursor()
    
    # Create a table to store failed SMILES
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS failed_smiles (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            smiles TEXT, 
            source TEXT,
            cause TEXT
        )
    """)
    # Insert the failed SMILES into the table
    cursor.execute("INSERT INTO failed_smiles (smiles, source, cause) VALUES (?,?,?)", (smiles,file,cause))
    conn.commit()
    conn.close()

def delete_nulls_table(db,table_name,column_name):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    cursor.execute(f"DELETE FROM {table_name} WHERE {column_name} IS NULL")
    conn.commit()
    conn.close()

def timeout_function(function, args=(),kwargs=None, timeout=60,on_timeout_args=None):
    if kwargs is None:
        kwargs = {}
    with concurrent.futures.ThreadPoolExecutor() as executor:
    #with concurrent.futures.ProcessPoolExecutor() as executor:
        future = executor.submit(function, *args, **kwargs)
        try:
            return future.result(timeout=timeout)
        except concurrent.futures.TimeoutError:
            write_failed_smiles_to_db(on_timeout_args[0],on_timeout_args[1],on_timeout_args[2],on_timeout_args[3])
            return None
        
def sybyl_to_gaff_mol2(input_mol2, output_mol2, sybyl_to_gaff_dict):
    with open(input_mol2, 'r') as infile, open(output_mol2, 'w') as outfile:
        atom_section = False
        for line in infile:
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                outfile.write(line)
                continue
            elif line.startswith("@<TRIPOS>BOND"):
                atom_section = False
                outfile.write(line)
                continue

            if atom_section and line.strip():
                parts = line.split()
                if len(parts) > 5:
                    sybyl_type = parts[5]
                    gaff_type = sybyl_to_gaff_dict.get(sybyl_type, sybyl_type)
                    parts[5] = gaff_type
                    line = "{:<7} {:<4} {:>10} {:>10} {:>10} {:<6} {:>5} {:<4}\n".format(
                        parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], parts[7]
                    ) if len(parts) >= 8 else " ".join(parts) + "\n"
            outfile.write(line)
            
def sybyl_to_gaff_mol2_moldf(sybyl_file,gaff_file,atoms_dict):
    
    # Read the mol2 file and make a copy for the gaff version
    mol_sybyl = moldf.read_mol2(sybyl_file)
    mol_gaff = copy.deepcopy(mol_sybyl)
    
    # Update atom types
    for idx, row in mol_sybyl['ATOM'].iterrows():
        sybyl_atom = row['atom_type']
        gaff_type = atoms_dict.get(sybyl_atom, sybyl_atom)
        mol_gaff['ATOM'].at[idx,'atom_type'] = gaff_type
        
    # Write the modified mol2 file
    moldf.write_mol2(mol_gaff,gaff_file)

    return gaff_file

def delete_smiles_row_from_table(smiles,db,table_name):
    conn = sqlite3.connect('your_database.db')
    cursor = conn.cursor()

    # Delete rows where column_name matches a specific value
    cursor.execute(f"DELETE FROM {table_name} WHERE SMILES = ?", (smiles,))

    conn.commit()
    conn.close()

def replace_charge_on_mol2_file(mol2_file,charge_array):
    # Read the .mol2 file
    mol = moldf.read_mol2(mol2_file)
    try:
        for idx, row in mol['ATOM'].iterrows():
            mol['ATOM'].loc[idx,['charge']] = charge_array[idx]
            
        # Replace the original file with the new one containing replaced charges
        moldf.write_mol2(mol, mol2_file)
    
        return mol2_file
    
    except Exception as error:
        print(Exception)