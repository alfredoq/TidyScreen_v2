import sys

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
                    column_widths = [7, 4, 4, 7, 8, 8, 7, 8]
                    # Format the line in accordance to .pdbqt 
                    formated_line = line_split[0].ljust(column_widths[0]) + line_split[1].rjust(column_widths[1]) + '  ' + line_split[2].ljust(column_widths[2]) + line_split[3].ljust(column_widths[3]) + destination_resnumber.ljust(column_widths[4]) + line_split[5].ljust(column_widths[5]) + line_split[6].ljust(column_widths[6]) + line_split[7].ljust(column_widths[7])
                    # Write the formatted line to the output file
                    output_file.write(f"{formated_line} \n")
                
                if len(line_split[2]) == 4:
                    # Parse a new line be written to the renamed .pdbqt file
                    column_widths = [7, 4, 5, 7, 8, 8, 7, 8]
                    # Format the line in accordance to .pdbqt 
                    formated_line = line_split[0].ljust(column_widths[0]) + line_split[1].rjust(column_widths[1]) + ' ' + line_split[2].ljust(column_widths[2]) + line_split[3].ljust(column_widths[3]) + destination_resnumber.ljust(column_widths[4]) + line_split[5].ljust(column_widths[5]) + line_split[6].ljust(column_widths[6]) + line_split[7].ljust(column_widths[7])
                    # Write the formatted line to the output file
                    output_file.write(f"{formated_line} \n")
                    
                if line_split[0] == "TER":
                    output_file.write(line)
                            
    