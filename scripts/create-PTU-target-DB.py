import csv
import os


def read_db_to_dict(file_path, key_column, value_column, delimiter=','):
    """
    Reads a file and returns list of unique values in value_column for each unique value of key_column
    Args:
        file_path (str): Path to the input file.
        key_column (str): Name of the column to use as the key.
        value_column (str): Name of the column to use as the value.
        delimiter (str): Delimiter for the file (default is tab-delimited).
        
    Returns:
        dict: A dictionary where keys are unique values from key_column and
              values are lists of unique items from value_column.
    """
    summary_dict = {}
    
    # Open and read the file
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter=delimiter)
        
        # Populate the dictionary
        for row in reader:
            key = row[key_column]
            value = row[value_column]
            if key not in summary_dict:
                summary_dict[key] = []
            if value not in summary_dict[key]:  # Avoid duplicates
                summary_dict[key].append(value)
    
    return summary_dict

def combine_target_dicts(targets_by_genera, genera_by_PTU):
    """
    Combines the dictionaries to create a new dictionary of targets by PTU.
    Args:
        targets_by_genera (dict): Dictionary with {genus: [list_of_targets]}.
        genera_by_PTU (dict): Dictionary with {PTU: [list_of_genera]}.
    
    Returns:
        dict: A dictionary with {PTU: [list_of_targets]}.
    """
    targets_by_PTU = {}
    
    for ptu, genera in genera_by_PTU.items():
        # Initialize an empty list for the current PTU
        targets_by_PTU[ptu] = []
        
        # Collect all targets for the genera associated with the current PTU
        for genus in genera:
            if genus in targets_by_genera:
                targets_by_PTU[ptu].extend(targets_by_genera[genus])
        
        # Remove duplicates from the list of targets
        targets_by_PTU[ptu] = list(set(targets_by_PTU[ptu]))
    
    return targets_by_PTU



def write_target_files(targets_by_PTU, output_directory):
    """
    Writes a file for each PTU with its list of targets.
    
    Args:
        targets_by_PTU (dict): Dictionary with {PTU: [list_of_targets]}.
        output_directory (str): Directory where the files will be written.
    """
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)
    
    for ptu, targets in targets_by_PTU.items():
        file_name = f"{ptu}_targets.txt"
        file_path = os.path.join(output_directory, file_name)
        # Write the targets to the file
        with open(file_path, 'w') as f:
            f.write("\n".join(targets))
    return

# Read in target db
targets_by_genera = read_db_to_dict('data/RM_target_db_deambiguated.csv' , 'genus', 'deambiguated_sequence')

# PTU database
genera_by_PTU = read_db_to_dict('data/plasmid_table.tsv', 'PTU', 'Genus', delimiter='\t')

# Now combine them 
targets_by_PTU = combine_target_dicts(targets_by_genera, genera_by_PTU)

# Write target files
write_target_files(targets_by_PTU, "data/targets-by-PTU/")

