import pandas as pd

ambiguity_codes =   {'A': ['A'],\
                    'G' : ['G'],\
                    'C'	: ['C'],\
                    'T' : ['T'],\
                    'Y'	: ['C', 'T'],\
                    'R'	: ['A','G'],\
                    'W'	: ['A','T'],\
                    'S'	: ['G','C'],\
                    'K'	: ['T','G'],\
                    'M'	: ['C','A'],\
                    'D'	: ['A','G','T'],\
                    'V'	: ['A','C','G'],\
                    'H'	: ['A','C','T'],\
                    'B'	: ['C','G','T'],\
                    'X'	: ['A','C','G','T'],\
                    'N'	: ['A','C','G','T'],\
                    '-' : ['-']}

def possibleSequences(dna_sequence):
    '''Takes a DNA sequence (possibly containing ambiguous bases) and returns
    a list of all possible sequences.
    Args:
        dna_sequence (str)
            String of DNA sequence (uppercase)
    Returns:
        possible_strings (list)
            List of strings - all possible DNA sequences
    '''
    # Get all possible bases at each position
    possible_bases = [ambiguity_codes[base] for base in dna_sequence]
    # Go through and add these on
    possible_strings = []
    for bases in possible_bases:
        if len(possible_strings)==0: # if no strings yet, use first base
            possible_strings = list(bases)
        elif len(bases)==1: # if just one possible base at position, add it on
            possible_strings = [x+bases[0] for x in possible_strings]
        else: # if multiple possibilities, add them on one by one
            additions = []
            for base in bases:
                additions.append([x+base for x in possible_strings])
            possible_strings = [x for addition in additions for x in addition]
    return possible_strings


target_db = pd.read_csv("data/RM_target_db.csv")
new_target_db = []
for i, sequence in enumerate(target_db['sequence']):
    new_sequences = possibleSequences(sequence)
    row = list(target_db.loc[i])
    for s in new_sequences:
        new_target_db.append(row+[s])

new_target_db_df = pd.DataFrame(new_target_db)
new_target_db_df.columns = ['genus', 'sequence', 'genome', 'total_genomes', 'normalized_count', 'deambiguated_sequence']
new_target_db_df = new_target_db_df.assign(length=[len(x) for x in new_target_db_df['deambiguated_sequence']])
new_target_db_df.to_csv("data/RM_target_db_deambiguated.csv", index=False)




