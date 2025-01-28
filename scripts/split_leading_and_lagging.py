
import argparse
import glob
import pandas as pd
import re

def get_options():
    parser = argparse.ArgumentParser(description='Orient plasmid by leading region then split into first and last X kb (assumes circular plasmid!')    
    parser.add_argument('--output_dir', help='output directory', type=str)
    parser.add_argument('--fasta_file_search', help='string to find files of plasmids with glob (can contain wildcards e.g. 06_plaspan/*/_fasta/*)', type=str)
    parser.add_argument('--leading_region_file', help='file with locations of leading region', type=str)
    parser.add_argument('--leading_region_size', help='leading region size (bp)', default=5000, required=False)
    return parser.parse_args()

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def get_plasmid_seq(input_fasta):
    dna = ''
    with open(input_fasta, 'r') as f:
        for line in f.readlines():
            if not line.startswith('>'):
                dna += line.strip('\n')
            else:
                pass
    return ''.join(char if char in 'ATGC' else 'N' for char in dna)

def orientate_by_leading_region(sequence, start_position=1, orientation="downstream"):
	return_sequence = ''
	if orientation=="downstream":
		return_sequence += sequence[(start_position-1):]
		return_sequence +=  sequence[0:(start_position-1)]
	elif orientation=="upstream":
		return_sequence += reverse_complement(sequence[0:(start_position-1)])
		return_sequence += reverse_complement(sequence[(start_position-1):])
	return return_sequence

def seq_with_leading_region_at_start(seq, 
								start_position,
								orientation):
	plasmid_seq = seq
	new_seq = orientate_by_leading_region(plasmid_seq, 
		int(start_position), orientation)
	return(new_seq)


def leading_and_lagging(plasmid_fasta, start_position, orientation, size=10000):
	plasmid_seq = get_plasmid_seq(plasmid_fasta)
	new_seq = seq_with_leading_region_at_start(plasmid_seq, start_position, orientation)
	leading = new_seq[:size]
	lagging = new_seq[-size:]
	return [leading, lagging]

if __name__=="__main__":
	args = get_options()
	fasta_files = glob.glob(args.fasta_file_search)
	leading_region_df = pd.read_csv(args.leading_region_file, index_col=0, sep='\t')
	size_region = int(args.leading_region_size)
	output_dir = args.output_dir
	for fasta_file in fasta_files:
		plasmid_name = re.sub(".fasta", "", re.sub(".*\\/", "", fasta_file))
		ptu_name = re.sub(".*06_plaspan\\/", "", re.sub("\\/_fasta.*", "", fasta_file))
		if plasmid_name in leading_region_df.index:
			print(plasmid_name, ptu_name)
			leading_region_string = leading_region_df.loc[plasmid_name]["Leading Region"]
			leading_region_start = re.sub(" in.*", "", re.sub("From ", "", leading_region_string))
			#print(leading_region_start)
			leading_region_orientation = re.sub(".*in ", "", re.sub(" direction.", "", leading_region_string))
			new_seqs = leading_and_lagging(fasta_file, leading_region_start, leading_region_orientation, size=size_region)
			with open(output_dir+'/'+plasmid_name+'_leading_'+str(size_region)+'.fa', 'w') as f:
				f.write('>test\n%s' % new_seqs[0])
			with open(output_dir+'/'+plasmid_name+'_lagging_'+str(size_region)+'.fa', 'w') as f:
				f.write('>test\n%s' % new_seqs[1])


