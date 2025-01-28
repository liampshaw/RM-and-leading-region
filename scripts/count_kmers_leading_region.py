# Count abundance of targets in leading region of plasmid vs. rest of plasmid

# Libraries
import argparse  
import itertools as iter
import glob
import pandas as pd
import re
from math import floor 
import kmer_functions as kf # Useful k-mer functions 


def get_options():
    parser = argparse.ArgumentParser(description='Orient plasmid by leading region i.e. write fastas so that the first base of the sequence is the start of the leading region. Assumes circular plasmid!')
    parser.add_argument('--output', help='output file', type=str)
    parser.add_argument('--fasta_file_search', help='string to find files of plasmids with glob (can contain wildcards e.g. 06_plaspan/*/_fasta/*)', type=str)
    parser.add_argument('--leading_region_file', help='file with locations of leading region', type=str)
    parser.add_argument('--leading_region_size', help='leading region size (bp)', default=5000, required=False)
    parser.add_argument('--sliding_window', help='sliding window size', default=False, required=False)
    parser.add_argument('--step_size', help='step size for sliding window', default=False, required=False)
    parser.add_argument('--k', help='palindrome length', default=False, required=False)
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

def seq_with_leading_region_at_start(plasmid_fasta, 
								start_position,
								orientation):
	plasmid_seq = get_plasmid_seq(plasmid_fasta)
	new_seq = orientate_by_leading_region(plasmid_seq, 
		int(start_position), orientation)
	return(new_seq)

def read_kmers_from_file(file):
	kmers = []
	with open(file, 'r') as f:
		for line in f.readlines():
			kmers.append(line.strip('\n'))
	return kmers

def count_targets(seq, target_kmers, k=6):
	all_kmers = []
	for kmer in iter.product('ACGT', repeat=k):
		all_kmers.append(''.join(kmer))
	lexicographical_indices_targets = [i for i in range(0, 4**k) if all_kmers[i] in target_kmers]
	total_count = 0
	kmer_counts = kf.count_kmers_seq(seq, k)
	for x in lexicographical_indices_targets:
			total_count += kmer_counts[x]
	return(total_count)


def count_with_sliding_window(plasmid_fasta,
							start_position,
								orientation,
								kmer_targets,
								window_size=5000,
								step_size=500, k=6):
	new_seq = seq_with_leading_region_at_start(plasmid_fasta,
		start_position, orientation)
	outputs = []
	for start in range(0, len(new_seq) - window_size + 1, step_size):
		end = start + window_size
		seq_section = new_seq[start:end]
		outputs.append([start, count_targets(seq_section, kmer_targets, k)])
	return outputs

def main():
	args = get_options()
	fasta_files = glob.glob(args.fasta_file_search)
	leading_region_df = pd.read_csv(args.leading_region_file, index_col=0, sep='\t')
	sliding_window_size = int(args.sliding_window)
	step_size = int(args.step_size)
	with open(args.output, 'w') as output_file:
		for f in fasta_files:
			plasmid_name = re.sub(".fasta", "", re.sub(".*\\/", "", f))
			ptu_name = re.sub(".*06_plaspan\\/", "", re.sub("\\/_fasta.*", "", f))
			print(plasmid_name, ptu_name)
			if plasmid_name in leading_region_df.index:
				leading_region_string = leading_region_df.loc[plasmid_name]["Leading Region"]
				leading_region_start = re.sub(" in.*", "", re.sub("From ", "", leading_region_string))
				leading_region_orientation = re.sub(".*in ", "", re.sub(" direction.", "", leading_region_string))
				PTU = leading_region_df.loc[plasmid_name]["PTU"]
				target_kmers = read_kmers_from_file('/Users/Liam/Downloads/trieste/'+PTU+'_targets.txt')
				results = count_with_sliding_window(f,
	leading_region_start,leading_region_orientation, kmer_targets=target_kmers, k=6)
				for r in results:
					output_file.write(",".join([plasmid_name, ptu_name, str(r[0])+","+str(r[1])])+"\n")


if __name__=="__main__":
	main()




