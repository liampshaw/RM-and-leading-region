# Orient plasmid by leading region

# Libraries
import argparse  
from itertools import product
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


def palindromes_in_seq(seq, k):
	kmers = kf.count_kmers_seq(seq, k)
	total_palindromes = 0
	for x in lexicographical_indices_palindromes:
			total_palindromes += kmers[x]
	return(total_palindromes)


def find_dna_palindromes(sequence, k):
    def reverse_complement(s):
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N':'N'}
        return ''.join(complement[base] for base in reversed(s))
    locations = []
    for i in range(len(sequence) - k + 1):
        segment = sequence[i:i+k]
        if segment == reverse_complement(segment):
            locations.append(i)
    return locations


def palindrome_density(sequence, k, window_size, step_size):
    palindrome_locations = find_dna_palindromes(sequence, k)
    densities = []
    for start in range(0, len(sequence) - window_size + 1, step_size):
        end = start + window_size
        count = sum(1 for loc in palindrome_locations if start <= loc < end)
        density = count / window_size
        densities.append((start, 
        	density, 
        	sequence[start:end].count('G')+sequence[start:end].count('C')))
    return densities

def leading_region(plasmid_fasta, 
								start_position,
								orientation):
	plasmid_seq = get_plasmid_seq(plasmid_fasta)
	new_seq = orientate_by_leading_region(plasmid_seq, 
		int(start_position), orientation)
	return(new_seq)

def palindromes_sliding_window(plasmid_fasta,
							start_position,
								orientation,
								window_size=5000,
								step=100, k=6):
	new_seq = leading_region(plasmid_fasta,
		start_position, orientation)
	# outputs = []
	# for x in range(0, len(new_seq)-window_size, step):
	# 	sequence = new_seq[x:(x+window_size)]
	# 	outputs.append([x+window_size/2, palindromes_in_seq(sequence, k)])
	outputs = palindrome_density(new_seq, k, window_size, step_size=step)
	return outputs



def palindromes_leading_region(plasmid_fasta, 
								start_position,
								orientation,
								leading_region_size=5000, k=6):	
	new_seq = leading_region(plasmid_fasta, start_position, orientation)
	# Palindrome info
	palindromes = kf.all_palindromes(k)
	kmers = [''.join(kmer) for kmer in product('ACGT', repeat=k)]
	lexicographical_indices_palindromes = [i for i in range(0, 4**k) if kmers[i] in palindromes]

	# Work out palindromes in first 5kb
	kmers_5kb = kf.count_kmers_seq(new_seq[0:leading_region_size], k)
	total_palindromes_5kb = 0
	for x in lexicographical_indices_palindromes:
			total_palindromes_5kb += kmers_5kb[x]

	kmers_rest = kf.count_kmers_seq(new_seq[leading_region_size:], k)
	total_palindromes_rest = 0
	for x in lexicographical_indices_palindromes:
			total_palindromes_rest += kmers_rest[x]

	#print(leading_region_size)
	#print(leading_region_size, 
	#	len(plasmid_seq)-leading_region_size,
	#	len(new_seq[leading_region_size:])
	if len(new_seq)>leading_region_size:
		return [len(new_seq), leading_region_size, 
					(new_seq[0:leading_region_size].count("G")+new_seq[0:leading_region_size].count("C"))/leading_region_size, 
					(new_seq[leading_region_size:].count("G")+new_seq[leading_region_size:].count("C"))/(len(new_seq)-leading_region_size),
					total_palindromes_5kb/leading_region_size, total_palindromes_rest/len(new_seq[leading_region_size:])]
	else:
		return [len(new_seq), leading_region_size,
		(new_seq[0:leading_region_size].count("G"))/leading_region_size, 
					new_seq[leading_region_size:].count("GC")/(len(new_seq)-leading_region_size), "NA", "NA"]


def main():
	args = get_options()
	fasta_files = glob.glob(args.fasta_file_search)
	leading_region_df = pd.read_csv(args.leading_region_file, index_col=0, sep='\t')
	if args.sliding_window==False:
		print("no sliding window")
		with open(args.output, 'w') as output_file:
			output_file.write("PTU,plasmid,plasmid_length,leading_region_taken,leading_region_GC,rest_of_plasmid_GC,leading_region_density,rest_plasmid_density\n")
			for f in fasta_files:
				plasmid_name = re.sub(".fasta", "", re.sub(".*\\/", "", f))
				ptu_name = re.sub(".*06_plaspan\\/", "", re.sub("\\/_fasta.*", "", f))
				print(plasmid_name)
				if plasmid_name in leading_region_df.index:
					if len(leading_region_df.loc[plasmid_name])!=2: # ignore double
						#print(leading_region_df.loc[plasmid_name])
						leading_region_string = leading_region_df.loc[plasmid_name]["Leading Region"]
						leading_region_start = re.sub(" in.*", "", re.sub("From ", "", leading_region_string))
						#print(leading_region_start)
						leading_region_orientation = re.sub(".*in ", "", re.sub(" direction.", "", leading_region_string))
						#print(leading_region_orientation)
						output = palindromes_leading_region(f, 
							int(leading_region_start),
							leading_region_orientation, int(args.leading_region_size))
						output_file.write("%s,%s,%s\n" % (ptu_name,plasmid_name, ",".join([str(x) for x in output])))
	elif args.sliding_window!=False:
		sliding_window_size = int(args.sliding_window)
		step_size = int(args.step_size)
		print(leading_region_df.index)
		with open(args.output, 'w') as output_file:
			for f in fasta_files:
				plasmid_name = re.sub(".fasta", "", re.sub(".*\\/", "", f))
				ptu_name = re.sub(".*06_plaspan\\/", "", re.sub("\\/_fasta.*", "", f))
				print(plasmid_name, ptu_name)
				if plasmid_name in leading_region_df.index:
					if len(leading_region_df.loc[plasmid_name])!=2: # ignore double oriT
						#print(leading_region_df.loc[plasmid_name])
						leading_region_string = leading_region_df.loc[plasmid_name]["Leading Region"]
						leading_region_start = re.sub(" in.*", "", re.sub("From ", "", leading_region_string))
						#print(leading_region_start)
						leading_region_orientation = re.sub(".*in ", "", re.sub(" direction.", "", leading_region_string))
						outputs = palindromes_sliding_window(f, 
							int(leading_region_start), leading_region_orientation, window_size=sliding_window_size, step=step_size, k=int(args.k))
						for o in outputs:
							output_file.write("%s,%s,%s\n" % (ptu_name,plasmid_name, ",".join([str(x) for x in o])))


if __name__=="__main__":
	main()






