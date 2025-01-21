# Convert pangenome files into core and accessory gene files

import argparse  
import pandas as pd
from math import floor
import re

def get_options():
	parser = argparse.ArgumentParser(description='Creates core and accessory gene files.',
                                     prog='make_core_accessory_fasta')
	parser.add_argument('--ffn_dir', help='Directory where .ffn files are stored', required=True)
	parser.add_argument('--gene_presence_absence', help='pan_genome_reference.Rtab (roary output)', required=True) 
	parser.add_argument('--clustered_proteins', help='file of clustered proteins (roary output)', required=True) 
	parser.add_argument('--core', help='threshold for core (default=0.8)', required=False, default=0.8) 
	parser.add_argument('--output', help='output prefix', required=False, default='output') 
	return parser.parse_args()


def read_clustered_proteins(clustered_proteins_file):
	family_dict = {}
	with open(clustered_proteins_file, 'r') as f:
		for line in f.readlines():
			if ":" in line:
				line = line.split(':')
				family_dict[line[0]] = line[1].strip('\n').split()
	return family_dict


def read_pan_genome_reference_fasta(fasta):
	seq_dict = {}
	current_family = ''
	for line in open(fasta, 'r').readlines():
		if line.startswith('>'):
			family_name = line.strip('\n').split(' ')[1]
			seq_dict[family_name] = ''
			current_family = family_name
		else:
			seq = line.strip('\n')
			seq_dict[current_family] += seq
	return(seq_dict)

def read_in_plasmid_genes(plasmid_ids, ffn_dir):
	seq_dict = {}
	for plasmid_id in plasmid_ids:
		for line in open(ffn_dir+'/'+plasmid_id+'.ffn', 'r').readlines():
			if line.startswith('>'):
				seq_id = line.strip('\n').split(' ')[0][1:]
				seq_dict[seq_id] = ''
				current_seq = seq_id
			else:
				seq = line.strip('\n')
				seq_dict[current_seq] += seq
	return seq_dict



def main():
	args = get_options()
	# read in p/a and calculate core threshold
	gene_pa_df = pd.read_csv(args.gene_presence_absence, sep='\t', index_col=0)
	number_of_genomes = len(gene_pa_df.columns)
	core_threshold_n = floor(number_of_genomes*float(args.core))

	# find core/accessory families
	row_sums = gene_pa_df.sum(axis=1).to_dict()
	core_families = [x for x in row_sums.keys() if row_sums[x]>core_threshold_n]
	accessory_families = [x for x in row_sums.keys() if row_sums[x]<=core_threshold_n]

	# read in clustered proteins
	protein_families = read_clustered_proteins(args.clustered_proteins)
	# and get the list of plasmid accessions
	all_proteins = [protein for family in protein_families.values() for protein in family]
	plasmids = list(set(['_'.join(x.split('_')[0:2]) for x in all_proteins]))

	# now read in those plasmid gene sequences
	gene_sequences = read_in_plasmid_genes(plasmids, args.ffn_dir)


	with open(args.output+'_core_genes.fa', 'w') as f:
		for family in core_families: # for each family
			for gene_id in protein_families[family]: # write all the gene sequences within it to file
				f.write('>%s|%s\n%s\n' % (gene_id, family, gene_sequences[gene_id]))
	with open(args.output+'_accessory_genes.fa', 'w') as f:
		for family in accessory_families:
			for gene_id in protein_families[family]:
				f.write('>%s|%s\n%s\n' % (gene_id, family, gene_sequences[gene_id]))


if __name__=="__main__":
	main()