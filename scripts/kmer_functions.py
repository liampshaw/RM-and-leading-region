import argparse
import itertools as iter  

# For reverse complement
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def get_options():
    parser = argparse.ArgumentParser(description='Counts k-mers in a fasta and returns in lexicographic order.',
                                     prog='countKmers')
    parser.add_argument('--fasta', help='fasta file', required=True)
    parser.add_argument('--k', help='value of k', required=True) 
    return parser.parse_args()


def reverse_complement(seq):
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def all_palindromes(N):
    '''Returns all palindromes of length N.'''
    palindromes = []
    for seq in iter.product('ATCG', repeat=N):
        seq = ''.join(seq)
        half_seq = seq[0:int(N/2)]
        end_seq = seq[int(N/2):int(N)]
        if half_seq==reverse_complement(end_seq):
            palindromes.append(seq)
    return palindromes


def get_seq(input_fasta):
    dna = ''
    with open(input_fasta, 'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                dna += 'Z'
            else:
                dna += line.strip('\n')
    return(dna)

# Using Matthew Fahrbach solution from Rosalind problem 
# K-mer counts are enumerated in lexicographical order
# of the k-mers
# e.g. AAAA comes first for k=4
def count_kmers(fasta, k, with_rev_comp=False):
    dna = get_seq(fasta)
    if with_rev_comp:
        rev_comp_dna = reverse_complement(dna)
        dna = dna+"Z"+rev_comp_dna
    table = {};
    for kmer in iter.product('ACGT', repeat=k):
        table[''.join(kmer)] = 0

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i + k]
        if all([base in ['A', 'T', 'C', 'G'] for base in kmer]): # only if kmer is ATCG-based, no other characters allowed
            table[kmer] = table[kmer] + 1
    output_list = []
    for kmer in iter.product('ACGT', repeat=k):
        output_list.append(table[''.join(kmer)])
        #print(, end=' ')
    return(output_list)

def count_kmers_seq(dna, k, with_rev_comp=False):
    if with_rev_comp:
        rev_comp_dna = reverse_complement(dna)
        dna = dna+"Z"+rev_comp_dna
    table = {};
    for kmer in iter.product('ACGT', repeat=k):
        table[''.join(kmer)] = 0

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i + k]
        if all([base in ['A', 'T', 'C', 'G'] for base in kmer]): # only if kmer is ATCG-based, no other characters allowed
            table[kmer] = table[kmer] + 1
    output_list = []
    for kmer in iter.product('ACGT', repeat=k):
        output_list.append(table[''.join(kmer)])
        #print(, end=' ')
    return(output_list)

def main():
    args = get_options()
    input_fasta = args.fasta
    k = int(args.k)
    kmers = count_kmers(input_fasta, k)  
    kmers = [str(k) for k in kmers] 
    print(' '.join(kmers))


if __name__=="__main__":
	main()