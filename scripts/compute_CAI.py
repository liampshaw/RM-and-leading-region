# Script to compute Codon Adaptation Index for genes across a plasmid

import math
from collections import defaultdict
from Bio.Data import CodonTable
import pandas as pd
from Bio import SeqIO
from collections import Counter
import sys

def codon_usage(seq):
    '''Codon usage for a given sequence.'''
    codons = [str(seq[i:i+3]) for i in range(0, len(seq), 3)
              if len(seq[i:i+3]) == 3]
    return Counter(codons)

def compute_cai(counts, weights):
    '''CAI as defined by Sharp and Li (1987) and used by Lucks et al. (2008)'''
    total_codons = sum(counts.values())
    log_sum = 0
    for codon, count in counts.items():
        if codon in weights and weights[codon] > 0:
            log_sum += count * math.log(weights[codon])
    return math.exp(log_sum / total_codons)

lucks_weights = pd.read_csv('/Users/Liam/Dropbox/_Projects/2023-Trieste/2026-new-data/lucks_ecoli_codon_weights.tsv', sep='\t', comment='#')
weights_dict = dict(zip(lucks_weights["codon"], lucks_weights["w"])) 

for record in SeqIO.parse(sys.argv[1], "fasta"):
    usage = codon_usage(record.seq)
    cai = compute_cai(usage, weights_dict)
    print(record.id, cai)
