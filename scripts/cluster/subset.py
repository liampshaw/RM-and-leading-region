# Subset a pangenome fasta (core/accessory) to one representative per family

import random
from collections import defaultdict
import sysi


def parse_fasta(filename):
    """Parse a FASTA file and yield (header, sequence) tuples."""
    with open(filename, "r") as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header, seq = line, []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)

def select_representatives(fasta_file, output_file):
    groups = defaultdict(list)

    # Group sequences by their group ID
    for header, sequence in parse_fasta(fasta_file):
        group_id = header.split("|")[-1]  # Extract group ID after '|'
        groups[group_id].append((header, sequence))

    # Select one random representative per group
    with open(output_file, "w") as f:
        for group_id, entries in groups.items():
            header, sequence = random.choice(entries)
            f.write(f"{header}\n{sequence}\n")

if __name__=="__main__":
    select_representatives(sys.argv[1], sys.argv[2])
