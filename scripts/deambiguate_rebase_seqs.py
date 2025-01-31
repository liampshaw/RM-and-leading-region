# Deambiguate REBASE sequences
# Created all_rebase_targets_250131 from:
# 1. Downloading https://rebase.neb.com/rebase/link_orgref
# 2. Extracted targets with grep "<5>""
# But they need cleaning and deambiguating, which this script does
import re 
import deambiguate_target_db as de

def clean_sequence(seq):
    return re.sub(r'[\(\)\d/^,\s?]', '', seq)


rebase_sequences = [clean_sequence(line.strip('\n')) for line in open('~/Downloads/trieste//all_rebase_targets_250131.txt')]

unambiguous_sequences = []
for s in rebase_sequences:
	if len(s)<7:
		new_seqs = de.possibleSequences(s)
		for new_s in new_seqs:
			unambiguous_sequences.append(new_s)

with open("../data/all_rebase_targets_250131_6bp_less_deambiguated.txt", "w") as f:
	for s in unambiguous_sequences:
		f.write("%s %d\n" % (s, len(s)))


