# Restriction-modification targets and the leading region of plasmids

A repository for scripts related to analysis for upcoming paper. 

The workflow is as follows:


We used a dataset of the 50 most prevalent plasmid taxonomic units (PTUs) in RefSeq200: 4,753 plasmids across 37 genera. These are available at `data/plasmid_table.csv`.
Within each PTU, we ran Prokka on each plasmid followed by Roary with a core threshold of 80% to generate pangenomes using scripts available [here](https://github.com/Adalijuanluo/Plasmid_pan). 
The pangenome data is available for download via Zenodo at [here](zenodo/pangenome_results). 
We then ran `scripts/make_core_accessory_fasta_from_roary.py` to generate core (in the paper we call this 'hard shell') and accessory fastas. 
These combined hard shell and accessory fastas for each PTU are available [here](zenodo/hard_shell_accessory_fastas). In these files the headers of the sequences give information on the genes with `seq_ID|gene_family` e.g. `>NC_019083.1_00042|group_17` 

We then used [RMES](https://forgemia.inra.fr/sophie.schbath/rmes) to analyse occurrences of short motifs (k=4,5,6). The wrapper scripts to do this are available in `scripts/cluster` (for running on Oxford BMRC slurm cluster).
The results from this analysis are available [here](zenodo/rmes_results)

We used [rmsFinder](https://github.com/liampshaw/rmsFinder) to identify Type II RM systems in the 37 genera the plasmids were from. We downloaded up to 100 genomes per genus from RefSeq complete genomes (scripts on BMRC cluster, `projects/R-M-trieste`).
The resulting database of target prevalences across the 37 genera is available at `data/RM_target_db.csv`

 
