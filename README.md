# Restriction-modification targets and the leading region of plasmids

A repository for scripts related to analysis for upcoming paper. 

The workflow is as follows:


We used a dataset of the 50 most prevalent plasmid taxonomic units (PTUs) in RefSeq200: 4,753 plasmids across 37 genera. These are available at `data/plasmid_table.csv`. The taxonomy for each genera is available at `data/genera_taxonomy.txt`
Within each PTU, we ran Prokka on each plasmid followed by Roary with a core threshold of 80% to generate pangenomes using scripts available [here](https://github.com/Adalijuanluo/Plasmid_pan). 
The pangenome data is available for download via Zenodo at [here](zenodo/pangenome_results). A summary of the pangenome stats is available at `data/PTU-pangenome-stats.csv` (Note that `PTU-B/O/K/Z` is renamed `PTU-BOKZ` and `PTU-L_M` becomes `PTU-LM` from this point forwards) 
We then ran `scripts/make_core_accessory_fasta_from_roary.py` to generate core (in the paper we call this 'hard shell') and accessory fastas. 
These combined hard shell and accessory fastas for each PTU are available [here](zenodo/hard_shell_accessory_fastas). In these files the headers of the sequences give information on the genes with `seq_ID|gene_family` e.g. `>NC_019083.1_00042|group_17` 

We then used [RMES](https://forgemia.inra.fr/sophie.schbath/rmes) to analyse occurrences of short motifs (k=4,5,6). The wrapper scripts to do this are available in `scripts/cluster` (for running on Oxford BMRC slurm cluster).
The results from this analysis are available [here](zenodo/rmes_results)

We used [rmsFinder](https://github.com/liampshaw/rmsFinder) to identify Type II RM systems in the 37 genera the plasmids were from. We downloaded up to 100 genomes per genus from RefSeq complete genomes (scripts on BMRC cluster, `projects/R-M-trieste`).
The resulting database of target prevalences across the 37 genera is available at `data/RM_target_db.csv`. This database contains ambiguity codes; it was de-ambiguated with `python scripts/deambiguate-target-db.py` to produce `data/RM_target_db_deambiguated.csv`.
The k-mers of a given length that are targeted by any RM system in the dataset can then be extracted with e.g.
```
awk -F ',' '$7=="4"' data/RM_target_db_deambiguated.csv | cut -d ',' -f 6 | sort -n  | uniq > 4mer_targets.txt
```

We computed discrepancies in RMES score between hard shell (aka core) vs. accessory with `scripts/rmes-score-discrepancy.py` with a wrapper script `scripts/run-rmes-score-discrepancy.py` that simply runs this on all RMES results files from the different PTUs. 

We created a database of targets by PTU (i.e. for a given PTU, the set of all motifs targeted by RM systems detected within genomes of its hosts) using `python scripts/create-PTU-target-DB.py`. The results are in `data/targets-by-PTU`. 
These can then be used to generate discrepancy scores for targets using e.g.
```
python ~/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/scripts/run-rmes-score-discrepancy.py --script ~/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/scripts/rmes-score-discrepancy.py --dir ./ --kmers per-PTU --output test_RM_targets.csv --k 6
```


# Python dependencies

The python scripts depend on `scipy, pandas`. 
