# Restriction-modification targets and the leading region of plasmids

A repository of analysis related to the paper:

**The leading region of many conjugative plasmids is depleted in restriction-modification targets** 
Liam P. Shaw, Arancha Peñil Celis, Manuel Ares Arroyo, Lijuan Luo, Tatiana Dimitriu, Fernando de la Cruz 

**Repository currently undergoing streamlining**

## Paper analysis

The main analysis presented in the paper is performed in `scripts/Paper-analysis.Rmd`. Associated python scripts used for data processing before then are in `scripts` and detailed below (dependencies: `scipy, pandas`).

The main figures are as follows:
* Figure 1: Summary of the 36 plasmid taxonomic units (PTUs) analysed in this paper.
* Figure 2: RM target depletion in the leading region of conjugative plasmids. This figure is an edited version (with schematic) of files generated in `results`: `6mer-leading-region-vs-rest-of-plasmid-density.pdf` and `Figure-2-6mer-deviation-overall.pdf`, both made by the `Paper-analysis.Rmd` Rmarkdown script. 
* Figure 3: GC-rich leading regions are more depleted in 6-bp palindromes.  
* Figure 4: Average depletion of within-range RM targets in the leading vs. lagging regions for conjugative PTUs. 

## Detailed workflow

We used a dataset of the 50 most prevalent plasmid taxonomic units (PTUs) in RefSeq200: 4,753 plasmids across 37 genera. These are available at `data/plasmid_table.csv`. The taxonomy for each genera is available at `data/genera_taxonomy.txt`
Within each PTU, we ran Prokka on each plasmid followed by Roary with a core threshold of 80% to generate pangenomes using scripts available [here](https://github.com/Adalijuanluo/Plasmid_pan). 
The pangenome data is available for download via Zenodo at [here](zenodo/pangenome_results). A summary of the pangenome stats is available at `data/PTU-pangenome-stats.csv` (Note that `PTU-B/O/K/Z` is renamed `PTU-BOKZ` and `PTU-L_M` becomes `PTU-LM` from this point forwards.) 
We then ran `scripts/make_core_accessory_fasta_from_roary.py` to generate core (in the paper we call this 'hard shell') and accessory fastas. 
These combined hard shell and accessory fastas for each PTU are available [here](zenodo/hard_shell_accessory_fastas). In these files the headers of the sequences give information on the genes with `seq_ID|gene_family` e.g. `>NC_019083.1_00042|group_17` 


We calculated and plotted overall pangenome statistics with `scripts/pangenome_statistics.R`

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
These can then be used to generate discrepancy scores for within-range targets using e.g.
```
python scripts/run-rmes-score-discrepancy.py --script scripts/rmes-score-discrepancy.py --dir ~/Downloads/trieste/ --kmers per-PTU --output results/rmes_discrepancies_targets_k6.csv  --k 6
# N.B. Wrapper script needs to be updated to pass correct arguments, given changes to rmes-score-discrepancy.py
```

Or one can run it like this:
```
while read f;                              
do
python scripts/rmes-score-discrepancy.py --file1 ~/Downloads/trieste/one_per_family/"$f"_core_one_per_family_k4.csv --file2 ~/Downloads/trieste/one_per_family/"$f"_accessory_one_per_family_subset_k4.csv --k 4 --kmers ~/Downloads/trieste/"$f"_targets.txt
done < test_PTUs.txt > test_one_per_family_within_range_targets_k4.csv
```


We identified leading regions of plasmids as outlined in the paper. Information on the leading region of PTUs is available at `data/leading_region_081123.tsv`. We then analyse this with `scripts/orient_leading_region.py`: 

```
python scripts/orient_leading_region.py --output test_output_window5000_step500.csv --fasta_file_search "/Users/Liam/Downloads/06_plaspan/*/_fasta/*" --leading_region_file data/leading_region_081123.tsv --leading_region_size 5000 --sliding_window 5000 --step_size 500 --k 6
```

We split the leading and lagging region (10kb) with 

```
 python scripts/split_leading_and_lagging.py --output_dir test_regions --fasta_file_search "/Users/Liam/Downloads/06_plaspan/*/_fasta/*" --leading_region_file data/leading_region_081123.tsv --leading_region_size 10000
```

We calculate the RMES scores for the leading vs. lagging region and then using those results calculate the difference with 

```
python scripts/rmes-score-discrepancy-leading-lagging.py --directory ~/Downloads/trieste/leading_lagging/ --kmers ~/Downloads/trieste/k6.txt --output test_6mer_leading_lagging --k 6
```

We then analyse the depletion in leading vs. lagging for within-range targets per PTU with `test_run_comparison_PTUs.sh` 
