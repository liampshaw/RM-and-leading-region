# Restriction-modification targets and the leading region of plasmids

A repository of analysis related to the paper:

**The leading region of many conjugative plasmids is depleted in restriction-modification targets**  
Liam P. Shaw, Arancha Peñil Celis, Manuel Ares Arroyo, Lijuan Luo, Tatiana Dimitriu, Fernando de la Cruz  
doi: [10.1101/2025.04.03.647016](https://doi.org/10.1101/2025.04.03.647016) 

Associated datasets can be downloaded at Zenodo: [10.5281/zenodo.15235227](https://doi.org/10.5281/zenodo.15235227)

The main analysis presented in the paper is performed in `scripts/Paper-analysis.Rmd`. Associated python scripts used for data processing before then are in `scripts` and detailed in `scripts/detailed-workflow.md` (dependencies: `scipy, pandas`).


## Data 

`1751_plasmids_oriented.tar.gz` - n=1751 plasmids oriented so that leading strand comes first
`1751_defense_systems.tsv` - results of `defense-finder -a` (with antidefense systems detected) on n=1751 plasmids
*to be added* - CDS of all n=1751 plasmids

## Palindrome depletion

Fasta files for plasmids have been correctly oriented so that the start of the fasta file corresponds to the leading strand. 

We then count palindromes using `scripts/rust/src/palindrome_density.rs`. This script can work for any value of k - we focus on k=4,6,8,10 because Type II RM targets are typically 4-8bp (and most are 5-6bp). We can use a sliding window or not, but we want to avoid autocorrelation for some statistics so in that case we don't use a sliding window. 

