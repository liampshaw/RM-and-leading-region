# Restriction-modification targets and the leading region of plasmids

**Note: repository currently undergoing streamlining**

A repository of analysis related to the paper:

**The leading region of many conjugative plasmids is depleted in restriction-modification targets**  
Liam P. Shaw, Arancha Peñil Celis, Manuel Ares Arroyo, Lijuan Luo, Tatiana Dimitriu, Fernando de la Cruz 
biorxiv: https://www.biorxiv.org/content/10.1101/2025.04.03.647016v1  

Associated datasets can be downloaded at Zenodo: [here](zenodo)

## Paper analysis

The main analysis presented in the paper is performed in `scripts/Paper-analysis.Rmd`. Associated python scripts used for data processing before then are in `scripts` and detailed below (dependencies: `scipy, pandas`).

The main figures are as follows:
* Figure 1: Summary of the 36 plasmid taxonomic units (PTUs) analysed in this paper.
* Figure 2: RM target depletion in the leading region of conjugative plasmids. This figure is an edited version (with schematic) of files generated in `results`: `6mer-leading-region-vs-rest-of-plasmid-density.pdf` and `Figure-2-6mer-deviation-overall.pdf`, both made by the `Paper-analysis.Rmd` Rmarkdown script. 
* Figure 3: GC-rich leading regions are more depleted in 6-bp palindromes.  
* Figure 4: Average depletion of within-range RM targets in the leading vs. lagging regions for conjugative PTUs. 

