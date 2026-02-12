# Plot the codon usage data for CDSs

codon_usage = read.csv('../../2026-new-data/codon_usage.tsv', 
         header=T, 
         stringsAsFactors = F,
         sep='\t',
         row.names = 1)

# Prepare data for ordination
codon.mat = as.matrix(codon_usage)
codon.mat.scaled = scale(codon.mat)

# Run PCA
pca <- prcomp(codon.mat.scaled, rank. = 10)


# Run umap
umap_res <- umap(
  codon.mat.scaled,
  n_components = 2,
  n_neighbors = 15,
  min_dist = 0.1,
  verbose = TRUE
)

