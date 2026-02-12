# Plot the codon usage data for CDSs
library(ggplot2)
library(compositions)
library(uwot) # for umap, if needed

plasmid_table = read.csv('/Users/Liam/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/data/plasmid_table.tsv', header=T, sep='\t', 
                         row.names = 1)
# Read in codon usage per CDS in all plasmids
codon_usage = read.csv('../../2026-new-data/codon_usage.tsv', 
          header=T, 
          stringsAsFactors = F,
          sep='\t',
          row.names = 1)
# rscu_res = read.csv('/Users/Liam/Downloads/new_1751/all_rscu.tsv', 
    #      header=T, 
    #      stringsAsFactors = F,
    #      sep='\t',
    #      row.names = 1)
row.names(codon_usage) = gsub(".*\\|", "", rownames(codon_usage))

# Prepare data for ordination
codon.mat = as.matrix(codon_usage)
codon.mat.norm = codon.mat / rowSums(codon.mat)
codon.mat.scaled = clr(codon.mat.norm)

# Run PCA
pca <- prcomp(codon.mat.scaled, rank. = 10)

# Create dataframe
pca.df = data.frame(pca$x)
pca.df$position = as.numeric(gsub(".*_", "", rownames(pca.df)))
pca.df$plasmid = gsub("\\.2_.*", "", rownames(pca.df))
pca.df$plasmid = gsub("\\.1_.*", "", pca.df$plasmid)
pca.df$PTU = plasmid_table[pca.df$plasmid, "PTU"]

ds_hits = read.csv('~/Downloads/new_1751/ds_results.tsv', header=T, sep='\t')
antidefense_hits = ds_hits[which(ds_hits$activity=="Antidefense"),]
defense_hits = ds_hits[which(ds_hits$activity=="Defense"),]

pca.df$antidefense = ifelse(rownames(pca.df) %in% antidefense_hits$sys_beg, "antidefense", "")
pca.df$defense = ifelse(rownames(pca.df) %in% defense_hits$sys_beg, "defense", "")
pca.df$type = ifelse(pca.df$antidefense=="antidefense", "antidefense", ifelse(pca.df$defense=="defense", "defense", "other"))


ggplot(pca.df, aes(PC1, PC2, colour=ds))+
  geom_point()+
  facet_wrap(~type)

ggplot(pca.df[which(pca.df$ds=="antidefense"),], aes(PC2, PC3, colour=ds))+
  geom_point()
ggplot(pca.df[which(pca.df$ds=="other"),], aes(PC2, PC3, colour=ds))+
  geom_point()

ggplot(pca.df[which(pca.df$plasmid %in% head(unique(pca.df$plasmid) ) & pca.df$position>=10),], aes(PC1, PC2))+
  geom_point()+
  theme_bw()+
  facet_wrap(~plasmid)+
  geom_point(data=pca.df[which(pca.df$plasmid %in% head(unique(pca.df$plasmid) ) & pca.df$position<10),], aes(fill=position), shape=21)+
  scale_fill_continuous(low="white", high="red")

# Or try umap?
# Run umap - takes several minutes
umap_res <- umap(
  codon.mat,
  n_components = 2,
  n_neighbors = 15,
  min_dist = 0.1,
  verbose = TRUE
)
umap_res_df = data.frame(umap_res)
umap_res_df$position = as.numeric(gsub(".*_", "", rownames(umap_res_df)))
umap_res_df$plasmid = gsub("\\.2_.*", "", rownames(umap_res_df))
umap_res_df$plasmid = gsub("\\.1_.*", "", umap_res_df$plasmid)
umap_res_df$PTU = plasmid_table[umap_res_df$plasmid, "PTU"]
ggplot(umap_res_df[which(umap_res_df$PTU=="PTU-BOKZ"  & umap_res_df$position>=10),], aes(X1, X2))+
  geom_point()+
  theme_bw()+
  geom_point(data=umap_res_df[which(umap_res_df$PTU=="PTU-BOKZ" & umap_res_df$position<10),], aes(fill=position), shape=21)+
  scale_fill_continuous(low="white", high="red")

ggplot(umap_res_df[which(umap_res_df$plasmid=="NZ_CP019253"),], aes(X1, X2, colour=position))+
  geom_point()+
  theme_bw()+
  geom_path(aes(group=1))
