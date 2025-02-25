# Get average size of PTU
ptu = read.csv('../data/PTU-pangenome-stats.csv')
for (each_ptu in ptu$PTU){
  gene_pa_file <- list.files(paste0("/Users/Liam/Downloads/06_plaspan/", each_ptu), 
                        pattern = "gene_presence_absence.Rtab$", full.names = TRUE, recursive = TRUE)
  gene_counts = read.csv(gene_pa_file, sep='\t', row.names = 1)
  ptu$median.genes[which(ptu$PTU==each_ptu)] = median(colSums(gene_counts))
  ptu$max.genes[which(ptu$PTU==each_ptu)] = max(colSums(gene_counts))
  ptu$min.genes[which(ptu$PTU==each_ptu)] = min(colSums(gene_counts))
}

# Plot these
ptu_for_plotting = ptu %>% pivot_longer(cols=c("n.core", "n.others"))

ptu = ptu[order(ptu$median.genes),]
ptu$PTU = ordered(ptu$PTU, 
                  levels=rev(ptu$PTU))
ptu$leading.region = ifelse(ptu$PTU %in% leading.region.PTUs,
                            "yes", "no")
pdf('../results/ptu-summary-pangenome-plot.pdf', width=9, height=5)
ggplot(ptu[which(ptu$median.genes>20),], aes(PTU, median.genes, ymin=min.genes, ymax=max.genes,
                                             fill=leading.region))+
  geom_bar(aes(y=n.core), stat="identity")+
  geom_errorbar()+
  geom_point(shape=21)+
  theme_bw()+
  xlab("")+
  ylab("Number of genes")+
  labs(fill="Leading region analysed")+
  coord_flip()+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  theme(axis.text=element_text(colour="black"))
dev.off()

