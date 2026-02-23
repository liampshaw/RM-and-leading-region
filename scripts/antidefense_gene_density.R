# Plot antidefense gene density
defensefinder = read.csv('~/Dropbox/_Projects/2023-Trieste/2026-new-data/defensefinder_results_prodigal_meta.tsv',
                         header=T,
                         sep='\t')
plasmid_table = read.csv('/Users/Liam/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/data/plasmid_table.tsv', header=T, sep='\t', 
                         row.names = 1)

antidefense = defensefinder[which(defensefinder$activity=="Antidefense"),]
  antidefense$position = as.numeric(gsub(".*_", "", antidefense$sys_beg))
  # Add PTU
  antidefense$plasmid = sub("_[0-9]+$", "", antidefense$sys_beg)
  antidefense$PTU = plasmid_table[antidefense$plasmid, "PTU"]
  p.all = ggplot(antidefense, aes(position, fill=type))+
    geom_histogram()+theme_bw()+
    xlab("Start of system (ORF position)")+
    theme(legend.position = c(0.8, 0.8))
  p.facet = ggplot(antidefense, aes(position, fill=type))+
  geom_histogram()+theme_bw()+
  xlab("Start of system (ORF position)")+
  facet_wrap(~PTU, scales="free_y")+
    theme(legend.position = "none")

  pdf('../results/figure-antidefense-gene-density.pdf', width=8, height=4)
cowplot::plot_grid(p.all, p.facet, nrow=1)
dev.off()