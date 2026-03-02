# Plot antidefense gene density
library(ggplot2)
library(dplyr)
defensefinder = read.csv('~/Dropbox/_Projects/2023-Trieste/2026-new-data/defensefinder_results_prodigal_meta.tsv',
                         header=T,
                         sep='\t')
defensefinder$plasmid = sub("_[0-9]+$", "", defensefinder$sys_beg)
defensefinder$position = as.numeric(gsub(".*_", "", defensefinder$sys_beg))
# Add PTU
defensefinder$plasmid = sub("_[0-9]+$", "", defensefinder$sys_beg)
defensefinder$PTU = plasmid_table[defensefinder$plasmid, "PTU"]

# Plasmid table
plasmid_table = read.csv('/Users/Liam/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/data/plasmid_table.tsv', header=T, sep='\t', 
                         row.names = 1)
# Subset to this study
plasmid_ids_this_study = read.csv('~/Dropbox/_Projects/2023-Trieste/2026-new-data/1751_plasmid_accs.txt', 
                                  header=F)$V1
plasmid_table = plasmid_table[plasmid_ids_this_study,]

# Values for plotting
colour.values = read.csv('../../2026-new-data/colour-values-antidefense.csv',
                         header=T)
colour.values = colour.values %>% arrange(Subtype)
antidefense = defensefinder[which(defensefinder$activity=="Antidefense"),]
# Change "Other" to "PsiAB" - only works for this dataset!
antidefense$type[antidefense$type=="Other"] = "psiAB"

# Exclude one Anti_Pycsar for purpose of plotting
antidefense.plot = antidefense[which(antidefense$type!="Anti_Pycsar"),]
  p.all = ggplot(antidefense.plot, aes(position, fill=type))+
    geom_histogram()+theme_bw()+
    xlab("Start of system (ORF position)")+
    ylab("Frequency")+
    theme(legend.position = c(0.8, 0.8),
          legend.background = element_rect(colour="black"))+
    scale_fill_manual(values=c("#fecc5c","#fd8d3c","#e31a1c"))

# Normalise facet by total number of plasmids within each PTU?
  # Add majority MOB type for PTUs
  mob.types = read.csv('../../2026-new-data/PTU-mob-types.csv',
                       header=T)
  antidefense.plot = antidefense.plot %>% left_join(mob.types, by="PTU")

  facetPlot = function(df){
    return(ggplot(df, aes(position, fill=type))+
      geom_histogram()+theme_bw()+
      xlab("Start of system (ORF position)")+
      facet_wrap(~PTU, scales="free_y")+
      theme(legend.position = "none")+
      theme(panel.grid=element_blank())+
      scale_fill_manual(values=c("#fecc5c","#fd8d3c","#e31a1c"))+
      ylab("Frequency"))
  }
  
p.facet.MOBF = facetPlot(antidefense.plot %>% filter(Majority.MOB=="MOBF"))
p.facet.MOBH = facetPlot(antidefense.plot %>% filter(Majority.MOB=="MOBH"))
p.facet.MOBFandH = facetPlot(antidefense.plot %>% filter(Majority.MOB=="MOBH" | Majority.MOB=="MOBF"))

p.facet.MOBP = facetPlot(antidefense.plot %>% filter(Majority.MOB=="MOBP"))
p.facet.all = facetPlot(antidefense.plot)
ggsave(p.all, file="../results/figure-facet-antidefense-all.pdf", width=6, height=4.5)

ggsave(p.facet.MOBF, file="../results/figure-facet-antidefense-MOBF.pdf", width=6, height=4)
ggsave(p.facet.MOBH, file="../results/figure-facet-antidefense-MOBH.pdf", width=6, height=4)
ggsave(p.facet.MOBFandH, file="../results/figure-facet-antidefense-MOBFandH.pdf", width=6, height=4) # to get ratio right

ggsave(p.facet.MOBP, file="../results/figure-facet-antidefense-MOBP.pdf", width=6, height=4)
# These are then combined and edited in fig-facet-antidefense-PTU.svg

# For each PTU, calculate the proportion of its plasmids 
# which carry an anti-defense system
plasmid_table$antidefense.detected = ifelse(rownames(plasmid_table) %in% unique(antidefense$plasmid),
                                                                                "yes", "no")
antidefense.PTU.count = plasmid_table %>% group_by(PTU) %>%
  summarise(n.antidefense=length(antidefense.detected[antidefense.detected=="yes"]),
            n.total=length(antidefense.detected)) %>%
  mutate(prop=n.antidefense/n.total * 100)


pdf('../results/figure-antidefense-gene-density.pdf', width=8, height=4)
cowplot::plot_grid(p.all, p.facet, nrow=1)
dev.off()


# Give genus counts for each plasmid
plasmid_table$Genus.simple = sapply(plasmid_table$Genus, 
                                    function(x) ifelse(x %in% head(names(sort(table(plasmid_table$Genus), 
                                                                             decreasing = TRUE)), n=6), x, "Other"))
table(plasmid_table$Genus.simple)
PTU.by.genus = plasmid_table %>% group_by(PTU, Genus.simple) %>%
  summarise(n=length(Genus.simple)) %>% mutate(prop=n/sum(n))
PTU.by.genus$Genus.simple = ordered(PTU.by.genus$Genus.simple,
                                    levels=c("Escherichia", "Shigella", "Salmonella",
                                             "Citrobacter", "Enterobacter", "Klebsiella", "Other"))
p.genus.props = ggplot(PTU.by.genus, aes(PTU, prop, fill=Genus.simple))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values=c("#2B8CBE", "#A6BDDB", "#A6D854",  "#66C2A5","#FC8D62", "#E78AC3", "grey"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()
ggsave(p.genus.props, file='../results/fig-prop-genus-counts.pdf', width=4, height=4)
