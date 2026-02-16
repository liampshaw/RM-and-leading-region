library(ggplot2)
library(dplyr)

# Read in CAI data
cai = read.csv('~/Dropbox/_Projects/2023-Trieste/2026-new-data/cai.tsv',
               header=T,
               sep='\t')
plasmid_table = read.csv('/Users/Liam/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/data/plasmid_table.tsv', header=T, sep='\t', 
                         row.names = 1)

cai$plasmid = gsub("\\.fna.*", "", gsub(".*\\/", "", cai$gene_id))
cai$position = as.numeric(gsub(".*_", "", cai$gene_id))
cai$PTU = plasmid_table[cai$plasmid, "PTU"]
cai$HostSpecies = plasmid_table[cai$plasmid, "Species"]

cai.ptu = cai %>% group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai))
# Add a length threshold for each PTU
cai.ptu.thresholds = cai %>%
  group_by(PTU, plasmid) %>%
  summarise(max=max(position)) %>%
  group_by(PTU) %>%
  summarise(percentile=quantile(max, 0.25)) 
cai.ptu = cai.ptu %>% left_join(cai.ptu.thresholds, by="PTU")

pdf('~/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/results/cai_results.pdf')
ggplot(cai.ptu %>% filter(position <= percentile), aes(position, mean.cai))+
  geom_point(size=0.5)+
  facet_wrap(~PTU, scales="free_x")
dev.off()

# Just Escherichia coli plasmids
cai.ecoli = cai %>% filter(HostSpecies=="Escherichia coli")
cai.ptu.ecoli = cai.ecoli %>%
  group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai))
ggplot(cai.ptu.ecoli, aes(position, mean.cai))+
  geom_point(size=0.5)+
  facet_wrap(~PTU, scales="free_x")


# Try inverting so that we have the lagging region going from -1 to -N
cai.ecoli.flipped <- cai.ecoli %>%
  group_by(plasmid) %>%
  mutate(position = -(max(position) + 1 - position))%>%
  ungroup()
cai.ptu.ecoli.flipped = cai.ecoli.flipped %>%
  group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai))
ggplot(cai.ptu.ecoli.flipped, aes(position, mean.cai))+
  geom_point(size=0.5)+
  facet_wrap(~PTU, scales="free_x")
cai.ecoli.both = rbind(cai.ecoli.flipped, cai.ecoli)
ggplot(cai.ecoli.both[which(cai.ecoli.both$PTU=="PTU-C"),], aes(position, cai))+
  geom_point(size=0.5)+
  facet_wrap(~plasmid)
cai.ptu.ecoli.both = cai.ecoli.both %>%
  group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai)) 

ggplot(cai.ptu.ecoli.both, aes(position, mean.cai))+
  geom_point(size=0.5)+
  facet_wrap(~PTU, scales="free_x")+
  xlim(c(-100, 100))
# Check maximum genes for each PTU and divide by two
cai.ecoli.both.max.positions = cai.ecoli.both %>%
  group_by(PTU, plasmid) %>%
  summarise(max=max(position)) %>%
  group_by(PTU) %>%
  summarise(percentile=quantile(max, 0.25)) # Choose a percentile threshold e.g. 0.2 means at this length 75% of plasmids in the PTU are at least this length
cai.ptu.ecoli.both = cai.ptu.ecoli.both %>%
  left_join(cai.ecoli.both.max.positions, by = "PTU")
# Filter out 
cai.ptu.ecoli.both.constrained = cai.ptu.ecoli.both %>%
  filter(position <= percentile/2 & position >= -percentile/2)

ggplot(cai.ptu.ecoli.both.constrained, aes(position, mean.cai))+
  geom_point(size=0.5)+
  facet_wrap(~PTU, scales="free_x")
  