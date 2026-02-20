library(ggplot2)
library(dplyr)
library(tidyr)

# Read in CAI data
cai = read.csv('~/Dropbox/_Projects/2023-Trieste/2026-new-data/cai_prodigal_meta.tsv',
               header=T,
               sep='\t')
defensefinder = read.csv('~/Dropbox/_Projects/2023-Trieste/2026-new-data/defensefinder_results_prodigal_meta.tsv',
                         header=T,
                         sep='\t')
plasmid_table = read.csv('/Users/Liam/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/data/plasmid_table.tsv', header=T, sep='\t', 
                         row.names = 1)
cai$id = gsub(".*\\|", "", cai$gene_id)
cai$plasmid = gsub("\\.fna.*", "", gsub(".*\\/", "", cai$gene_id))
cai$position = as.numeric(gsub(".*_", "", cai$gene_id))
cai$PTU = plasmid_table[cai$plasmid, "PTU"]
cai$HostGenus = plasmid_table[cai$plasmid, "Genus"]

# check just included plasmids
plasmid_table_this_study = plasmid_table[unique(cai$plasmid),]


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
  geom_point(size=0.5, aes(colour=mean.cai))+
  facet_wrap(~PTU, scales="free_x")+
  scale_colour_continuous(low="white", high="red")+
  xlab("ORF position")+
  ylab("Mean CAI")+
  theme_bw()
dev.off()
# And average over all PTUs at each positino
cai.ptu %>% group_by(position) %>%
  summarise(mean.cai=mean(mean.cai)) %>%
  ggplot(aes(position, mean.cai))+
  geom_point(size=0.5, aes(colour=mean.cai))+
  scale_colour_continuous(low="white", high="red")+
  xlab("ORF position")+
  ylab("Mean CAI")+
  theme_bw()+
  xlim(c(0,100))

# Flip them as well
cai.flipped <- cai %>%
  group_by(plasmid) %>%
  mutate(position = -(max(position) + 1 - position))%>% # 1 -> -N, N -> -1
  ungroup()
cai.both = rbind(cai, cai.flipped)
cai.both.max.positions = cai.both %>%
  group_by(PTU, plasmid) %>%
  summarise(max=max(position)) %>%
  group_by(PTU) %>%
  summarise(percentile=quantile(max, 0.25)) # Choose a percentile threshold e.g. 0.2 means at this length 75% of plasmids in the PTU are at least this length
cai.ptu.both = cai.both %>%
  group_by(PTU, position) %>%
  summarise(mean.cai=mean(cai)) 
cai.ptu.both = cai.ptu.both %>%
  left_join(cai.both.max.positions, by = "PTU")
# Filter out 
cai.ptu.both.filter = cai.ptu.ecoli.both %>%
  filter(position <= percentile/2 & position >= -percentile/2)

ggplot(cai.ptu.both.filter, aes(position, mean.cai))+
  geom_point(size=0.5)+
  facet_wrap(~PTU, scales="free_x",nrow=3)

# Make plots for 
MOBF.PLASMIDS = c("PTU-N1", "PTU-FE", "PTU-FS")
MOBP.PLASMIDS = c("PTU-I1", "PTU-I2", "PTU-BOKZ", "PTU-LM", 
                  "PTU-X1", "PTU-X3", "PTU-E81", "PTU-E84")
MOBH.PLASMIDS = c("PTU-HI1A", "PTU-C")

p.mobF = ggplot(cai.ptu.both.filter %>% filter(PTU %in% MOBF.PLASMIDS), aes(position, mean.cai))+
  geom_point(size=0.5, aes(colour=mean.cai))+
  facet_wrap(~PTU, scales="free_x", nrow=1)+
  scale_colour_continuous(low="white", high="red")
p.mobP = ggplot(cai.ptu.both.filter %>% filter(PTU %in% MOBP.PLASMIDS), aes(position, mean.cai))+
  geom_point(size=0.5, aes(colour=mean.cai))+
  facet_wrap(~PTU, nrow=1,scales="free_x")+
  scale_colour_continuous(low="white", high="red")
p.mobH = ggplot(cai.ptu.both.filter %>% filter(PTU %in% MOBH.PLASMIDS), aes(position, mean.cai))+
  geom_point(size=0.5, aes(colour=mean.cai))+
  facet_wrap(~PTU,nrow=1, scales="free_x")+
  scale_colour_continuous(low="white", high="red")


##############################
# Just Escherichia coli plasmids
##############################
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


# Analysing defensefinder results
antidefense = defensefinder[which(defensefinder$activity=="Antidefense"),]
defense = defensefinder[which(defensefinder$activity=="Defense"),]

antidefense_proteins <- antidefense %>%
  separate_rows(protein_in_syst, sep = ",") %>%
  pull(protein_in_syst) %>%
  unique()
defense_proteins <- defense %>%
  separate_rows(protein_in_syst, sep = ",") %>%
  pull(protein_in_syst) %>%
  unique()
cai <- cai %>%
  mutate(protein_type = if_else(id %in% antidefense_proteins, "antidefense", 
                               ifelse(id %in% defense_proteins, "defense", "other")),
         leading_region = if_else(position < 20, " In first 25 ORFs", "Beyond first 25 ORFs"))
# Plot of antidefense/defense/other genes
pdf('~/Dropbox/_Projects/2023-Trieste/RM-and-leading-region/results/Figure-CAI-gene-categories.pdf')
ggplot(cai, aes(leading_region, cai, group=interaction(leading_region, protein_type),
                colour=protein_type))+
  geom_boxplot()+
  xlab("")+
  ylab("Codon Adaptation Index (CAI)")+
  labs(colour="Type")+
  theme_bw()+
  scale_color_brewer(palette="Set2")
dev.off()

# Plot of antidefense genes by position
ggplot(cai %>% filter(protein_type=="antidefense") %>% 
         group_by(PTU, position) %>%
         summarise(median_cai=median(cai)), aes(position, median_cai))+
  geom_point()+
  facet_wrap(~PTU)

# Test 99th percentile for CAI with Chi-squared test
cai$antidefense = ifelse(cai$protein_type=="antidefense", "antidefense", "other")
chisq.test(table(cai$antidefense, cai$cai>quantile(cai$cai, 0.99)))

# Test whether this holds per PTU?
ggplot(cai %>% group_by(PTU) %>%
  summarise(median_cai_diff=median(cai[protein_type=="antidefense"])-median(cai[protein_type=="other"])),
  aes(PTU, median_cai_diff)) +
  geom_bar(stat="identity")


# ANtidefense *type*
antidefense_by_protein = antidefense %>% 
  separate_rows(protein_in_syst, sep=",") %>%
  mutate(id=protein_in_syst)
cai_antidefense = cai %>% filter(protein_type=="antidefense") %>%
 left_join(antidefense_by_protein, by="id") %>%
  group_by(type, subtype) %>%
  summarise(median_cai=median(cai), 
            n=length(cai), 
            n.PTU=length(unique(PTU)),
            min_cai=min(cai),
            max_cai=max(cai)) %>%
  arrange(median_cai)
  
cai_all_antidefense = cai %>% filter(protein_type=="antidefense") %>%
  left_join(antidefense_by_protein, by="id")
