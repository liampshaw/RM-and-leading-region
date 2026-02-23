# Palindrome density
d = read.csv('../../2026-new-data/palindrome_density.tsv', header=T, sep='\t')
d$PTU = plasmid_table_this_study[d$contig,"PTU"]
library(ggplot2)

d$start.factor.kb = as.factor(d$start/1000)
pdf('../results/Figure-palindrome-density-einverted.pdf', width=8, height=4)
ggplot(d[which(d$start<50000),], 
       aes(start.factor.kb, count))+
  geom_boxplot()+facet_wrap(~PTU, scales="free_y")+
  xlab("Position (kb)")+
  ylab("Number of palindromes per 5kb")+
  theme_bw()
dev.off()


# Potential modelling
library(lme4)
d$leading.region = d$start < 20000
model <- glmer(
  count ~ leading.region + PTU + (1 | contig),
  data = d,
  family = poisson(link = "log")
)