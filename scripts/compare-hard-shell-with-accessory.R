# Compare hard shell genes vs. accessory
library(ggplot2)

# Location of hard shell (core) and accessory fasta files
FASTA_DIR = '~/Downloads/trieste/'

# Read in PTU stats and only keep those with >9 hard shell genes
ptu.stats = read.csv('../data/PTU-pangenome-stats.csv', 
                     header=T)
SUBSET_PTUS = ptu.stats$PTU[ptu.stats$pangenome.category=="include"]
  

compareHardShellVsAccessory = function(PTU_name, k){
  ptu_df_hard_shell = read.csv(paste0(FASTA_DIR, 
                                      PTU_name,
                                      '_core_genes_k', k, '.csv'), 
                               header=T)
    ptu_df_accessory = read.csv(paste0(FASTA_DIR, 
                                     PTU_name,
                                     '_accessory_genes_k', k, '.csv'), 
                              header=T)
  # Normalise scores by maximum (for easier comparison)
  ptu_df_hard_shell$score.norm = ptu_df_hard_shell$score/max(ptu_df_hard_shell$score)
  ptu_df_accessory$score.norm = ptu_df_accessory$score/max(ptu_df_accessory$score)
  # There would be a better scaling here...
  
  ptu_df = merge(ptu_df_hard_shell, ptu_df_accessory, by="word")
  
  palindromes = tolower(read.csv(paste0('~/Dropbox/_Projects/2023-Trieste/2023-Trieste-RM/palindromes-k', k, '.txt'),
                                 header=F)$V1)
  ptu_df$palindrome = ifelse(ptu_df$word %in% palindromes, "palindrome", "non.palindrome")
  ptu_df$score.discrepancy = ptu_df$score.x-ptu_df$score.y
  
  # ggplot(ptu_df, aes(score.norm.x, score.norm.y))+
  #   geom_point()+
  #   theme_minimal()+
  #   facet_wrap(~palindrome)
  p = ggplot(ptu_df, aes(palindrome, score.discrepancy))+
    ggbeeswarm::geom_quasirandom()+
    #geom_boxplot()+
    ggsignif::geom_signif(comparisons = list(c("palindrome", "non.palindrome")),
                          test = "wilcox.test")+
    theme_bw()+
    ggtitle(PTU_name)+
    stat_summary(fun="median",geom="errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, colour="red")+
    ylab("Discrepancy in avoidance scores")+
    xlab("")
  return(p)
  
}


for (PTU_name in SUBSET_PTUS[1:length(SUBSET_PTUS)]){
  k = 6
  p = compareHardShellVsAccessory(PTU_name, k=k)
  ggsave(p, file=paste0('../figures/PTU-plots/2024-10-25-', PTU_name, '-palindromes-k-', k, '-hard-shell.pdf'),
         width=6, height=4)
  k = 4
  p = compareHardShellVsAccessory(PTU_name, k=k)
  ggsave(p, file=paste0('../figures/PTU-plots/2024-10-25-', PTU_name, '-palindromes-k-', k, '-hard-shell.pdf'),
         width=6, height=4)
}

hard_shell_list <- list.files(path = "~/Dropbox/_Projects/2023-Trieste/2023-Trieste-RM/output/", pattern = "core_genes_k6\\.csv$", full.names = TRUE)
accessory_list <- list.files(path = "~/Dropbox/_Projects/2023-Trieste/2023-Trieste-RM/output/", pattern = "accessory_genes_k6\\.csv$", full.names = TRUE)

ptu_info = read.csv('../data/top50ptus_info.csv', header=T, row.names = 1)

median_discrepancies = c()
wilcox_test_results = c()
for (i in 1:length(hard_shell_list)){
  ptu_df_hard_shell = read.csv(hard_shell_list[i],
                               header=T)
  ptu_df_accessory = read.csv(accessory_list[i],
                              header=T)
  ptu_df = merge(ptu_df_hard_shell, ptu_df_accessory, by="word")
  palindromes = tolower(read.csv(paste0('~/Dropbox/_Projects/2023-Trieste/2023-Trieste-RM/palindromes-k6.txt'),
                                 header=F)$V1)
  ptu_df$palindrome = ifelse(ptu_df$word %in% palindromes, "palindrome", "non.palindrome")
  ptu_df$score.discrepancy = ptu_df$score.x-ptu_df$score.y
  median_discrepancy = median(ptu_df$score.discrepancy[which(ptu_df$palindrome=="palindrome")])
  median_discrepancies = c(median_discrepancies, median_discrepancy)
  wilcox_test_result = wilcox.test(x=ptu_df$score.discrepancy[which(ptu_df$palindrome=="palindrome")],
                                   y=ptu_df$score.discrepancy[which(ptu_df$palindrome=="non.palindrome")])$p.value
  wilcox_test_results = c(wilcox_test_results, wilcox_test_result)
}

results_df = data.frame(discrepancy=median_discrepancies, 
                        wilcox.test.result = wilcox_test_results,
                        PTU=gsub("_core.*", "", gsub(".*output\\/\\/", "", hard_shell_list)))
results_df$host.range = ptu_info[results_df$PTU,"host.range"]

results_df$PTU = ordered(results_df$PTU,
                         levels=results_df$PTU[order(results_df$discrepancy)])

ggplot(results_df, aes(PTU, discrepancy, colour=wilcox.test.result<0.05))+
  geom_point()+
  stat_summary(fun="median", aes(group=1), geom="errorbar")+
  facet_wrap(~host.range)
