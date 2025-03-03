# Analyse palindrome densities in hard shell vs. accessory
library(ggplot2)
palindromes_k4 = tolower(read.csv('~/Downloads/trieste/k4.txt', header=F)$V1)
palindromes_k6 = tolower(read.csv('~/Downloads/trieste/k6.txt', header=F)$V1)

palindromeDensity = function(rmes_table, k=6){
  rmes_df = read.csv(rmes_table)
  if (k==4){
    palindrome_count = sum(rmes_df[which(rmes_df$word %in% palindromes_k4),"count"])
  }
  if (k==6){
    palindrome_count = sum(rmes_df[which(rmes_df$word %in% palindromes_k6),"count"])
  }
  return(palindrome_count/sum(rmes_df$count) * 1000)
}

# Go through PTUs
PTUs_of_interest = read.csv('../test_PTUs.txt', header=F)$V1
results_df_k4 = data.frame()
results_df_k6 = data.frame()
for (PTU in PTUs_of_interest){
  results_df_k6 = rbind(results_df_k6, c(PTU, palindromeDensity(paste0('~/Downloads/trieste/', PTU, '_core_genes_k6.csv')),
    palindromeDensity(paste0('~/Downloads/trieste/', PTU, '_accessory_genes_k6.csv'))))
  results_df_k4 = rbind(results_df_k4, c(PTU, palindromeDensity(paste0('~/Downloads/trieste/', PTU, '_core_genes_k4.csv'), 4),
                                      palindromeDensity(paste0('~/Downloads/trieste/', PTU, '_accessory_genes_k4.csv'), 4)))
  
}
colnames(results_df_k6) = c("PTU", "core", "accessory")
results_df_k6$core = as.numeric(results_df_k6$core)
results_df_k6$accessory = as.numeric(results_df_k6$accessory)
results_df_k6$diff = results_df_k6$core-results_df_k6$accessory
ggplot(results_df_k6, aes(diff))+
  geom_boxplot()

colnames(results_df_k4) = c("PTU", "core", "accessory")
results_df_k4$core = as.numeric(results_df_k4$core)
results_df_k4$accessory = as.numeric(results_df_k4$accessory)
results_df_k4$diff = results_df_k4$core-results_df_k4$accessory

median(results_df_k4$diff)
median(results_df_k6$diff)

wilcox.test((results_df_k6$diff))
wilcox.test((results_df_k4$diff))


