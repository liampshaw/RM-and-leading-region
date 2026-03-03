# Assumes that plasmid_table_this_study is read in

scores_average_all = NULL
for (PTU in unique(plasmid_table_this_study$PTU)){
  for (k in c(5,6)){
    scores_df = read.csv(paste0('../results/targets_leading_lagging/', PTU, '_k', k, '_scores.csv'), header=T)
    n_targets = nrow(scores_df)
    scores_average <- scores_df %>%
      rowwise() %>%
      mutate(score_avg = mean(c_across(starts_with("score"))),
             sd=sd(c_across(starts_with("score")))) %>%
      ungroup() %>%
      select(kmer, score_avg, sd)
    counts_df = read.csv(paste0('../results/targets_leading_lagging/', PTU, '_k', k, '_counts.csv'), header=T)
    counts_average <- counts_df %>%
      rowwise() %>%
      mutate(count_avg = mean(c_across(starts_with("count")))) %>%
      ungroup() %>%
      select(kmer, count_avg)
    scores_average = scores_average %>% left_join(counts_average, by="kmer")
    scores_average$PTU = PTU
    scores_average$k = nchar(scores_average$kmer)
    scores_average_all = rbind(scores_average_all, scores_average)
      }
}

scores_average_all = data.frame(scores_average_all)
scores_average_all = scores_average_all %>% left_join(mob.types, by="PTU")
scores_average_all$Majority.MOB = ordered(scores_average_all$Majority.MOB,
                                          levels=c("MOBF", "MOBP", "MOBH"))
p.k.5 = ggplot(scores_average_all %>% filter(k==5), aes(Majority.MOB, group=PTU, score_avg, fill=Majority.MOB))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_boxplot(position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=friendly_pal("contrast_three"))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid = element_blank())+
  ylab("Relative representation score\n(leading 10kb vs. lagging 10kb)")+
  xlab("")
p.k.6 = ggplot(scores_average_all %>% filter(k==6), aes(Majority.MOB, group=PTU, score_avg, fill=Majority.MOB))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_boxplot(position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=friendly_pal("contrast_three"))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid = element_blank())+
  ylab("Relative representation score\n(leading 10kb vs. lagging 10kb)")+
  xlab("")

p.targets.combined = plot_grid(p.k.5+ggtitle("(a) 5-bp within-range targets"), 
          p.k.6+ggtitle("(b) 6-bp within-range targets"), nrow=1)
ggsave(file="../results/fig-within-range-targets.pdf", width=10, height=4, p.targets.combined)



colnames(results) = c("PTU", "k", "n.targets", "median.score", "median.count")
results$median.score= as.numeric(results$median.score)
results$median.count = as.numeric(results$median.count)

# Order by median score for k=6
order.k.6.scores = results[which(results$k==6), "PTU"][order(results[which(results$k==6), "median.score"])]
results$PTU = ordered(results$PTU,
                      levels=rev(order.k.6.scores))

p.figure4 = ggplot(results, aes(median.score, as.factor(PTU)))+
  geom_vline(xintercept = 0, size=1)+
  geom_bar(stat="identity", position="dodge")+
  facet_wrap(~k, scales="free_x")+
  ylab("")+
  xlab("Median depletion score of within-range RM targets\n(leading vs. lagging, 10kb)")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour="black", size=12),
        axis.title.x=element_text(size=24),
        axis.text.x=element_text(size=18))+
  theme(axis.text.y=element_text(hjust=0))
ggsave(p.figure4, file='../results/Figure-4-figure-depletion-within-range-targets.pdf',
       width=10, height=6)