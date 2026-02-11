library(data.table)
library(zoo)
library(ggplot2)
library(dplyr)

plasmid_table = read.csv('../../data/plasmid_table.tsv', header=T, sep='\t', 
                         row.names = 1)

# Read in results for palindromes of length k=4,6,8
# that have been produced with rust pal_density script
average_density = function(results_filename){
  dt <- fread(results_filename,
              col.names = c("plasmid_long", "position", "density"))
  dt$plasmid = gsub("_leading.fa", "", gsub(".*\\/", "", dt$plasmid_long))
  dt$PTU = plasmid_table[dt$plasmid, "PTU"]
  
  
  pos_avg <- dt[
    ,
    .(
      mean_density = mean(density),
      N = .N
    ),
    by = position
  ]
  
  setorder(pos_avg, position)
  window_bins <- 10   # e.g. 10 x 1000 bp = 10kb
  
  pos_avg[, smoothed_mean :=
            rollapply(mean_density,
                      width = window_bins,
                      FUN = mean,
                      align = "center",
                      partial = TRUE, # shorter window at edge
                      fill = NA)]
  
  pos_avg[, smoothed_N :=
            rollapply(N,
                      width = window_bins,
                      FUN = mean,
                      align = "center",
                      partial = TRUE, # shorter window at edge
                      fill = NA)]
  return(pos_avg)
}

plotDensity = function(results, x.limit=50000){
  ggplot(results, aes(position, mean_density, alpha=smoothed_N))+
    geom_line()+
    xlim(c(0,x.limit))
}

gc = fread("avg_gc.tsv",
                 col.names = c("plasmid_long", "position", "gc"))
res_4 = average_density("avg_4_step500_window5000.tsv")
res_6 = average_density("avg_6_step500_window5000.tsv")
res_8 = average_density("avg_8_step500_window5000.tsv")
res_10 = average_density("avg_10_step500_window5000.tsv")
res_12 = average_density("avg_12_step500_window5000.tsv")

p.4 = plotDensity(res_4)
p.6 = plotDensity(res_6)
p.8 = plotDensity(res_8)
p.10 = plotDensity(res_10)
p.12 = plotDensity(res_12)

# GC content
pos_avg_gc <- gc[
  ,
  .(
    mean_gc = mean(gc),
    N = .N
  ),
  by = position
]
pos_avg_gc[, smoothed_mean :=
          rollapply(mean_gc,
                    width = window_bins,
                    FUN = mean,
                    align = "center",
                    partial = TRUE, # shorter window at edge
                    fill = NA)]
p.gc = ggplot(pos_avg_gc, aes(position, mean_gc))+
  geom_line()+
  xlim(c(0,50000))

pdf('../../results/Figure-GC-with-pal-densities.pdf', width=4, height=10)
cowplot::plot_grid(p.gc+ggtitle("GC content"), 
                   p.4+theme(legend.position = "none")+ggtitle("4-bp palindromes"), 
                   p.6+theme(legend.position = "none")+ggtitle("6-bp palindromes"),
                   p.8+theme(legend.position = "none")+ggtitle("8-bp palindromes"), 
                   p.10+theme(legend.position = "none")+ggtitle("10-bp palindromes"), 
                   nrow=5)
dev.off()

gc_k4 = data.frame(gc=pos_avg_gc$mean_gc, 
                   mean_density=res_4$mean_density,
                   position=pos_avg_gc$position)
p.gc.4 = gc_k4[gc_k4$position<50000,] %>% arrange(position) %>% ggplot(aes(gc, mean_density, colour=position))+
  geom_point()+
  geom_path(size=1)+
  theme_bw()+
  theme(panel.grid = element_blank())
gc_k6 = data.frame(gc=pos_avg_gc$mean_gc, 
           mean_density=res_6$mean_density,
           position=pos_avg_gc$position,
           gc_smoothed=pos_avg_gc$smoothed_mean,
           mean_density_smoothed=res_6$smoothed_mean)
p.gc.6 = gc_k6[gc_k6$position<50000,] %>% arrange(position) %>% ggplot(aes(gc, mean_density, colour=position))+
  geom_point()+
  geom_path(size=1)+
    theme_bw()+
  theme(panel.grid = element_blank())
gc_k8 = data.frame(gc=pos_avg_gc$mean_gc, 
                   mean_density=res_8$mean_density,
                   position=pos_avg_gc$position)
p.gc.8 = gc_k8[gc_k8$position<50000,] %>% arrange(position) %>% ggplot(aes(gc, mean_density, colour=position))+
  geom_point()+
  geom_path(size=1)+
  theme_bw()+
  theme(panel.grid = element_blank())

# Checking palindrome densities without smoothing
gc = fread("avg_gc_1k_1k.tsv",
           col.names = c("plasmid_long", "position", "gc"))
# GC content
pos_avg_gc <- gc[
,
.(
mean_gc = mean(gc),
N = .N
),
by = position
]
pos_avg_gc[, smoothed_mean :=
rollapply(mean_gc,
width = window_bins,
FUN = mean,
align = "center",
partial = TRUE, # shorter window at edge
fill = NA)]
# Read in data without smoothing
res.6.non = average_density("avg_6_step1000_window1000.tsv")
# Combine data
gc_k6 = data.frame(gc=pos_avg_gc$mean_gc,
            mean_density=res.6.non$mean_density,
            position=pos_avg_gc$position,
            gc_smoothed=pos_avg_gc$smoothed_mean,
            mean_density_smoothed=res.6.non$smoothed_mean)
# Plot first 50kb
p.gc.no.smoothing = gc_k6[which(gc_k6$position<50000),] %>% arrange(position) %>% ggplot(aes(gc, mean_density, colour=position))+
        geom_point()+
        geom_path(size=1)+
        theme_bw()+
        theme(panel.grid = element_blank())
# Model and look at residuals
# theoretical expectation is
theoreticalPalindromeExpectation = function(p, k){# p is GC content, k is k-mer length
  return ((p**2 - p + 1/2)**(k/2))
}
gc_k6$theoretical_density =theoreticalPalindromeExpectation(gc_k6$gc, 6) 
#plot(lm(mean_density ~ theoretical_density, data=gc_k6[which(gc_k6$position<50000),]))
ggplot(gc_k6[which(gc_k6$position<50000),], aes(theoretical_density, mean_density, colour=position))+
         geom_point()+
         geom_path(aes(group=1))
ggplot(gc_k6[which(gc_k6$position<50000),], aes(position, mean_density-theoretical_density))+
  geom_point()+
  geom_path(aes(group=1))
  
#
# Read in data without smoothing
theoreticalDepletionPlot = function(filename, k, x.limit=c(0,50000)){ # filename is results from rust script, k is palindrome length
  obs = average_density(filename)
  # Combine data
  gc_obs = data.frame(gc=pos_avg_gc$mean_gc,
                     mean_density=obs$mean_density,
                     position=pos_avg_gc$position,
                     N=pos_avg_gc$N)
  gc_obs$theoretical_density =theoreticalPalindromeExpectation(gc_obs$gc, k) 
  return(ggplot(gc_obs, aes(position, (mean_density-theoretical_density) *1000))+
    geom_point()+
    geom_path(aes(group=1)))+
    xlim(x.limit )
}
p.4 = theoreticalDepletionPlot("avg_4_step1000_window1000.tsv", 4)
p.8 = theoreticalDepletionPlot("avg_8_step1000_window1000.tsv", 8)
p.10 = theoreticalDepletionPlot("avg_10_step1000_window1000.tsv", 10)

# Combine data
gc_k8 = data.frame(gc=pos_avg_gc$mean_gc,
                   mean_density=res.4.non$mean_density,
                   position=pos_avg_gc$position,
                   gc_smoothed=pos_avg_gc$smoothed_mean,
                   mean_density_smoothed=res.4.non$smoothed_mean)
gc_k4$theoretical_density =theoreticalPalindromeExpectation(gc_k4$gc, 4) 
ggplot(gc_k4[which(gc_k4$position<50000),], aes(position, (mean_density-theoretical_density) *1000))+
  geom_point()+
  geom_path(aes(group=1))



# Average over just one PTU
plotForPTU = function(data.table, my_ptu){
  dt_PTU = data.table[PTU==my_ptu,]
  print(nrow(dt_PTU))
  pos_avg_PTU = dt_PTU[
    ,
    .(
      mean_density = mean(density),
      N = .N
    ),
    by = position
  ]
  pos_avg_PTU[, smoothed_mean :=
               rollapply(mean_density,
                         width = window_bins,
                         FUN = mean,
                         align = "center",
                         partial = TRUE, # shorter window at edge
                         fill = NA)]
  
  pos_avg_PTU[, smoothed_N :=
               rollapply(N,
                         width = window_bins,
                         FUN = mean,
                         align = "center",
                         partial = TRUE, # shorter window at edge
                         fill = NA)]
  return(ggplot(pos_avg_PTU, aes(position, mean_density, alpha=smoothed_N))+
    geom_line()+ggtitle(PTU))
}

