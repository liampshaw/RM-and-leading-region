library(data.table)
library(zoo)
library(ggplot2)

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

cowplot::plot_grid(p.4, p.6, p.8, p.10, p.12, nrow=5)


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

