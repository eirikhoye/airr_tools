args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
library(poweRlaw)
library(RColorBrewer)

################################################################################################

                           # Get File Paths and load degree distribution

################################################################################################

file_path <- args[1]
out_path <- args[2]

deg_dist <- read_delim(file_path, col_types = cols(), col_names=FALSE, delim = ' ')
deg_dist <- deg_dist$X2
deg_dist <- sort(deg_dist, decreasing=TRUE)

################################################################################################

                           # Estimate parameters and perform goodness of fit test

################################################################################################


m_pl = displ$new(deg_dist[deg_dist > 0])
est = estimate_xmin(m_pl)
m_pl$setXmin(est$xmin)
est = estimate_pars(m_pl)
m_pl$setPars(est$pars)

bs_p = bootstrap_p(m_pl, no_of_sims = 1000, threads=12, distance='ks', seed=123)

#print(bs_p)
write_tsv(bs_p$bootstraps, paste(out_path,'bootstraps.tsv', sep=''))

write_tsv(tibble(
  'parameter'=c('p_val', 'gof', 'sim_time', 'seed', 'distance'),
  'value' = c(bs_p$p, bs_p$gof, bs_p$sim_time, bs_p$seed, bs_p$distance)
), paste(out_path, 'parameters.tsv', sep=''))
