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




m_exp = disexp(deg_dist[deg_dist > 0])
m_exp$setXmin(m_pl$getXmin())
est = estimate_pars(m_exp)
m_exp$setPars(est$pars)

m_ln = dislnorm(deg_dist[deg_dist > 0])
m_ln$setXmin(m_pl$getXmin())
est = estimate_pars((m_ln))
m_ln$setPars(est$pars)

m_pois = dispois(deg_dist[deg_dist > 0])
m_pois$setXmin(m_pl$getXmin())
est = estimate_pars(m_pois)
m_pois$setPars(est$pars)

comp_pl_vs_exp  = compare_distributions(m_pl, m_exp)
comp_pl_vs_ln   = compare_distributions(m_pl, m_ln)
comp_pl_vs_pois = compare_distributions(m_pl, m_pois)

stat_annot <- paste(
  'PL alpha=', round(m_pl$pars, 2), '\n',
  'PL xmin=',  m_pl$xmin, '\n',
  'PLvsLN:    R=', round(comp_pl_vs_ln$test_statistic, 1), ' p=', round(comp_pl_vs_ln$p_two_sided, 5), '\n',
  'PLvsEXP:  R=', round(comp_pl_vs_exp$test_statistic, 1), ' p=', round(comp_pl_vs_exp$p_two_sided, 5), '\n',
  'PLvsPOIS: R=', round(comp_pl_vs_pois$test_statistic, 1), ' p=', round(comp_pl_vs_pois$p_two_sided, 5),
  sep=''
  )

#print(stat_annot)














