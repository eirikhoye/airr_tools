suppressPackageStartupMessages(library(tidyverse))
library(ggpubr)

metadata <- snakemake@input[['meta']]
hill_div <- snakemake@input[['in_path']]



my_comparisons <- list(c('no_chemo', 'short_chemo'),
                       c('no_chemo', 'long_chemo'),
                       c('short_chemo', 'long_chemo'))

p <- ggboxplot(Diversity_AUC, x='group', y='Diversity_AUC',
               color='group', palette='jco',
               add='jitter')

p + stat_compare_means(comparisons = my_comparisons, method='t.test')
