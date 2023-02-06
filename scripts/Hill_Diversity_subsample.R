suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(tidyverse))

# Get paths from snakemake script
meta_path <- snakemake@input[['meta']]
file_path <- snakemake@input[['rear']]
outfile <- snakemake@output[['hill_div']]

# Read dataframes
metadata <- read_tsv(meta_path, col_types = cols())
Rearrangement <- read_tsv(file_path, col_types = cols())

# Wrangle and merge
Rearrangement <- Rearrangement %>%
dplyr::select(sample_name, amino_acid, templates, productive_frequency) %>%
dplyr::rename('sample_id'=sample_name, 'clone_id'=amino_acid,
'seq_count'=templates, 'seq_freq'=productive_frequency)

# Create Hill diversity and evenness profliles from 0 to 10, step=0.2, no resampling
hill_div_prof <- alphaDiversity(Rearrangement, group="sample_id",
                               uniform=TRUE,
                               min_q=0, max_q=10, step_q=0.2,
                               ci=0.95, nboot=200)

# write outfiles
write_tsv(hill_div_prof@diversity, file = outfile)
