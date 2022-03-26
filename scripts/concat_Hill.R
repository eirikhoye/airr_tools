suppressPackageStartupMessages(library(tidyverse))

#set paths
hill_out <- snakemake@output[['hill_out']]
auc_out <- snakemake@output[['auc_out']]
file_paths <- snakemake@input[['in_path']]

### CONCATENATE ALL HILL DIVERSITIES ###
df <- file_paths %>%
    map(read_tsv, col_types=cols()) %>%
    reduce(rbind)
write_tsv(df, file = hill_out)

#### MAKE AUC CURVES ###
library(Bolstad2)
samples <- unique(df$sample_id)
auc_df <- matrix(ncol=1, nrow=0)

for (i in samples){
  x <- df %>% filter(sample_id == i) %>% pull(q)
  y <- df %>% filter(sample_id == i) %>% pull(e)
  auc <- sintegral(x, y)$int
  auc_df[i] <- auc
}

auc_df <- data.frame(auc_df)
auc_df['sample_id'] <- rownames(auc_df)
auc_df['clonality'] = 10 - auc_df['auc_df']
auc_df <- auc_df[, c("sample_id", "auc_df", "clonality")]
colnames(auc_df) <- c("sample_id", "auc_evenness", "clonality")
auc_df <- as_tibble(auc_df)
write_tsv(auc_df, file=auc_out)
