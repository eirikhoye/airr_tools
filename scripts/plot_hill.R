library(tidyverse)

"""
Use this script to plot hill diversity and evenness profiles from 
Hill diversity and evenness profile files created with:

Hill_diversity_v2.R
"""

# Load metadata with NACT interval
metadata <- readxl::read_excel("/path/to/metadata.xlsx")

# Set path to directory of hill diversity profiles
path <- '/path/to/hill_div/'
list_of_files <- list.files(path = path,
                            recursive = TRUE,
                            pattern = "\\.tsv$",
                            full.names = TRUE)

# create dataframe and merge with metadata, then exclude samples with coverage < 5
df <- readr::read_tsv(list_of_files)
df <- df %>% left_join(metadata[, c("Sample", "COMET_ID")], by=c("sample_id"="Sample")) %>%
  left_join(coverage) %>% filter(Coverage >= 5)

# calculate confidence intervalls:
get_CI_half_width <- function(x, prob) {
  n <- length(x)
  z_t <- qt(1 - (1 - prob) / 2, df = n - 1)
  z_t * sd(x) / sqrt(n)
}

lower <- function(x, prob = 0.95) {
  mean(x) - get_CI_half_width(x, prob)
}

upper <- function(x, prob = 0.95) {
  mean(x) + get_CI_half_width(x, prob)
}

no_chemo_ <- df %>% 
  merge(metadata, by.x = 'COMET_ID', by.y = 'COMET_ID') %>% 
  as_tibble() %>% filter(group == 'no_chemo') %>%
  group_by(q) %>% summarize_at(vars(c('d','e')), funs(mean, sd, min, max, lower, upper)) %>%
  mutate(group='no_chemo')

short_chemo_ <- df %>% 
  merge(metadata, by.x = 'COMET_ID', by.y = 'COMET_ID') %>% 
  as_tibble() %>% filter(group == 'short_chemo') %>%
  group_by(q) %>% summarize_at(vars(c('d','e')), funs(mean, sd, min, max, lower, upper)) %>%
  mutate(group='short_chemo')

long_chemo_ <- df %>% 
  merge(metadata, by.x = 'COMET_ID', by.y = 'COMET_ID') %>% 
  as_tibble() %>% filter(group == 'long_chemo') %>%
  group_by(q) %>% summarize_at(vars(c('d','e')), funs(mean, sd, min, max, lower, upper)) %>%
  mutate(group='long_chemo')

df_hill <- bind_rows(bind_rows(no_chemo_, short_chemo_), long_chemo_) %>%
  mutate(group = factor(group, levels=c("no_chemo", "short_chemo", "long_chemo")))
  
# Use ggplot2 to plot hill diversity
p1 <- df_hill %>%
  mutate(`NACT-group` = case_when(
    group == 'no_chemo' ~ 'No-NACT',
    group == 'short_chemo' ~ 'Short-interval',
    group == 'long_chemo' ~ 'Long-interval'
  )) %>%
  mutate(`NACT-group` = factor(`NACT-group`, levels = c('No-NACT', 'Short-interval', 'Long-interval'))) %>%
  ggplot(aes(x=q, y=d_mean, group=`NACT-group`, color=`NACT-group`)) +
  geom_line() +
  geom_ribbon(aes(ymin = d_lower, ymax = d_upper, fill=`NACT-group`), alpha = 0.1, color=NA) +
  theme_classic() +
  xlab("q") + ylab('Hill diversity') + 
  theme(axis.title.x = element_text(face='italic')) +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + theme(legend.title = element_blank())
p1 <- set_palette(p1, 'jco')  

# Use ggplot2 to plot hill evenness
p2 <- df_hill %>%
  ggplot(aes(x=q, y=e_mean, group=group, color=group)) +
  geom_line() +
  geom_ribbon(aes(ymin = e_lower, ymax = e_upper, fill=group), alpha = 0.1, color=NA) +
  theme_classic() +
  xlab("q") + ylab('Hill evenness') + 
  theme(axis.title.x = element_text(face='italic')) +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
p2 <- set_palette(p2, 'jco')  

# merge to one plot with the same legend, then save 
leg <- get_legend(p1)
p3 <- ggarrange(p1, p2, ncol=2, nrow=1,
          legend.grob = leg, 
          legend='right')
ggsave(filename = '/path/to/hill_div_even.pdf',
       plot = p3, width = 9, height=3)
