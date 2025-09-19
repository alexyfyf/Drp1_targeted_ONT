
library(FLAMES)
library(tidyverse)
library(ggh4x)

setwd('~/davidson_longread/yan.a/Max_2nd/fastq/merged/ncbi/dedup_bambu/')

group <- c( "CERA","CERA","CERA", "LV",  "LV",  "LV",  "LV", "CL2", "CL2", "CL2")

human_anova <- read_csv("human_anova.csv")
txid <- human_anova$...1
# gtf <- '../ref/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf'

#########
# transcriptome bam
#########

bams <- list.files('tx_bam', pattern = '_sorted.bam$', full.names = T)[1:20]

covs <- lapply(bams, function(x) {
  get_coverage(x, min_counts = 1 #, 
               # remove_UTR = T, 
               # annotation = gtf
               )
})

saveRDS(covs, 'tx_bam/tx_cov.rds')

# First block: Initial pivot and join across samples
df <- lapply(1:20, function(x) {
  covs[[x]] %>% 
    filter(transcript %in% txid) %>%
    pivot_longer(-c(transcript, read_counts, tr_length),
                             names_to = "name",      # Explicitly name this column
                             values_to = "value")    # Explicitly name this column
}) %>% reduce(full_join, by = c('transcript', 'name', 'read_counts', 'tr_length')) %>%
  mutate(position = str_remove(name, 'coverage_') %>% as.numeric()) # Removed redundant as.character()

# This assumes 'bams' is a vector of 10 elements corresponding to your 10 samples
# This renames the 'value', 'value.y', 'value.x.x' etc. columns (created by full_join)
# to the desired sample names.
# Added length(bams) for robustness if bams list size changes
colnames(df)[5:(4 + length(bams))] <- basename(bams) %>% str_remove('_sorted.bam')

# Second block: Pivot the sample columns, categorize length, and summarize
df <- df %>%
  # Use cols = c(...) instead of 5:14 for robustness, explicitly naming the columns
  # This pivots the now-named sample columns into a 'sample' column and a single 'coverage_value' column
  pivot_longer(
    cols = all_of(basename(bams) %>% str_remove('_sorted.bam')), # Recommended: use actual names for robustness
    names_to = 'sample',      # Column for sample names (e.g., 'sample1', 'sample2')
    values_to = 'coverage_value' # Column for the actual coverage values
  ) %>%
  mutate(
    # Add labels for better plot legend
    tx_length = cut(tr_length, c(0, 2000, 5000, Inf),
                labels = c("0-2000 bp", "2001-5000 bp", ">5000 bp"),
                right = TRUE, include.lowest = TRUE)
  ) %>%
  group_by(sample, tx_length, position) %>%
  summarise(mean = mean(coverage_value, na.rm = TRUE), .groups = "drop") %>% # Use coverage_value and add .groups = "drop" for good practice
  mutate(group = rep(group, each = 200)) %>%
  separate(sample, c('sample', 'method'))

# Third block: Plotting
ggplot(df %>% filter(method == 'dedup'),
       aes(x = position, y = mean, color = group)) + # <--- CRITICAL CHANGE: color by 'group'
  geom_line(aes(group = sample)) +                                   # 'group' aesthetic is implicitly handled by 'color'
  facet_grid(. ~ group) +                        # Facet by 'sample'
  labs(
    title = "Average Coverage by Group of Samples in Deduplicated data",
    x = "Relative Position",
    y = "Mean Coverage",
    color = "Group" # Clear legend title
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set2") +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(c(3,3,3), 'in'))
ggsave('plot/coverage_transcript_dedup.pdf', width = 14, height = 5)

ggplot(df %>% filter(method == 'sorted'),
       aes(x = position, y = mean, color = group)) + # <--- CRITICAL CHANGE: color by 'group'
  geom_line(aes(group = sample)) +                                   # 'group' aesthetic is implicitly handled by 'color'
  facet_grid(. ~ group) +                        # Facet by 'sample'
  labs(
    title = "Average Coverage by Group of Samples in Raw data",
    x = "Relative Position",
    y = "Mean Coverage",
    color = "Group" # Clear legend title
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set2") +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(c(3,3,3), 'in'))
ggsave('plot/coverage_transcript_raw.pdf', width = 14, height = 5)

#########
# genome bam with deeptools metagene
#########

library(data.table)
library(tidyverse)

# ---- 1. Read the data ----
# IndividualValues.tab is tab-delimited with a header, no comment lines
dat <- read_tsv('plot_coverage/coverage_matrix_dedup.tab', comment = "#", col_names = F, skip = 3)

df <- data.frame(value = dat[1,] %>% as.numeric(),
                 sample = rep(1:11, each = 60),
                 group = rep(c(group, 'other'), each = 60 ),
                 position = rep(1:60, 11))

dat2 <- read_tsv('plot_coverage/coverage_matrix_sorted.tab', comment = "#", col_names = F, skip = 3)

df2 <- data.frame(value = dat2[1,] %>% as.numeric(),
                 sample = rep(1:11, each = 60),
                 group = rep(c(group, 'other'), each = 60 ),
                 position = rep(1:60, 11))

df <- rbind(df, df2) %>%
  mutate(method = rep(c('dedup','sorted'), each = 660))

df %>% 
  filter(group != 'other', method == 'dedup') %>%
  ggplot(aes(x = position, y = value, color = group )) + # <--- CRITICAL CHANGE: color by 'group'
  geom_line(aes(group = sample)) +                                   # 'group' aesthetic is implicitly handled by 'color'
  facet_grid(. ~ group, scales = 'free_y') +
  labs(
    title = "Coverage by binned genome in Deduplicated data",
    x = "Bins",
    y = "Coverage",
    color = "Group" # Clear legend title
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set2") +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(c(3,3,3), 'in'))
ggsave('plot/coverage_genomic_dedup.pdf', width = 14, height = 5)

df %>% 
  filter(group != 'other', method == 'dedup') %>%
  ggplot(aes(x = position, y = value, color = group )) + # <--- CRITICAL CHANGE: color by 'group'
  geom_line(aes(group = sample)) +                                   # 'group' aesthetic is implicitly handled by 'color'
  facet_grid(. ~ group, scales = 'free_y') +
  labs(
    title = "Coverage by binned genome in Raw data",
    x = "Bins",
    y = "Coverage",
    color = "Group" # Clear legend title
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set2") +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(c(3,3,3), 'in'))
ggsave('plot/coverage_genomic_raw.pdf', width = 14, height = 5)


