library(tidyverse)
library(data.table)
library(ggh4x)

setwd("/vast/projects/davidson_longread/yan.a/20240326_stVincent_Max_MitocDNAPool/ncbi/dedup_bambu")

files <- list.files('saturation_results_mouse/', pattern = '*curve*', full.names = T)

group <- rep(c('Brain','Lung','Spleen','Heart','Liver','Kidney','Muscle'), each = 2)

df <- lapply(files, function(x) {
  read.csv(x) %>%
    mutate(sample = basename(x) %>% str_remove(pattern = "_saturation_curve.csv")) %>%
    separate(sample, c('batch', 'barcode'), remove = F)
}) %>% 
  rbindlist()

df <- df %>%
  mutate(group = rep(rep(group, each = 10), 2),
         saturation = 1 - (unique_umi_reads/total_reads)) %>% 
  filter(group != 'Liver')

# impute 0 
df <- list(data.frame(fraction = 0, 
           saturation = 0, 
           total_reads = 0,
           unique_umi_reads = 0,
           sample = df$sample %>% unique(),
           group = group[group!='Liver']) %>%
  separate(sample, c('batch','barcode'), remove = F),
  df) %>% rbindlist(use.names = T)

ggplot(df, aes(x = total_reads, y = unique_umi_reads, 
               group = sample, col = group)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Saturation analysis on chr16",
    x     = "Total reads count (millions)",
    y     = "Unique UMI count (millions)"
  ) +
  facet_wrap(. ~ batch, scales = 'free') +
  scale_color_brewer(palette = "Set2") +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(c(3,3), 'in'))
ggsave('plot/saturation_curve.pdf', width = 9, height = 4)

ggplot(df, aes(x = total_reads, y = saturation, 
               group = sample, col = group)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Saturation analysis on chr16",
    x     = "Total reads count (millions)",
    y     = "Saturation (1-UMI/Total)"
  ) +
  geom_hline(yintercept = 1, linetype = 'dashed', col = 'grey70') +
  facet_wrap(. ~ batch, scales = 'free') +
  scale_color_brewer(palette = "Set2") +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(c(3,3), 'in'))
ggsave('plot/saturation_curve_pct.pdf', width = 9, height = 4)
