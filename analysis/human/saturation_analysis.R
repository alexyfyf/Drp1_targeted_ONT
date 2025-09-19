library(tidyverse)
library(data.table)
library(ggh4x)

setwd('/home/users/allstaff/yan.a/davidson_longread/yan.a/Max_2nd/fastq/merged/ncbi/dedup_bambu')

files <- list.files('saturation_results/', pattern = '*curve*', full.names = T)[1:10]

group <- c( "CERA","CERA","CERA", "LV",  "LV",  "LV",  "LV", "CL2", "CL2", "CL2")

df <- lapply(files, function(x) {
  read.csv(x) %>%
    mutate(sample = basename(x) %>% str_remove(pattern = "_saturation_curve.csv")) 
}) %>% 
  rbindlist()

df <- df %>%
  mutate(group = rep(group, each = 10),
         saturation = 1 - (unique_umi_reads/total_reads))

# impute 0 
df <- list(data.frame(fraction = 0, 
                      saturation = 0, 
                      total_reads = 0,
                      unique_umi_reads = 0,
                      sample = df$sample %>% unique(),
                      group = group) ,
           df) %>% rbindlist(use.names = T)

ggplot(df, aes(x = total_reads/1e6, y = unique_umi_reads/1e6, 
               group = sample, col = group)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Saturation analysis on chr12",
    x     = "Total reads count (millions)",
    y     = "Unique UMI count (millions)"
  ) +
  # facet_wrap(. ~ batch, scales = 'free') +
  scale_color_brewer(palette = "Set2") +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(3, 'in'))
ggsave('plot/saturation_curve.pdf', width = 5, height = 4)


ggplot(df, aes(x = total_reads/1e6, y = saturation, 
               group = sample, col = group)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Saturation analysis on chr12",
    x     = "Total reads count (millions)",
    y     = "Saturation (1-UMI/Total)"
  ) +
  geom_hline(yintercept = 1, linetype = 'dashed', col = 'grey70') +
  # facet_wrap(. ~ batch, scales = 'free') +
  scale_color_brewer(palette = "Set2") +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(3, 'in'))
ggsave('plot/saturation_curve_pct.pdf', width = 5, height = 4)
