library(data.table)
library(ggh4x)
library(tidyverse)
library(ggridges)

setwd("/vast/projects/davidson_longread/yan.a/20240326_stVincent_Max_MitocDNAPool/ncbi/dedup_bambu")

# 1. Define the barcodes, their sample names, and their groups
group <- rep(c('Brain','Lung','Spleen','Heart','Liver','Kidney','Muscle'), each = 2)

group_map <- setNames(
  rep(group, 2),
    1:28
  )
  
sample_map <- setNames(
  c(paste(group, c('S1_0326','S2_0326'), sep = '_'),
    paste(group, c('S1_0328','S2_0328'), sep = '_')), # rename by replicate (S) 
  1:28
)

batch_map <- setNames(
  c(
    rep('20240326', 14), 
    rep('20240328', 14)),  
  1:28
)

# 2. List only the files for mouse
all_files <- c(
  list.files(
    "../../length_dist_mouse/",
    pattern    = "_length_distribution\\.txt$",
    full.names = TRUE
  ))

# 3. Read & bind, adding barcode, sample and group
df <- rbindlist(
  lapply(seq_along(all_files), function(i) {
    f  <- all_files[i]
    dt <- fread(f)
    dt[, sample  := sample_map[i]]
    dt[, group   := group_map[i]]
    dt[, batch   := batch_map[i]]
    dt
  }),
  fill = TRUE
)

df <- df[group != 'Liver']

# 4. If you want `group` as a factor in the specific order:
# df[, group := factor(group, levels = c("Cera","LV","CL2"))]

# Inspect
head(df)
df %>% group_by(group) %>% summarise(mean = mean(length))

library(dplyr)
library(ggplot2)

df %>%
  filter(length < 3000) %>%
  ggplot(aes(x = length, fill = group)) +
  geom_histogram(
    aes(y = after_stat(count / 1e6)),  # <- counts in millions
    binwidth = 100,       # adjust as needed
    position = "identity",
    alpha    = 0.4,      # make bars semi‐transparent
    color    = "black"   # optional: bar outlines
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x     = "Fragment length (<3000)",
    y     = "Count (million)",
    fill  = "Group",
    title = "Overlayed length distributions by group"
  ) +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(3, 'in'))
ggsave('plot/read_length_dist.pdf', width = 6, height = 5)

df %>%
  filter(length < 3000) %>%
  ggplot(aes(x = length, fill = group)) +
  geom_histogram(
    aes(y = after_stat(count / 1e6)),  # <- counts in millions
    binwidth = 100,       # adjust as needed
    position = "identity",
    alpha    = 0.4,      # make bars semi‐transparent
    color    = "black"   # optional: bar outlines
  ) +
  facet_wrap(group ~ ., ncol = 3) + 
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x     = "Fragment length (<3000)",
    y     = "Count (million)",
    fill  = "Group",
    title = "Length distributions by group"
  ) +
  force_panelsizes(rows = unit(c(3,3), 'in'),
                   cols = unit(c(3,3,3), 'in'))
ggsave('plot/read_length_dist_facet.pdf', width = 12, height = 8)

df %>%
  filter(length < 3000) %>%
  ggplot(aes(x = length, fill = group)) +
  geom_histogram(
    aes(y = after_stat(count / 1e6)),  # <- counts in millions
    binwidth = 100,       # adjust as needed
    position = "identity",
    alpha    = 0.4,      # make bars semi‐transparent
    color    = "black"   # optional: bar outlines
  ) +
  facet_wrap(sample ~ ., ncol = 4) + 
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x     = "Fragment length (<3000)",
    y     = "Count (million)",
    fill  = "Group",
    title = "Length distributions by sample"
  ) +
  force_panelsizes(rows = unit(c(3,3,3,3,3,3,3), 'in'),
                   cols = unit(c(3,3,3,3), 'in'))
ggsave('plot/read_length_dist_sample.pdf', width = 14, height = 25)

df %>%
  filter(length < 3000 & length > 1000) %>%
  ggplot(aes(x = length, y = group, fill = group)) +
  geom_density_ridges(
    alpha = 0.7,
    color = "black",
    scale = 1.1,
    rel_min_height = 0.001
  ) +
  theme_ridges() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Fragment length (1000-3000)",
    y = NULL,
    title = "Length distributions by group (density)"
  ) +
  force_panelsizes(rows = unit(0.5 * 6, 'in'),
                   cols = unit(3, 'in'))
ggsave('plot/read_length_dist_group_ridges.pdf', width = 6, height = 5)

df %>%
  filter(length < 3000 & length > 1000) %>%
  ggplot(aes(x = length, y = sample, fill = group)) +
  geom_density_ridges(
    alpha = 0.7,
    color = "black",
    scale = 1.1,
    rel_min_height = 0.001
  ) +
  theme_ridges() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Fragment length (1000-3000)",
    y = NULL,
    title = "Length distributions by sample (density)"
  ) +
  force_panelsizes(rows = unit(0.5 * 24, 'in'),
                   cols = unit(3, 'in'))
ggsave('plot/read_length_dist_sample_ridges.pdf', width = 7, height = 15)

library(forcats)
library(tidyverse)

# 1. Summarise: count rows per sample
counts <- df %>%
  dplyr::count(sample, group, batch, name = "reads") %>%
  arrange(batch) %>%
  mutate(sample = fct_inorder(sample))

# 2. Bar plot
ggplot(counts, aes(x = sample, y = reads/1e6, fill = group)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Number of Reads per Sample",
    x     = "Sample",
    y     = "Read Count (million)",
  ) +
  facet_grid(. ~ batch, scales = 'free_x') + 
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  force_panelsizes(rows = unit(4, 'in'),
                   cols = unit(5, 'in'))
ggsave('plot/number_reads.pdf', width = 12, height = 6)
