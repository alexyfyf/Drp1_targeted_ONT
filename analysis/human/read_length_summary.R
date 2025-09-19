library(data.table)
library(ggh4x)
library(tidyverse)
library(ggridges)

setwd('/home/users/allstaff/yan.a/davidson_longread/yan.a/Max_2nd/fastq/merged/ncbi/dedup_bambu')
# 1. Define the barcodes, their sample names, and their groups

sample_map <- setNames(
  c(
    "Cera_A", "Cera_B", "Cera_C", "LV_3160","LV5138","LV6086","LV3149","CM1",   "CM2",   "CM4"), # batch 2
  1:10
)
group_map <- setNames(
  c(
    "CERA","CERA","CERA", "LV",  "LV",  "LV",  "LV", "CL2", "CL2", "CL2"),   # batch 2
  1:10
)

batch_map <- setNames(
  c(
    rep('batch2', 10)),   # batch 2 human
  1:10
)

# 2. List only the files for barcode01–barcode10
all_files <- c(
  list.files(
  "../../length_dist/",
  pattern    = "_length_distribution\\.txt$",
  full.names = TRUE
)[1:10])

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
    y     = "Count (millions)",
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
    y     = "Count (millions)",
    fill  = "Group",
    title = "Length distributions by group"
  ) +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(c(3,3,3), 'in'))
ggsave('plot/read_length_dist_facet.pdf', width = 14, height = 5)

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
  facet_wrap(sample ~ ., nrow = 4) + 
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x     = "Fragment length (<3000)",
    y     = "Count (millions)",
    fill  = "Group",
    title = "Length distributions by sample"
  ) +
  force_panelsizes(rows = unit(c(3,3,3), 'in'),
                   cols = unit(c(3,3,3,3), 'in'))
ggsave('plot/read_length_dist_sample.pdf', width = 14, height = 16)

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
  force_panelsizes(rows = unit(1.5, 'in'),
                   cols = unit(3, 'in'))
ggsave('plot/read_length_dist_group_ridges.pdf', width = 5, height = 3)

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
  force_panelsizes(rows = unit(5, 'in'),
                   cols = unit(3, 'in'))
ggsave('plot/read_length_dist_sample_ridges.pdf', width = 5, height = 8)

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
    y     = "Read Count (millions)"
  ) +
  # facet_grid(. ~ batch, scales = 'free_x') + 
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  force_panelsizes(rows = unit(4, 'in'),
                   cols = unit(5, 'in'))
ggsave('plot/number_reads.pdf', width = 7, height = 6)



