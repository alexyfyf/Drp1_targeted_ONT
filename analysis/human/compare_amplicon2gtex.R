## compare with GTeX using same Bambu pipeline
library(ggrepel)
library(ggpubr)
library(scales)
library(ggh4x)
library(tidyverse)

setwd('/home/users/allstaff/yan.a/davidson_longread/yan.a/Max_2nd/fastq/merged/ncbi/dedup_bambu')

# load human amplicon from bambu

amp <- read.table('Drp1_human_dge/filter/counts_transcript.txt', header = T)[, c(1, 6:9)]
txid <- amp$TXNAME
amp <- amp[, -1]
rownames(amp) <- txid
# amp <- amp[tx_order, ]

# load gtex bambu

gtex <- read.table('../../../../../gtex_ncbi/bambu/Drp1_human_dge/filter/counts_transcript.txt', header = T)[, -2]
txid <- gtex$TXNAME
gtex <- gtex[, -1]
rownames(gtex) <- txid
# gtex <- gtex[tx_order, ]

stopifnot(rownames(gtex) == rownames(amp))

df_plot <- data.frame(gtex_bambu = log10(rowMeans(gtex) + 1), 
                      amplicon_bambu = log10(rowMeans(amp) + 1),
                      txid = rownames(gtex)) 

p1 <- df_plot %>%
  ggplot(aes(x = gtex_bambu, y = amplicon_bambu, label = txid)) +
  geom_point(size = 2) +
  geom_text_repel(size = 3) +
  stat_cor(method = "pearson",
           label.x.npc = "left", label.y.npc = 0.90, size = 3.5) +
  labs(
    title = "Correlation (amplicon vs bulk)",
    x     = "GTEx log10(count+1)",
    y     = "Amplicon log10(count+1)"
  ) +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(3, 'in'))

df_plot2 <- data.frame(gtex_bambu = rowMeans(gtex), 
                      amplicon_bambu = rowMeans(amp),
                      txid = rownames(gtex)) 

p2 <- ggplot(df_plot2, aes(x = gtex_bambu, y = amplicon_bambu, label = txid)) +
  geom_point(size = 2) +
  geom_text_repel(size = 3) +
  stat_cor(method = "pearson",
           label.x.npc = "left", label.y.npc = 0.90, size = 3.5) +
  scale_x_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 1),
    name  = "GTEx (pseudo‑log scale)"
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 1),
    name  = "Amplicon (pseudo‑log scale)"
  ) +
  labs(
    title = "Correlation (amplicon vs bulk)"
  ) +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(3, 'in'))

cowplot::plot_grid(p1,p2)

ggsave('plot/DNM1L_gtex2amplicon_transcript_exp.pdf', width = 8, height = 5)


## compare human vs mouse amplicon for conserved transcripts

# load human amplicon from bambu

# human LV
hamp <- read.table('Drp1_human_dge/filter/counts_transcript.txt', header = T)[, c(1, 6:9)]
txid <- hamp$TXNAME
hamp <- hamp[, -1]
rownames(hamp) <- txid
# amp <- amp[tx_order, ]

# mouse heart
mamp <- read.table('~/davidson_longread/yan.a/20240326_stVincent_Max_MitocDNAPool/ncbi/dedup_bambu/Drp1_mouse_dge/merged/filter/counts_transcript.txt', header = T)[, c(1, 9,10)]
txid <- mamp$TXNAME
mamp <- mamp[, -1]
rownames(mamp) <- txid
# amp <- amp[tx_order, ]

data.frame(mouse =  log10(rowMeans(mamp[c('NM_001025947.3', #'NM_001405252.1',
                                          'NM_001360007.2','NM_001405259.1'),]) + 1),
           human =  log10(rowMeans(hamp[c('NM_005690.5', #'NM_012063.4', 
                                          'NM_001278464.2', 'NM_001278465.2'),]) + 1), 
           # mouse_id = c('b', #'d',
           #              'e','o'),
           # human_id = c('3', #'2', 
           #              '5', '8'),
           mouse_id = c('NM_001025947.3', #'NM_001405252.1',
                        'NM_001360007.2','NM_001405259.1'),
           human_id = c('NM_005690.5', #'NM_012063.4', 
                        'NM_001278464.2', 'NM_001278465.2')
           ) %>%
  unite(join_id, mouse_id:human_id, sep = ":", remove = F) %>%
  ggplot(aes(x = mouse, y = human)) + 
  geom_point() +
  # geom_smooth(method = 'lm', se = F) + 
  geom_text_repel(aes(label = join_id)) +
  stat_cor(method = "pearson",
           label.x.npc = "left", label.y.npc = 0.90, size = 3.5) +
  # stat_cor(label.x = 3)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  # stat_regline_equation(label.x = 7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + # Adds y = x line
  labs(
    title = "Selected conserved isoform",
    x = "Mouse heart tissue (log10(count+1))",
    y = "Human left ventricle tissue (log10(count+1))"
  ) +
  xlim(c(1, 6)) + ylim(c(1,6)) +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(3, 'in'))
ggsave('plot/mouse2human_exp.pdf', width = 5, height = 5)



