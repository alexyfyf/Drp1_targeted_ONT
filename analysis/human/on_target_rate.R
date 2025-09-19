library(readr)
library(ggh4x)
library(ggplot2)
library(tidyverse)

ontarget_analysis <- read_csv("ontarget_analysis_20250912_094606.csv")[1:10,]

group <- factor(c('CERA','CERA','CERA','LV','LV','LV','LV','CL2','CL2','CL2'))

ontarget_analysis %>% 
  mutate(group = group,
         On_Target_Rate = Region_Alignments/Primary_Mapped_Alignments*100) %>%
  ggplot(aes(x = group, y = On_Target_Rate, col = group)) +
  # geom_boxplot() +
  geom_point(# position = position_jitter(width = 0.05, height = 0),
    size = 4) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "On Target Rate\n(DNM1L reads/total mapped reads)",
    x = "Group",
    y = "Percentage",
    color = "Group" # Clear legend title
  ) + 
  scale_y_continuous(limits = c(0, NA)) +
  # geom_hline(yintercept = 1, linetype = 'dashed', col = 'grey70') +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(3, 'in')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")
ggsave('plot/on_target_rate.pdf', width = 5, height = 5)

ontarget_analysis %>% 
  mutate(group = group,
         On_Target_Rate = Region_Alignments/Primary_Mapped_Alignments*100) %>%
  write.csv('plot/on_target_rate.csv')
