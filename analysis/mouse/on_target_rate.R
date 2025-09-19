library(readr)
library(ggh4x)
library(ggplot2)
library(tidyverse)

ontarget_analysis <- read_csv("ontarget_analysis_20250912_094726.csv") 

group_human <- factor(rep(c('CM','CM','EC','EC','SMC','CF','CF','SMC'), each = 2))

group_mouse <- factor(rep(c('Brain','Lung','Spleen','Heart','Liver',
                      'Kidney','Muscle'), each = 4))

ontarget_analysis %>% 
  mutate(group = c(group_human, group_mouse),
         On_Target_Rate = Region_Alignments/Primary_Mapped_Alignments*100) %>%
  filter(group %in% c(c('Brain','Lung','Spleen','Heart',
                        'Kidney','Muscle'))) %>% 
  group_by(group, Barcode) %>%
  summarise(Total_Reads = sum(Total_Reads), 
            Primary_Mapped_Alignments = sum(Primary_Mapped_Alignments),
            Region_Alignments = sum(Region_Alignments)) %>%
  ungroup() %>%
  mutate(On_Target_Rate = Region_Alignments/Primary_Mapped_Alignments * 100) %>%
  ggplot(aes(x = group, y = On_Target_Rate, col = group)) +
  # geom_boxplot() +
  geom_point(# position = position_jitter(width = 0.05, height = 0),
    size = 4) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "On Target Rate flowcell merged\n(Dnm1l reads/total mapped reads)",
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
ggsave('plot/on_target_rate.pdf', width = 5, height = 6)

ontarget_analysis %>% 
  mutate(group = c(group_human, group_mouse),
         On_Target_Rate = Region_Alignments/Primary_Mapped_Alignments*100) %>%
  filter(group %in% c(c('Brain','Lung','Spleen','Heart',
                        'Kidney','Muscle'))) %>% 
  ggplot(aes(x = group, y = On_Target_Rate, col = group)) +
  # geom_boxplot() +
  geom_point(# position = position_jitter(width = 0.05, height = 0),
    size = 4) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "On Target Rate flowcell merged\n(Dnm1l reads/total mapped reads)",
    x = "Group",
    y = "Percentage",
    color = "Group" # Clear legend title
  ) + 
  facet_grid(.~Date) + 
  scale_y_continuous(limits = c(0, NA)) +
  # geom_hline(yintercept = 1, linetype = 'dashed', col = 'grey70') +
  force_panelsizes(rows = unit(3, 'in'),
                   cols = unit(3, 'in')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")
ggsave('plot/on_target_rate_facet.pdf', width = 10, height = 6)

ontarget_analysis %>% mutate(group = c(group_human, group_mouse),
                             On_Target_Rate = Region_Alignments/Primary_Mapped_Alignments*100) %>%
  write.csv('plot/on_target_rate.csv')
