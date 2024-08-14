setwd("/proj/milovelab/wu/scLDAseq")
library(ggplot2)
library(dplyr)

########### adjusted Rand Index Plot for Single Patient Simulation ###############
dat <- read.csv("res/adjRandIndex_single_sample_benchmark.csv")

dat_long <- dat %>% 
  dplyr::select(-sizeFactor) %>%
  gather(method, adjR, scSTM_nC_P_cluster:monocle3_cluster, factor_key=TRUE) %>%
  drop_na()  %>%
  group_by(nCellType, level, control, method) %>%
  summarise(mean_adjR = mean(adjR, na.rm = TRUE)) %>%
  mutate(line_type = ifelse(grepl("scSTM", method), "solid", "dashed"))
  
# scSTM <- dat_long[grep("scSTM", dat_long$method),]
ggplot(dat_long, aes(x = level, y = mean_adjR, color = method, group = method, linetype = line_type)) +
  geom_line() +
  facet_wrap(control ~ nCellType) +
  labs(x = "Level", y = "Mean Adjusted R", title = "Mean Adjusted R by Level and Method, Faceted by Cell Type") +
  theme_minimal()

# first demonstrate different version of scSTM
res_stats <- res_stats %>% #filter(grepl("^scSTM", method)) %>%
  mutate(line_type = ifelse(grepl("stm_cluster", method), "solid", "dashed")) %>%
  mutate(control = sapply(strsplit(level, "_"), `[`, 1),
         level = sapply(strsplit(level, "_"), `[`, 2))

# res_stats$control <- factor(res_stats$control)
png("res/figure1_scSTM_adjR.png", width = 2500, height = 1500, res = 300)
ggplot(dat_long, aes(x = level, y = mean_adjR, color = method, group = method, linetype = line_type)) +
  geom_line(linewidth = 1) +
  labs(title = "Mean Adjusted Rand Index by scSTM Methods and Level of Noise",
       x = "Level of Noise",
       y = "Mean adjRandIndex") +
  facet_grid(~control) +
  #  geom_errorbar(aes(ymin = mean_adjR - sd, ymax = mean_adjR + sd), width = 0.2) +
  theme_minimal() +
  # scale_color_brewer(palette = "Set5") +
  # scale_color_viridis_d() + 
  scale_color_brewer(palette = "Set2") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
  theme(axis.text.x = element_text(hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),  # Increase size of y-axis text
        axis.title.x = element_text(size = 12, face = "bold"),  # Increase size of x-axis title
        axis.title.y = element_text(size = 12, face = "bold"),  # Increase size of y-axis title
        plot.title = element_text(size = 12, face = "bold"),  # Increase size and bold the plot title
        strip.text = element_text(size = 12, face = "bold"))
dev.off()
