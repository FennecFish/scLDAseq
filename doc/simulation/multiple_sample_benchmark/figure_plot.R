setwd("/proj/milovelab/wu/scLDAseq")
library(ggplot2)
library(dplyr)

########### adjusted Rand Index Plot for Single Patient Simulation ###############
dat <- read.csv("res/adjRandIndex_multiple_sample_benchmark.csv")

dat_long <- dat %>% 
  #dplyr::select(-sctransform_cluster) %>%
  gather(method, adjR, scSTM_nC_P_cluster:scSTM_C_P_cluster, factor_key=TRUE) %>%
  drop_na()  %>%
  group_by(nCellType, level, control, method) %>%
  summarise(mean_adjR = mean(adjR, na.rm = TRUE)) %>%
  mutate(line_type = ifelse(grepl("scSTM", method), "solid", "dashed"))
  
dat_long$nCellType <- factor(dat_long$nCellType, levels = c("c5", "c9", "c13"))

#png("res/single_patient_simulation_V2/adjustedRandIndex_single_sample_benchmark_V2_figure.png", width = 2500, height = 2500, res = 300)
ggplot(dat_long, aes(x = level, y = mean_adjR, color = method, group = method)) +
  geom_line(aes(linetype = line_type)) +
  facet_wrap(~ control + nCellType) +
  labs(
    x = "Noise Level",
    y = "Mean Adjusted Rand Index",
    title = "Benchmarking Clustering Accuracy Using Mean Adjusted Rand Index"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  guides(linetype = "none") # Remove linetype from the legend
#dev.off()
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
