setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(compositions)
library(dplyr)
library(tidyverse)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(CompDTUReg)
library(stats)

all.Y <- readRDS("allY_V2_sims.rds")
compDTU <- function(Y){
  Group <- Y$time
  Group <- factor(Group)
  Y <- Y[,3:ncol(Y)]
  # browser()
  Y <- compositions::ilr(Y)
  res <- data.frame()
  temp <- CompDTUReg(genename = "sample", Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL)
  res <- rbind(res, temp)
  return(res)
}

res <- lapply(all.Y, FUN = compDTU)
res <- Map(function(df, name) {
  df$seed <- str_split(name, "_", simplify = TRUE)[1]
  df$control <- str_split(name, "_", simplify = TRUE)[2]
  df$level <- str_split(name, "_", simplify = TRUE)[3]
  df$nCellType <- str_split(name, "_", simplify = TRUE)[4]
  df
}, res, names(res))

res <- do.call(rbind, res) %>%
  as.data.frame() %>%
  mutate(truth = ifelse(control == "neg", 0, 1),
         gene_id = paste0(control, "_", level, "_", nCellType))

ggplot(res, aes(x = pval_CompDTU)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  facet_wrap(~nCellType) +
  labs(title = "Histogram of pval directly on true simulated data",
       x = "pval",
       y = "Count")

dat <- res %>%
  # mutate(padj = p.adjust(pval_CompDTU, method = "fdr")) %>% # not using adjusted pvalue
  mutate(sig_change = ifelse(pval_CompDTU < 0.05, 1, 0))

typeI <- dat %>%
  filter(truth == 0) %>%  # Subset where the null hypothesis is false
  group_by(gene_id) %>%
  summarise(
    total_cases = n(),  # Total cases where truth = 1
    false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
    alpha = false_rejections / total_cases  # Proportion of successful rejections
  )

typeIerror_plot <- function(dat, threshold = 0.05,
                            title){
  dat <- dat %>%
    # mutate(padj = p.adjust(pval_CompDTU, method = "fdr")) %>% # not using adjusted pvalue
    mutate(sig_change = ifelse(pval_CompDTU < threshold, 1, 0))
  
  typeI <- dat %>%
    filter(truth == 0) %>%  # Subset where the null hypothesis is false
    group_by(gene_id) %>%
    summarise(
      total_cases = n(),  # Total cases where truth = 1
      false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
      alpha = false_rejections / total_cases  # Proportion of successful rejections
    )
  
  t1 <- ggplot(typeI, aes(x = gene_id, y = alpha, group = 1)) +
    geom_line(color = "blue") +  # Connect points with lines
    geom_point(size = 3, color = "red") +  # Highlight each point
    labs(title = title,
         x = "id",
         y = "Type I Error") +
    theme_minimal()
  
  return(t1)
}
typeIerror_plot(res, title = "Type I error plot for CompDTU")
