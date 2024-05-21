setwd("/proj/milovelab/wu/scLDAseq")
library(compositions)
library(Matrix)
library(dplyr)
library(tidyverse)
library(Rcpp)
library(SingleCellExperiment)
library(CompDTUReg)
library(stats)
library("gridExtra")
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/scSTM_combat_f_nc/")
# 
# all.Y <- vector(mode = "list")
# for (file_name in files){
#   scSTMobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/scSTM_combat_f_nc/", file_name))
#   set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file_name)
# 
#   time <- scSTMobj$settings$covariates$X[,-1]
#   names(time) <- scSTMobj$DocName
# 
#   theta <- scSTMobj$theta
#   rownames(theta) <- scSTMobj$DocName
# 
#   t1 <- theta[match(names(time)[time==1], rownames(theta)),]
#   t2 <- theta[match(names(time)[time==2], rownames(theta)),]
# 
#   names(scSTMobj$sampleID) <- scSTMobj$DocName
#   sampleID <- unique(scSTMobj$sampleID)
# 
#   # the following code is to calculate the composition mean for each patient
#   res <- data.frame()
#   for (sample in sampleID){
#     x1 <- t1[rownames(t1) %in% names(which(scSTMobj$sampleID==sample)),]
#     x2 <- t2[rownames(t2) %in% names(which(scSTMobj$sampleID==sample)),]
#     Y <- rbind(colMeans(x1),colMeans(x2))
#     cat(file_name, "\n")
#     res <- rbind(res, Y)
#   }
#   all.Y[[set_level]] <- res
# }
# saveRDS(all.Y, file = "allY_V4_scSTM.rds")

all.Y <- readRDS("allY_V4_scSTM.rds")
power_error_plot <- function(dat, threshold = 0.05,
                             title1, title2){
  dat <- dat %>%
    # mutate(padj = p.adjust(pval_CompDTU, method = "fdr")) %>% # not using adjusted pvalue
    mutate(sig_change = ifelse(pval_CompDTU < threshold, 1, 0))

  power <- dat %>%
    filter(truth == 1) %>%  # Subset where the null hypothesis is false
    group_by(gene_id) %>%
    summarise(
      total_cases = n(),  # Total cases where truth = 1
      successful_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
      power = successful_rejections / total_cases  # Proportion of successful rejections
    )

  typeI <- dat %>%
    filter(truth == 0) %>%  # Subset where the null hypothesis is false
    group_by(gene_id) %>%
    summarise(
      total_cases = n(),  # Total cases where truth = 1
      false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
      alpha = false_rejections / total_cases  # Proportion of successful rejections
    )

  par(mfrow = c(1, 2))
  p <- ggplot(power, aes(x = gene_id, y = power, group = 1)) +
    geom_line(color = "blue") +  # Connect points with lines
    geom_point(size = 3, color = "red") +  # Highlight each point
    labs(title = title1,
         x = "Setup",
         y = "Power") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  t1 <- ggplot(typeI, aes(x = gene_id, y = alpha, group = 1)) +
    geom_line(color = "blue") +  # Connect points with lines
    geom_point(size = 3, color = "red") +  # Highlight each point
    labs(title = title2,
         x = "Setup",
         y = "Type I Error") +
    theme_minimal()

  return(grid.arrange(p, t1))
}

##### CompDTU ########
compDTU <- function(Y){
  # Y <- Y[,2:ncol(Y)]
  # browser()
  Y <- compositions::ilr(Y)
  Group <- rep(c(1,2), times = nrow(Y)/2)
  Group <- factor(Group)
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

power_error_plot(res, title1 = "V4 Power Plot for CompDTU", title2 = "V4 Type I error plot for CompDTU")

