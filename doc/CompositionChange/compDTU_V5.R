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
# 
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/scSTM/")
# 
# all.Y <- vector(mode = "list")
# for (file_name in files){
#   scSTMobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/scSTM/", file_name))
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
# saveRDS(all.Y, file = "allY_V5_scSTM.rds")

all.Y <- readRDS("allY_V5_scSTM.rds")
power_error_plot <- function(dat, threshold = 0.05,
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

power_error_plot(res, title = "Type I error plot for CompDTU on V5 scSTM data")

ggplot(res, aes(x = pval_CompDTU)) +
  geom_histogram() +
  facet_wrap(~ gene_id) +
  theme_minimal() +
  labs(title = "Histogram of pval_CompDTU grouped by setup",
       x = "pval_CompDTU",
       y = "Count")

ggplot(res_filtered, aes(x = pval_CompDTU)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of pval_CompDTU",
       x = "pval_CompDTU",
       y = "Count")

hist(res$pval_CompDTU)

######## directly apply on sims ########
files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", pattern = "c5")
all.Y <- vector(mode = "list")
for (file_name in files){
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", file_name))
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

  dat <- colData(sims) %>%
          as.data.frame() %>%
          dplyr::count(Batch, Group, time) %>%
          dplyr::group_by(Batch, time) %>%
          dplyr::mutate(Proportion = n / sum(n)) %>%
          dplyr::select(Batch, time, Proportion) %>%
          dplyr::group_by(time, Batch) %>%
          dplyr::mutate(id = row_number()) %>%
          tidyr::pivot_wider(
            names_from = id,
            values_from = Proportion,
            names_prefix = "K_")
  all.Y[[set_level]] <- dat
  cat(file_name, "\n")
}
saveRDS(all.Y, file = "allY_V5_c5_sims.rds")

all.Y <- readRDS("allY_V5_c9_sims.rds")
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
  labs(title = "Histogram of pval for cell Type 9 on sims",
       x = "pval",
       y = "Count")

### test ###
save_ps <- length(1000)
save_ts <- length(1000)
for ( i in 1:1000 ){
  condition_A <- cbind(rnorm(10,10,0.001), rnorm(10,20,0.001), rnorm(10,30,0.001))
  condition_A <- condition_A/rowSums(condition_A)
  condition_B <- cbind(rnorm(10,10,0.001), rnorm(10,20,0.001), rnorm(10,30,0.001))
  condition_B <- condition_B/rowSums(condition_B)
  
  Y <- rbind(condition_A, condition_B)
  Group <- rep(c(1,2), each = 10)
  Group <- factor(Group)
  Y <- compositions::ilr(Y)
  test <- CompDTUReg(genename = "sample", Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL)
  save_ps[i] <- test$pval_CompDTU
  save_ts[i] <- test$FStat
}
hist(save_ps)
