####### This script is the perform permutation test on time #####
################### ad-hoc to scSTM ###########################

setwd("/proj/milovelab/wu/scLDAseq")
set.seed(1)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(mclust)
library(tibble)
library(tidyr)
library(stats)
library(data.table)

res <- read.csv("res/res_permutation.csv")
dat <- res %>%
  mutate(sig_change = ifelse(FDR < 0.05, 1, 0))

power <- dat %>%
  filter(truth == 1) %>%  # Subset where the null hypothesis is false
  group_by(level) %>%
  summarise(
    total_cases = n(),  # Total cases where truth = 1
    successful_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
    power = successful_rejections / total_cases  # Proportion of successful rejections
  ) 

typeI <- dat %>%
  filter(truth == 0) %>%  # Subset where the null hypothesis is false
  group_by(level) %>%
  summarise(
    total_cases = n(),  # Total cases where truth = 1
    false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
    alpha = false_rejections / total_cases  # Proportion of successful rejections
  ) 

par(mfrow = c(1, 2))
p <- ggplot(power, aes(x = level, y = power, group = 1)) +
  geom_line(color = "blue") +  # Connect points with lines
  geom_point(size = 3, color = "red") +  # Highlight each point
  labs(title = "Power plot under permutation test",
       x = "level",
       y = "Power") +
  theme_minimal()

t1 <- ggplot(typeI, aes(x = level, y = alpha, group = 1)) +
  geom_line(color = "blue") +  # Connect points with lines
  geom_point(size = 3, color = "red") +  # Highlight each point
  labs(title = "Type I plot under permutation test",
       x = "level",
       y = "Type I Error") +
  theme_minimal()

grid.arrange(p, t1)
  
thetaPosteriorSample <- function(model, nsims=100) {
  lambda <- model$eta
  nu <- model$nu
  out <- vector(mode="list",length=nrow(lambda)) 
  for (i in 1:length(out)) {
    sigma <- nu[[i]]
    choleskydecomp <- chol(sigma)
    mat <- rmvnorm(nsims, lambda[i,],nu[[i]],choleskydecomp)
    mat <- cbind(mat, 0)
    out[[i]] <- exp(mat - row.lse(mat))
  }
  return(out)
}

rmvnorm<-function(n,mu,Sigma,chol.Sigma=chol(Sigma)) {
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol.Sigma) +c(mu))
}

row.lse <- function(mat) {
  matrixStats::rowLogSumExps(mat)
}

permutation_test <- function(scSTMobj, bysample = FALSE, n_permutations = 1000) {
  
  K <- scSTMobj$settings$dim$K
  
  # extract metadata
  metadata <- data.frame(scSTMobj$settings$covariates$X) %>% select(-X.Intercept.)
  colnames(metadata) <- sub(".*\\.", "", colnames(metadata))
  metadata$sampleID <- scSTMobj$sampleID
  
  # extract proportion data
  prop <- data.frame(scSTMobj$theta)
  colnames(prop) <- paste0("topic", 1:K)
  metadata <- cbind(metadata, prop) %>% as.data.frame()
  
  # calculate permutation by sample
  if(bysample) {
    # Calculate the proportion difference between time for each topic and sample
    prop_diff <- metadata %>%
      pivot_longer(cols = starts_with("topic"), names_to = "topic", values_to = "value") %>%
      group_by(time, sampleID, topic) %>%
      summarise(proportion = sum(value) / n(), .groups = 'drop') %>% ## calculate topic proportion by time and sample
      arrange(sampleID, topic, time) %>%
      group_by(sampleID, topic) %>%
      summarise(
        log2_diff = log2(proportion[time == 2] / proportion[time == 1]),.groups = 'drop') %>%
      ungroup()
    
    ## Permutation test.
    unique_sampleIDs <- unique(metadata$sampleID)
    n_samples <- length(unique_sampleIDs)
    
    perm_res <- vector(mode = "list", length = length(unique(metadata$sampleID)))
    
    ## initialize bootstrapping CI 
    boot_perm <- vector(mode = "list", length = length(unique(metadata$sampleID)))
    
    
    for(i in seq_along(unique_sampleIDs)) {
      perm_res[[i]] <- matrix(NA, nrow = K, ncol = n_permutations) # Using NA to initialize
      boot_perm[[i]] <-  matrix(NA, nrow = K, ncol = n_permutations)
      rownames(boot_perm[[i]]) <- sort(colnames(prop))
    }
    
    for (i in seq_len(n_permutations)) {
  
      metadata_permuted <- metadata %>%
        group_by(sampleID) %>%
        mutate(time_permuted = sample(time)) %>%
        ungroup()
      
      perm_log2_diff <- metadata_permuted %>%
        pivot_longer(cols = starts_with("topic"), names_to = "topic", values_to = "value") %>%
        group_by(time_permuted, sampleID, topic) %>%
        summarise(proportion = sum(value) / n(), .groups = 'drop') %>% 
        arrange(sampleID, topic, time_permuted) %>%
        group_by(sampleID, topic) %>%
        summarise(
          log2_diff = log2(proportion[time_permuted == 2] / proportion[time_permuted == 1]),.groups = 'keep') %>%# log2 difference between time 
        ungroup()
      
      # bootstrap CI
      bootdata <- metadata %>% 
        pivot_longer(cols = starts_with("topic"), names_to = "topic", values_to = "value") %>%
        group_by(time, sampleID, topic) %>%
        mutate(boot_value = sample(value, size = n(), replace = TRUE)) %>%
        group_by(time, sampleID, topic) %>%
        summarise(boot_proportion = sum(boot_value) / n(), .groups = 'drop') %>% 
        arrange(topic, time) %>%
        group_by(topic, sampleID) %>%
        summarise(
          log2_diff = log2(boot_proportion[time == 2] / boot_proportion[time == 1]),.groups = 'keep') %>%# log2 difference between time 
        ungroup()
      
      # store permutated log2_diff by sample
      for (s in seq_along(unique_sampleIDs)){
        perm_sample <- perm_log2_diff %>% filter(sampleID == unique_sampleIDs[s])
        perm_res[[s]][,i] <- perm_sample[["log2_diff"]]
        
        boot_sample <- bootdata %>% filter(sampleID == unique_sampleIDs[s])
        boot_perm[[s]][, i] <- boot_sample[["log2_diff"]]
      }
      if(i %% 100 ==0) {cat("permutation ", i, "\n")}
    }
    
    res <- vector(mode = "list", length = n_samples)
    for (s in 1:n_samples) {
      res[[s]] <- data.frame()
      sample <- unique_sampleIDs[s]
      sample_dat_log2_diff <- prop_diff %>% filter(sampleID == sample)
      
      increased <- rowSums(apply(perm_res[[s]], 2, function(x) sample_dat_log2_diff[["log2_diff"]] <= x))
      increased <- (increased + 1) / (n_permutations + 1)
      
      decreased <- rowSums(apply(perm_res[[s]], 2, function(x) sample_dat_log2_diff[["log2_diff"]] >= x))
      decreased <- (decreased + 1) / (n_permutations + 1)
      
      # calcualte adjusted p-value
      sample_dat_log2_diff <- sample_dat_log2_diff %>%
        mutate(pval = ifelse(log2_diff > 0, increased[row_number()], decreased[row_number()])) %>%
        mutate(FDR = p.adjust(pval, method = "fdr"))
      
      boot_mean <- rowMeans(boot_perm[[s]], na.rm = TRUE)
      boot_ci <- as.data.frame(t(apply(boot_perm[[s]], 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))))
      colnames(boot_ci) <- c("CI_2.5", "CI_97.5")
      
      sample_dat_log2_diff <- sample_dat_log2_diff %>% 
        mutate(log2_diff_boot_mean = boot_mean,
               CI_2.5_boot = boot_ci$CI_2.5,
               CI_97.5_boot = boot_ci$CI_97.5)
      
      res[[s]] <- rbind(res[[s]], sample_dat_log2_diff)
    }
  } else { # combined sample
    
    # Calculate the proportion difference between time for each topic
    prop_diff <- metadata %>%
      pivot_longer(cols = starts_with("topic"), names_to = "topic", values_to = "value") %>%
      group_by(time, topic) %>%
      summarise(proportion = sum(value) / n(), .groups = 'drop') %>% ## calculate topic proportion by time and sample
      arrange(topic, time) %>%
      group_by(topic) %>%
      summarise(
        log2_diff = log2(proportion[time == 2] / proportion[time == 1]) # log2 difference between time 
      ) %>%
      ungroup()
    
    ## Permutation test.
    unique_sampleIDs <- unique(metadata$sampleID)
    n_samples <- length(unique_sampleIDs)
    
    perm_res <- matrix(NA, nrow = K, ncol = n_permutations) # Using NA to initialize
    
    ## initialize bootstrapping CI 
    boot_perm <- matrix(NA, nrow = K, ncol = n_permutations)
    rownames(boot_perm) <- sort(colnames(prop))
    
    for (i in seq_len(n_permutations)) {

      metadata_permuted <- metadata %>%
        group_by(sampleID) %>%
        mutate(time_permuted = sample(time)) %>%
        ungroup()
      
      perm_log2_diff <- metadata_permuted %>%
        pivot_longer(cols = starts_with("topic"), names_to = "topic", values_to = "value") %>%
        group_by(time_permuted, topic) %>%
        summarise(proportion = sum(value) / n(), .groups = 'drop') %>% 
        arrange(topic, time_permuted) %>%
        group_by(topic) %>%
        summarise(
          log2_diff = log2(proportion[time_permuted == 2] / proportion[time_permuted == 1]),.groups = 'keep') %>%# log2 difference between time 
        ungroup()
      
      perm_res[,i] <- perm_log2_diff[["log2_diff"]]
      
      bootdata <- metadata %>% 
        pivot_longer(cols = starts_with("topic"), names_to = "topic", values_to = "value") %>%
        group_by(time, topic) %>%
        mutate(boot_value = sample(value, size = n(), replace = TRUE)) %>%
        group_by(time, topic) %>%
        summarise(boot_proportion = sum(boot_value) / n(), .groups = 'drop') %>% 
        arrange(topic, time) %>%
        group_by(topic) %>%
        summarise(
          log2_diff = log2(boot_proportion[time == 2] / boot_proportion[time == 1]),.groups = 'keep') %>%# log2 difference between time 
        ungroup()
      
      boot_perm[, i] <- bootdata[["log2_diff"]]
      if(i %% 100 ==0) {cat("permutation ", i, "\n")}
    }
    

    increased <- rowSums(apply(perm_res, 2, function(x) prop_diff[["log2_diff"]] <= x))
    increased <- (increased + 1) / (n_permutations + 1)
    
    decreased <- rowSums(apply(perm_res, 2, function(x) prop_diff[["log2_diff"]] >= x))
    decreased <- (decreased + 1) / (n_permutations + 1)
    
    prop_diff <- prop_diff %>%
      mutate(pval = ifelse(log2_diff > 0, increased[row_number()], decreased[row_number()])) %>%
      mutate(FDR = p.adjust(pval, method = "fdr"))
    res <- prop_diff
    
    boot_mean <- rowMeans(boot_perm, na.rm = TRUE)
    boot_ci <- as.data.frame(t(apply(boot_perm, 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))))
    colnames(boot_ci) <- c("CI_2.5", "CI_97.5")
    
    res <- res %>% 
      mutate(log2_diff_boot_mean = boot_mean,
             CI_2.5_boot = boot_ci$CI_2.5,
             CI_97.5_boot = boot_ci$CI_97.5)
    
  }
  return(res)
}

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", pattern = "sims")
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

dat <- data.frame()
for (file_name in files){
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  seed <- sub("(.*)_.*$", "\\1", set_level)
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", file_name))
  
  tolerance <- 0.05
  truth <- colData(sims) %>%
    as.data.frame() %>%
    count(Batch, Group, time) %>%
    pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
    mutate(ratio = `1` / `2`) %>% #calcualte ratio between time for each group
    select(Batch, Group, ratio) 
  
  truth <- truth %>%
    group_by(Batch) %>%
    mutate(truth = if (all(abs(ratio - 1) <= tolerance)) 0 else 1) %>%
    ungroup() %>%
    group_by(Batch) %>%
    slice(1)
  
  stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/scSTM_combat_filterGenes_noContent_",set_level,".rds"))
  res <- permutation_test(stmobj,bysample = TRUE)
  res <- do.call(rbind,res)
  L <- sub(".*_", "", set_level)
  res$level = L
  res$seed = seed
  res$truth <- truth$truth[match(res$sampleID,truth$Batch)]
  dat <- rbind(dat,res)
}

write.csv(dat, file = "res/res_permutation.csv")


