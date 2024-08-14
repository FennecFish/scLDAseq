setwd("/proj/milovelab/wu/scLDAseq")
set.seed(1)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(mclust)
library(lme4)
library(lmerTest)
library(dplyr)
library(tibble)

source("doc/permutation_test.R")
source("doc/limma.R")
source("doc/LinearMixedModel.R")

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/", pattern = "3cellTypes_sims.rds")
sim_name <- sub("\\_sims.rds$", "", files)
truth <- vector(mode = "list", length = length(sim_name))
res_perm <- vector(mode = "list", length = length(sim_name))
res_lm <- vector(mode = "list", length = length(sim_name))
for (i in sim_name){
  
  dat <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/scLDAseq_",i,".rds"))
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",i,"_sims.rds"))
  truth_diff <- colData(sims) %>%
    as.data.frame() %>%
    group_by(time, Group) %>%
    summarise(count = n(), .groups = 'drop') %>% # Counts the number of rows per group
    ungroup() %>% # Removes grouping
    group_by(time) %>%
    mutate(total_count = sum(count), # Calculates total counts across all groups
           proportion = count / total_count) %>%
    group_by(Group) %>%
    summarise(
      log2_diff = log2(proportion[time == 2] / proportion[time == 1]),.groups = 'drop') %>%
    ungroup()
  truth[[i]] <- truth_diff
  perm <- permutation_test(dat)
  # limma_res <- limma_test(dat, "time")
  lmr <- lm_est(dat, "time")
  res_lm[[i]] <- lmr
  res_perm[[i]] <- perm
  rm(dat)
  rm(sims)
  cat(i, "\n")
}
