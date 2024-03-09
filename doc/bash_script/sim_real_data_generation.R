# This script is design to generate simulation results from real data
# the parameters are estimated from parameter_estimate.R
setwd("/proj/milovelab/wu/scLDAseq/scLDAseq")

library("Seurat")
library(dplyr)
library(scuttle)
library(tidyverse)
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(mclust)
library(Seurat)
library(RaceID)
library(cidr)
library(cluster)

set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
cat(file_name)

swap = TRUE # to flip two patients

# files <- list.files(path = "data/", pattern = "BIOKEY_.*_params.rds")
# file_name <- "BIOKEY_8_sims.rds"
file_name <- basename(file_name)

samp <- sub("\\_params.rds$", "", file_name)

params <- readRDS(paste0("data/", file_name))

params <- setParams(params, group.prob = c(0.3, 0.3, 0.4),
                    de.prob = c(0.2, 0.2, 0.2), 
                    nGenes = 5000, batchCells=c(2000,2000,2000))
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.3, 0.5, 0.2))
pre_count <- cell_count %*% t(pre_prp)
colnames(pre_count) <- paste0("Group",1:length(unique(sims$Group)))
rownames(pre_count) <- paste0("Batch", 1:length(unique(sims$Batch)))

sampled_data <- colData(sims) %>%
  data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = 2) %>% 
  ungroup() 

for (i in 1:nrow(pre_count)) {
  batch_name <- rownames(pre_count)[i]
  for (j in 1:ncol(pre_count)) {
    group_name <- colnames(pre_count)[j]
    sampled_data <- sampled_data %>%
      group_by(Group, Batch) %>%
      mutate(time = ifelse(
        Group == group_name & Batch == batch_name & 
          row_number() %in% sample(row_number(), min(pre_count[i,j], n())), 1, time)) %>%
      ungroup() 
  }
}
sims$time <- sampled_data$time

cat("Generation Completed.\n")

#### QC ######
sims <- quickPerCellQC(sims)
#### feature selection #####
sims <- scuttle::logNormCounts(sims)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=2500)
sims.sub <- sims[p2.chosen,]

nsample <- length(unique(sims.sub$Batch))
ngroup <- length(unique(sims.sub$Group))

cat("QC Completed.\n")

saveRDS(sims, file = paste0("data/", samp, "_", nsample, "samples_", ngroup, "cellTypes_sims.rds"))

if (swap){
  rand.samp <- sample(unique(sims.sub$Batch), size = 2, replace = FALSE)
  sampled_data <- colData(sims.sub) %>%
    data.frame() %>%
    mutate(new_time = ifelse(time == 1 & Batch %in% rand.samp, 2, 1)) %>%
    mutate(time = ifelse(Batch %in% rand.samp, new_time, time))
  sims.sub$time <- sampled_data$time
  cat("Flip Completed.\n")
  saveRDS(sims.sub, file = paste0("data/", samp, "_flip_", nsample, "samples_", ngroup, "cellTypes_sims.rds"))
}
