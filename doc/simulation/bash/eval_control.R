setwd("/proj/milovelab/wu/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(mclust)
library(Seurat)
library(RaceID)
library(cidr)
library(cluster)

set.seed(1)

cat("start loading args \n")
args <- commandArgs(trailingOnly = TRUE)

file_name <- args[1]
cat(file_name)

# files <- list.files(path = "data/", pattern = "BIOKEY_.*_sims.rds")
swap = FALSE
# file_name <- "BIOKEY_8_sims.rds"
file_name <- basename(file_name)
sim_name <- sub("\\.rds$", "", file_name)
source("doc/functions.R")

cat(paste0("data/",file_name))
sims <- readRDS(paste0("data/control/",file_name))
# 
# #### QC ######
# sims <- quickPerCellQC(sims)
# #### feature selection #####
# sims <- scuttle::logNormCounts(sims)
# dec.p2 <- modelGeneVar(sims)
# # feature selection
# p2.chosen <- getTopHVGs(dec.p2, n=2000)
# sims <- sims[p2.chosen,]

dat <- sc_methods(sims)
write.csv(dat, file = paste0("res/colData_", sim_name, ".csv"))

if (swap) {
  samp <- sub("\\_sims$", "", sim_name)
  sampled_data <- colData(sims) %>%
    data.frame() %>%
    mutate(new_time = ifelse(time == 1 & Batch %in% c("Batch1", "Batch2"), 2, 1)) %>%
    mutate(time = ifelse(Batch %in% c("Batch1", "Batch2"), new_time, time))
  sims$time <- sampled_data$time
  samp <- paste0(samp,"_oppo")
  saveRDS(sims, paste0("data/",samp,"_sims.rds"))
  dat <- sc_methods(sims)
  write.csv(dat, file = paste0("res/colData_", sim_name, ".csv"))
}


#### evaluating all methods #####
source("doc/functions.R")
files <- list.files(path = "data/control", pattern = "rds")
sim_name <- sub("\\.rds$", "", files)
res <- data.frame()

for (i in sim_name){
  dat <- read.csv(paste0("res/colData_",i,".csv"))
  sims <- readRDS(paste0("data/control/",i,".rds"))
  temp <- sc_eval(sims, dat = dat)
  res <- rbind(res, temp)
  rm(dat)
  rm(sims)
  cat(i)
}

rownames(res) <- sim_name
write.csv(res, file = "res/adjR_silhoutte_controls.csv")

