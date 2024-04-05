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
library(cluster)

set.seed(1)

cat("start loading args \n")
args <- commandArgs(trailingOnly = TRUE)

file_name <- args[1]
cat(file_name, "\n")

file_name <- basename(file_name)
sim_name <- sub("\\_sims.rds$", "", file_name)

cat(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",file_name), "\n")
# file_name <- "BIOKEY_22_flip_3samples_3cellTypes_sims.rds"
sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",file_name))

source("doc/scLDAseq_functions.R")
dat <- scLDAseq(sims)
saveRDS(dat, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/scLDAseq_", sim_name, ".rds"))

