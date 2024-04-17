setwd("/proj/milovelab/wu/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(scater)
library(scran)
library(SC3)

# args <- commandArgs(trailingOnly = TRUE)
# file_name <- args[1]
# file_name <- basename(file_name)
# cat(file_name, "\n")

file_name <- "sims_1712873849_L8.rds"
set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)


sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", file_name))
rowData(sims)$feature_symbol <- rownames(sims)
# remove features with duplicated names
sims <- sims[!duplicated(rowData(sims)$feature_symbol), ]

sims <- quickPerCellQC(sims)
sims <- scuttle::logNormCounts(sims)

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

sims <- sc3(sims, ks = ngroup-1:ngroup+1, biology = TRUE)

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                                   "sc3_", set_level, ".rds"))