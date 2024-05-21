setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(countsplit)
library(tidyverse)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(scater)
library(scran)
library(sva)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", file_name))
# quick qc
sims <- quickPerCellQC(sims, filter=TRUE)
### remove genes with count 0
sims <- sims[rowSums(counts(sims)) != 0,]
#### feature selection #####
sims <- scuttle::logNormCounts(sims)

dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=2000)
sims <- sims[p2.chosen,]


# batch effect removal using combat seq
adjusted_counts <- ComBat_seq(counts(sims), batch=sims$Batch, group=NULL)
counts(sims) <- adjusted_counts
cat("Batch Effect Removed \n")

###### split count ############
split <- countsplit(adjusted_counts)
Xtrain <- split[[1]]
Xtest <- split[[2]]

saveRDS(Xtest, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/test/",
                             "test_", set_level, ".rds"))
