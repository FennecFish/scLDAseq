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
library(scuttle)
library(TopicScore)
library(sva)
library(stm)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

# file_name <- "sims_1716366084_neg_L1_c5.rds"

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/sims/", file_name))

Xtest <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/test/",
                              "test_", set_level, ".rds"))
sims <- sims[rownames(sims) %in% rownames(Xtest), colnames(sims) %in% colnames(Xtest)]
Xtrain <- counts(sims)-Xtest

train <- sims
counts(train) <- Xtrain
saveRDS(train, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/train/",
                             "train_", set_level, ".rds"))

scSTMobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/scSTM/scSTM_", set_level, ".rds"))
seed <- scSTMobj$settings$seed
initial <- scSTMobj$settings$init$mode

nsample <- length(unique(train$Batch))
ngroup <- length(unique(train$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()

scSTMobj_train <- multi_stm(train, K = ngroup, prevalence=NULL, content=NULL, sample = "Batch",
                                  init.type=initial, #seed=seed,
                                  max.em.its=20, emtol=1e-5,
                                  verbose=TRUE) 

msg <- sprintf("Completed scSTMseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(scSTMobj_train, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/scSTM_train/",
                           "scSTM_", set_level, ".rds"))