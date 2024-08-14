setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(MASS)
library(Rcpp)
library(scran)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(doParallel)
library(foreach)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sims/", pattern = "sims*")
file_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sims/",files[1])
sce <- readRDS(file_name)
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")
scSTM.mod <- selectModel_parallel(sce = sce,
                                  K = length(unique(sce$Group)), 
                                  prevalence = ~Batch, content = ~Batch,
                                  N = 3, ts_runs = 5, random_run = 5,
                                  max.em.its = 100 , sample = "Sample", emtol=1e-6,
                                  gc = 5)
saveRDS(scSTM.mod, file = "parallel_test_res.rds")