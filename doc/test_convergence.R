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

sims <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/sims_1712873779_L3.rds")
### remove genes with count 0 
sims <- sims[rowSums(counts(sims)) != 0,]


nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()

scSTM.mod <- selectModel(sce = sims,
                         K = ngroup, prevalence = ~time, content = NULL,
                         sample = "Batch", N = 5, ts_runs = 10, random_run = 20)
# sims <- quickPerCellQC(sims, filter=TRUE)
#  plotColData(sims.qc, x = "sum", y="detected") 


