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