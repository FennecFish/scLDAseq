setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("No arguments supplied. Usage: Rscript test.R <filename>", call. = FALSE)
}

file_name <- args[1]

sim_name <- sub("\\.rds$", "", file_name)
source("doc/functions.R")
sims <- readRDS(file_name)

dat <- sc_methods(sims)
write.csv(dat, file = paste0("res/colData_", sim_name, ".csv"))

res <- sc_eval(sims, dat = dat)
rownames(res) <- sim_name
write.csv(res, file = paste0("res/adjR_silhoutte_", sim_name, ".csv"))

