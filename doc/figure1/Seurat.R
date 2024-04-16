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
library(scater)
library(scran)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)


sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", file_name))

sims <- scuttle::logNormCounts(sims)
seurat.sims <- as.Seurat(sims, counts = "counts", data = "logcounts")
# Visualize QC metrics as a violin plot
# VlnPlot(seurat.sims, features = c("nFeature_originalexp", "nCount_originalexp"), ncol = 2)
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 500)
all.genes <- rownames(seurat.sims)
seurat.sims <- ScaleData(seurat.sims, features = all.genes)
seurat.sims <- RunPCA(seurat.sims, , features = VariableFeatures(object = seurat.sims))

seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)

saveRDS(seurat.sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                           "seurat_", set_level, ".rds"))