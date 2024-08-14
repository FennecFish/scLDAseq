setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(splatter)
library(SingleCellExperiment)
library(Seurat)
library(sctransform)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)


sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sims/", file_name))
seurat.sims <- CreateSeuratObject(counts = counts(sims))
seurat.sims <- SCTransform(seurat.sims, verbose = TRUE, variable.features.n = 800)

seurat.sims <- RunPCA(seurat.sims, verbose = TRUE)
seurat.sims <- FindNeighbors(seurat.sims, dims = 1:30, verbose = TRUE)
seurat.sims <- FindClusters(seurat.sims, verbose = TRUE, resolution = seq(0.5,2,0.5))

saveRDS(seurat.sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sctransform/", 
                                   "sctransform_", set_level, ".rds"))