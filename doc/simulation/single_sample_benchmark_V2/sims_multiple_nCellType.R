# this script is to generate simulated data for figure 1
# with only three levels of noise
# but different number of cell types
# this is in addition to sims.R
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library("scater")
library(SingleCellExperiment)
library(MASS)
library(VariantAnnotation)
library(checkmate)

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))

r.file <- paste0("doc/rev_splatter/",list.files("doc/rev_splatter/"))
sapply(r.file, source)

vcf <- mockVCF(n.samples = 1)
gff <- mockGFF(n.genes = 1500)

sims_scRNA <- function(nsample = 1, nCellType, de.prob, 
                       de.facLoc, seed, batch.facLoc){
  
  first.prob <- rep(1/nCellType, nCellType)
  change <- c(rep(1/nCellType, nCellType)[1]/2, -rep(1/nCellType, nCellType)[1]/2)
  change <- rep(change, each = (ceiling(nCellType/2) - 1))
  second.prob <- rep(1/nCellType, nCellType) + c(change, 0)
  
  nGenes <- 1000
  de.facScale = 0.1
  
  batchCells <- c(300, 300)
  batch.size = 1
  batch.facScale = 0.1
  
  params <- newSplatPopParams(
    # number of genes
    nGenes = nGenes,
    # number of cells & time info
    batchCells = batchCells,
    batch.size = batch.size,
    batch.facLoc = batch.facLoc,
    batch.facScale = batch.facScale,
    
    # group prob
    group.prob = list(first.prob, second.prob),
    de.prob = de.prob,
    de.facLoc = de.facLoc,
    de.facScale = de.facScale,
    
    seed = seed
  )
  
  sim <- splatPopSimulate(
    vcf = vcf,
    gff = gff,
    params = params,
    sparsify = FALSE,
    verbose = T
  )
  
  return(sim)
}

nsample <- 1
#nCellType <- c(5, 9, 13)
nCellType <- c(3, 7, 11, 15, 19, 25, 31)

nparam <- expand.grid(nsample, nCellType)
colnames(nparam) <- c("nsample", "nCellType")

for(i in 1:dim(nparam)[1]){
  
  nsample <- nparam[i,1]
  nCellType <- nparam[i,2]
  
  #### level 2 ####
  # both cell type and treatment account for variation
  de.prob <- 0.5
  de.facLoc <- 1
  batch.facLoc <- 0.2
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc)
  
  # prop.table(table(sims[[2]]$Group,sims[[2]]$Sample, sims[[2]]$Batch),
  #            margin = c(2, 3))
  # sims <- logNormCounts(sims)
  # sims <- runPCA(sims,ncomponents = 10)
  # plotPCA(sims, colour_by = "Group", shape_by = "Batch")
  # 
  saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/sims/",
                                   "sims_multiCellType_", seed, "_L2_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L2.\n")
  rm(sims)
  
  #### level 3 ####
  de.prob <- 0.3
  de.facLoc <- 1
  batch.facLoc <- 0.5
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc)
  
  # prop.table(table(sims[[2]]$Group,sims[[2]]$Sample, sims[[2]]$Batch),
  #            margin = c(2, 3))
  # sims <- logNormCounts(sims)
  # sims <- runPCA(sims,ncomponents = 10)
  # plotPCA(sims, colour_by = "Group", shape_by = "Time")
  # 
  saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/sims/",
                                   "sims_multiCellType_", seed, "_L3_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L3.\n")
  rm(sims)
  
  #### level 4 ####
  de.prob <- 0.3
  de.facLoc <- 0.5
  batch.facLoc <- 0.5
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc)

  # 
  saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/sims/",
                                   "sims_multiCellType_", seed, "_L4_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L4.\n")
  rm(sims)
}


