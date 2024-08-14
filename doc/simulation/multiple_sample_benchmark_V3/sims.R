# this script is to generate simulated data for figure 1
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

sims_scRNA <- function(nsample, nCellType, de.prob, 
                       de.facLoc, seed, batch.facLoc, similarity.scale){
  
  first.prob <- rep(1/nCellType, nCellType)
  neg.prob <- rep(1/nCellType, nCellType)
  
  if(nCellType == 5){
    pos.prob <- first.prob + c(-0.15, 0.1, 0.2, -0.15, 0)
  }
  if(nCellType == 9){
    pos.prob <- first.prob + c(0.1, -0.05, -0.08, 0.03, -0.05, 0.02, 0.02, 0.01, 0)
  }
  if(nCellType == 13){
    pos.prob <- first.prob + c(0.03, -0.01, -0.02, 0.02, 0.02, -0.03, -0.01, 0.1, -0.02, -0.03, -0.01,
                               -0.04, 0)
  }

  nGenes <- 1000
  de.facScale = 0.1
  
  batchCells <- c(500, 500)
  batch.size = nsample
  batch.facScale = 0.1
  
  params.neg <- newSplatPopParams(
  similarity.scale = similarity.scale,
    # number of genes
    nGenes = nGenes,
    # number of cells & time info
    batchCells = batchCells,
    batch.size = batch.size,
    batch.facLoc = batch.facLoc,
    batch.facScale = batch.facScale,

    # group prob
    group.prob = list(first.prob,neg.prob),
    de.prob = de.prob,
    de.facLoc = de.facLoc,
    de.facScale = de.facScale,

    seed = seed
  )

  sim.neg <- splatPopSimulate(
    vcf = vcf,
    gff = gff,
    params = params.neg,
    sparsify = T,
    verbose = F
  )
  sim.neg$Time <- gsub("Batch", "Time", sim.neg$Batch)

  params.pos <- newSplatPopParams(
    similarity.scale = similarity.scale,
    # number of genes
    nGenes = nGenes,
    # number of cells & time info
    batchCells = batchCells,
    batch.size = batch.size,
    batch.facLoc = batch.facLoc,
    batch.facScale = batch.facScale,
    
    # group prob
    group.prob = list(first.prob,pos.prob),
    de.prob = de.prob,
    de.facLoc = de.facLoc,
    de.facScale = de.facScale,
    
    seed = seed
  )
  
  sim.pos <- splatPopSimulate(
    vcf = vcf,
    gff = gff,
    params = params.pos,
    sparsify = T,
    verbose = F
  )
  sim.pos$Time <- gsub("Batch", "Time", sim.pos$Batch)
  # return(sim.pos)
 return(list(sim.neg, sim.pos))
  # df <- colData(s) %>% 
  #   as.data.frame() %>% 
  #   count(Batch, Group, time) %>%
  #   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0))
  # total_1 <- sum(df$`1`)
  # total_2 <- sum(df$`2`)
  # result <- df %>%
  #   mutate(
  #     Proportion_1 = `1` / total_1,
  #     Proportion_2 = `2` / total_2,
  #     Ratio = Proportion_1 / Proportion_2
  #   )
}

nsample <- c(3, 6, 12)
nCellType <- c(5, 9, 13)

nparam <- expand.grid(nsample, nCellType)
colnames(nparam) <- c("nsample", "nCellType")

for(i in 1:dim(nparam)[1]){
  
  nsample <- nparam[i,1]
  nCellType <- nparam[i,2]
  
  vcf <- mockVCF(n.samples = nsample)
  gff <- mockGFF(n.genes = 1000)
  
  # #### level 1 ####
  # # cell Type  variaiton dominates
  # de.prob <- 0.5
  # de.facLoc <- 1
  # batch.facLoc <- 0.01
  # similarity.scale <- 10
  # sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc, similarity.scale)
  # 
  # # prop.table(table(sims[[2]]$Group,sims[[2]]$Sample, sims[[2]]$Batch),
  # #            margin = c(2, 3))
  # # sims <- logNormCounts(sims[[1]])
  # # sims <- runPCA(sims,ncomponents = 10)
  # # plotPCA(sims, colour_by = "Group", shape_by = "Time")
  # # plotPCA(sims, colour_by = "Sample", shape_by = "Time")
  # # 
  # saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
  #                                  "sims_", seed, "_neg_L1_c", nCellType, "_nsample", nsample, ".rds"))
  # cat("Generation Completed for L1-neg.\n")
  # saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
  #                                  "sims_", seed, "_pos_L1_c", nCellType, "_nsample", nsample, ".rds"))
  # cat("Generation Completed for L1-pos.\n")
  # rm(sims)
  # 
  # #### level 2 ####
  # # same setting with increased sample variation
  # de.prob <- 0.5
  # de.facLoc <- 1
  # batch.facLoc <- 0.01
  # similarity.scale <- 6
  # sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc, similarity.scale)
  # 
  # # prop.table(table(sims$Group,sims$Sample, sims$Batch),
  # #            margin = c(2, 3))
  # # sims <- logNormCounts(sims[[1]])
  # # sims <- runPCA(sims,ncomponents = 10)
  # # plotPCA(sims, colour_by = "Group", shape_by = "Sample")
  # # plotPCA(sims, colour_by = "Group", shape_by = "Time")
  # 
  # saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
  #                                  "sims_", seed, "_neg_L2_c", nCellType, "_nsample", nsample, ".rds"))
  # cat("Generation Completed for L2-neg.\n")
  # saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
  #                                  "sims_", seed, "_pos_L2_c", nCellType, "_nsample", nsample, ".rds"))
  # cat("Generation Completed for L2-pos.\n")
  # rm(sims)
  # 
  # #### level 3 ####
  # # Sample variation dominant compared to group and treatment, increased treatment difference
  # de.prob <- 0.5
  # de.facLoc <- 1
  # batch.facLoc <- 0.01
  # similarity.scale <- 2
  # sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc, similarity.scale)
  # 
  # # prop.table(table(sims[[2]]$Group,sims[[2]]$Sample, sims[[2]]$Batch),
  # # #            margin = c(2, 3))
  # # sims <- logNormCounts(sims)
  # # sims <- runPCA(sims,ncomponents = 10)
  # # plotPCA(sims, colour_by = "Group", shape_by = "Sample")
  # # plotPCA(sims, colour_by = "Sample", shape_by = "Time")
  # 
  # saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
  #                                  "sims_", seed, "_neg_L3_c", nCellType, "_nsample", nsample, ".rds"))
  # cat("Generation Completed for L3-neg.\n")
  # saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
  #                                  "sims_", seed, "_pos_L3_c", nCellType, "_nsample", nsample, ".rds"))
  # cat("Generation Completed for L3-pos.\n")
  # rm(sims)
  
  #### level 4 ####
  # both Sample variation and treatment variance are stronger than cell type 
  de.prob <- 0.2
  de.facLoc <- 0.5
  batch.facLoc <- 0.1
  similarity.scale <- 6
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc, similarity.scale)
  
  # sims <- logNormCounts(sims[[1]])
  # sims <- runPCA(sims,ncomponents = 10)
  # plotPCA(sims, colour_by = "Sample", shape_by = "Group")
  # plotPCA(sims, colour_by = "Time", shape_by = "Group")
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
                                   "sims_", seed, "_neg_L4_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L4-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
                                   "sims_", seed, "_pos_L4_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L4-pos.\n")
  rm(sims)
  
  #### level 5 ####
  # treatment factors in a very strong variation
  de.prob <- 0.2
  de.facLoc <- 0.5
  batch.facLoc <- 1
  similarity.scale <- 6
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc, similarity.scale)
  
  # sims <- logNormCounts(sims)
  # umap <- calculateUMAP(sims)
  # umap$Sample <- sims$Sample
  # umap$Time <- sims$Time
  # umap$Group <- sims$Group
  # 
  # sims <- runPCA(sims,ncomponents = 10)
  # plotPCA(sims, colour_by = "Time", shape_by = "Group")
  # plotPCA(sims, colour_by = "Sample", shape_by = "Time")
  # 
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
                                   "sims_", seed, "_neg_L5_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L5-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
                                   "sims_", seed, "_pos_L5_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L5-pos.\n")
  rm(sims)
  
  #### level 6 ####
  # treatment factors in a very strong variation
  de.prob <- 0.2
  de.facLoc <- 0.5
  batch.facLoc <- 2
  similarity.scale <- 6
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, seed, batch.facLoc, similarity.scale)
  
  # sims <- logNormCounts(sims)
  # umap <- calculateUMAP(sims)
  # umap$Sample <- sims$Sample
  # umap$Time <- sims$Time
  # umap$Group <- sims$Group
  # 
  # sims <- runPCA(sims,ncomponents = 10)
  # plotPCA(sims, colour_by = "Time", shape_by = "Group")
  # plotPCA(sims, colour_by = "Sample", shape_by = "Time")
  # 
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
                                   "sims_", seed, "_neg_L6_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L5-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
                                   "sims_", seed, "_pos_L6_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L6-pos.\n")
  rm(sims)
}


