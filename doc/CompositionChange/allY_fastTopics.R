setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/test_fastTopics/")[-1]
dat <- data.frame()

all.Y <- vector(mode = "list")
for (file_name in files){
  ft_obj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/test_fastTopics/", file_name))
  set_level <- sub("fastTopics_([^.]*)\\.rds", "\\1",  file_name)
  
  L <- ft_obj$L
  
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/sims/sims_", set_level, ".rds"))
  time <- sims$time
  names(time) <- sims$Cell
  time <- time[names(time) %in% rownames(L)]
  
  t1 <- L[match(names(time)[time==1], rownames(L)),]
  t2 <- L[match(names(time)[time==2], rownames(L)),]
  
  sampleID <- unique(sims$Batch)
  res <- data.frame()
  for (sample in sampleID){
    x1 <- t1[rownames(t1) %in% sims$Cell[which(sims$Batch==sample)],]
    x2 <- t2[rownames(t2) %in% sims$Cell[which(sims$Batch==sample)],]
    Y <- rbind(colMeans(x1), colMeans(x2))
    res <- rbind(res, Y)
  }
  all.Y[[set_level]] <- res
  cat(file_name, "\n")
}

saveRDS(all.Y, file = "allY_neg_V3_fastTopics.rds")
