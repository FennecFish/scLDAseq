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
sc_eval <- function(sims, dat) {
  
  res <- data.frame()
  # compute silhouette score
  dist.matrix <- dist(t(counts(sims)))
  scSTM.sil <- silhouette(as.numeric(as.factor(dat$scSTM_cluster)), dist.matrix)
  seurat.sil <- silhouette(as.numeric(as.factor(dat$seurat_cluster)), dist.matrix)
  raceid.sil <- silhouette(as.numeric(as.factor(dat$raceID_cluster)), dist.matrix)
  # cidr.sil <- silhouette(as.numeric(as.factor(dat$cidr_cluster)), dist.matrix)
  
  res <- data.frame(
    scSTM_adjR = adjustedRandIndex(dat$scSTM_cluster,sims$Group),
    Seurat_adjR = adjustedRandIndex(dat$seurat_cluster, sims$Group), 
    raceID_adjR = adjustedRandIndex(dat$raceID_cluster,sims$Group),
    # CIDR_adjR = adjustedRandIndex(dat$cidr_cluster,sims$Group)#,
        scSTM_sil = mean(scSTM.sil[,3]),
        seurat_sil = mean(seurat.sil[,3]),
        raceID_sil = mean(raceid.sil[,3]))

  return(res)
}

files <- list.files(path = "data", pattern = "sims.rds")
sim_name <- sub("\\_sims.rds$", "", files)
res <- data.frame()

for (i in sim_name){
  dat <- read.csv(paste0("res/colData_",i,".csv"))
  sims <- readRDS(paste0("data/",i,"_sims.rds"))
  temp <- sc_eval(sims, dat = dat)
  res <- rbind(res, temp)
  rm(dat)
  rm(sims)
  cat(i)
}

rownames(res) <- sim_name
write.csv(res, file = "res/res_eval_methods.csv")
res <- tibble::rownames_to_column(res, "sim")
# change the format from wide to long
library(tidyr)
data_long <- gather(res, method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
ggplot(data_long, aes(x=method, y=adjR)) + 
  geom_boxplot()
