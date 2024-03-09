setwd("/proj/milovelab/wu/scLDAseq/scLDAseq")
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

cat("start loading args \n")
args <- commandArgs(trailingOnly = TRUE)

file_name <- args[1]
cat(file_name)

# files <- list.files(path = "data/", pattern = "BIOKEY_.*_sims.rds")
swap = FALSE
# file_name <- "BIOKEY_8_sims.rds"
file_name <- basename(file_name)
sim_name <- sub("\\.rds$", "", file_name)
source("doc/functions.R")

cat(paste0("data/",file_name))
sims <- readRDS(paste0("data/control/",file_name))
# 
# #### QC ######
# sims <- quickPerCellQC(sims)
# #### feature selection #####
# sims <- scuttle::logNormCounts(sims)
# dec.p2 <- modelGeneVar(sims)
# # feature selection
# p2.chosen <- getTopHVGs(dec.p2, n=2000)
# sims <- sims[p2.chosen,]

dat <- sc_methods(sims)
write.csv(dat, file = paste0("res/colData_", sim_name, ".csv"))

if (swap) {
  samp <- sub("\\_sims$", "", sim_name)
  sampled_data <- colData(sims) %>%
    data.frame() %>%
    mutate(new_time = ifelse(time == 1 & Batch %in% c("Batch1", "Batch2"), 2, 1)) %>%
    mutate(time = ifelse(Batch %in% c("Batch1", "Batch2"), new_time, time))
  sims$time <- sampled_data$time
  samp <- paste0(samp,"_oppo")
  saveRDS(sims, paste0("data/",samp,"_sims.rds"))
  dat <- sc_methods(sims)
  write.csv(dat, file = paste0("res/colData_", sim_name, ".csv"))
}



# 
# res <- sc_eval(sims, dat = dat)
# rownames(res) <- sim_name
# write.csv(res, file = paste0("res/adjR_silhoutte_", sim_name, ".csv"))
# 
# dat <-read.csv("res/colData_BIOKEY_14_sims.csv")
# sims <- readRDS("data/BIOKEY_14_sims.rds")
# #### QC ######
# sims <- quickPerCellQC(sims)
# #### feature selection #####
# sims <- scuttle::logNormCounts(sims)
# dec.p2 <- modelGeneVar(sims)
# # feature selection
# p2.chosen <- getTopHVGs(dec.p2, n=2000)
# sims <- sims[p2.chosen,]
# 
# swap <- grepl("oppo", index)
# 
# max_indices <- apply(res.stm$theta, 1, which.max)
# colnames(res.stm$theta) <- paste0("topic_", 1:ncol(res.stm$theta))
# rownames(res.stm$theta) <- colnames(res.stm$mu$mu)
# res_cluster <- colnames(res.stm$theta)[max_indices]
# names(res_cluster) <- rownames(res.stm$theta)
# res$scSTM_cluster <- res_cluster[match(names(res_cluster), res$Cell)]
# 
# res.adj <- data.frame()
# res.adj <- data.frame(
#   scSTM_adjR = adjustedRandIndex(res$scSTM_cluster,sims$Group),
#   Seurat_adjR = adjustedRandIndex(res$seurat_cluster, sims$Group), 
#   raceID_adjR = adjustedRandIndex(res$raceID_cluster,sims$Group)
#   # CIDR_adjR = adjustedRandIndex(dat$cidr_cluster,sims$Group),
#   # scSTM_sil = mean(scSTM.sil[,3]),
#   # seurat_sil = mean(seurat.sil[,3]),
#   # raceID_sil = mean(raceid.sil[,3]),
#   # CIDR_sil = mean(cidr.sil[,3])
# )
