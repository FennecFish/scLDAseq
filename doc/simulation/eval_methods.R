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
  # dist.matrix <- dist(t(counts(sims)))
  # scSTM.sil <- silhouette(as.numeric(as.factor(dat$scSTM_cluster)), dist.matrix)
  # seurat.sil <- silhouette(as.numeric(as.factor(dat$seurat_cluster)), dist.matrix)
  # raceid.sil <- silhouette(as.numeric(as.factor(dat$raceID_cluster)), dist.matrix)
  # cidr.sil <- silhouette(as.numeric(as.factor(dat$cidr_cluster)), dist.matrix)
  
  res <- data.frame(
    scSTM_adjR = adjustedRandIndex(dat$scSTM_cluster,sims$Group),
    Seurat_adjR = adjustedRandIndex(dat$seurat_cluster, sims$Group), 
    raceID_adjR = adjustedRandIndex(dat$raceID_cluster,sims$Group),#,
    CIDR_adjR = adjustedRandIndex(dat$cidr_cluster,sims$Group))#,
        # scSTM_sil = mean(scSTM.sil[,3]),
        # seurat_sil = mean(seurat.sil[,3]),
        # raceID_sil = mean(raceid.sil[,3]))

  return(res)
}

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/", pattern = "3cellTypes_sims.rds")
sim_name <- sub("\\_sims.rds$", "", files)
res <- data.frame()

for (i in sim_name){
  dat <- read.csv(paste0("res/clustering_benchmark/colData_",i,".csv"))
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",i,"_sims.rds"))
  temp <- sc_eval(sims, dat = dat)
  res <- rbind(res, temp)
  rm(dat)
  rm(sims)
  cat(i, "\n")
}

rownames(res) <- sim_name
res <- tibble::rownames_to_column(res, "sim")
res$NumCellType <- sub(".*samples_(\\d+)cellTypes", "\\1", res$sim)
res$NumSample <- sub(".*_(\\d+samples).*", "\\1", res$sim)
write.csv(res, file = "res/res_eval_clustering_4methods.csv")

# change the format from wide to long
library(tidyr)
res <- read.csv("res/res_eval_methods.csv")
wcidr <- read.csv("res/res_eval_clustering_4methods.csv")

wcidr_long <- wcidr %>% gather(method, adjR, scSTM_adjR:CIDR_adjR, factor_key=TRUE)
# gather(res, method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
ggplot(wcidr_long, aes(x=method, y=adjR)) + 
  geom_boxplot() + 
  ggtitle("Comparison of Clustering Accuracy among Methods with 3 cell Types")

res <- res %>%
  mutate(num_samples = as.numeric(gsub(".*_(\\d+)samples.*", "\\1", sim)),
         num_cell_types = as.numeric(gsub(".*_(\\d+)cellTypes", "\\1", sim)))
data_long <- res %>% filter(num_cell_types == 5) %>% gather(method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
# gather(res, method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
ggplot(data_long, aes(x=method, y=adjR)) + 
  geom_boxplot() + 
  ggtitle("Comparison of Clustering Accuracy among Methods with 5 cell Types")


# looking at why scSTM is low

sim_name <- "BIOKEY_11_10samples_5cellTypes"
sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",sim_name,"_sims.rds"))
