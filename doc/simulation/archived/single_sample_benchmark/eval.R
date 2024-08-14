setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(Seurat)
library(cluster)
library(monocle3)

process_scSTM <- function(scSTMobj) {
  if(class(scSTMobj) == "selectModel") {
    all_values <- unlist(scSTMobj$bound)
    max_value <- max(all_values)
    max_position_in_vector <- which(all_values == max_value)
    scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
  }
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- colnames(scSTMobj$mu$mu)
  if(is.null(rownames(scSTMobj$theta))){
    rownames(scSTMobj$theta) <- colnames(scSTMobj$settings$sce)
  }
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sims/", pattern = "sims*")
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", pattern = "scSTM*")

res <- data.frame()

for(file_name in files){
  # create a new dataframe to calculate adjRandIndex
  dat <- data.frame()
  
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sims/sims_", set_level, ".rds"))
  dat <- colData(sims) %>% data.frame() 
  
  # scSTM_nC_P
  scSTM_nC_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/scSTM_noContent_Prevalance/scSTM_",set_level,".rds")
  if(file.exists(scSTM_nC_P_name)){
    scSTM_nC_P <- readRDS(scSTM_nC_P_name)
    scSTM_nC_P_cluster <- process_scSTM(scSTM_nC_P)
    dat$scSTM_nC_P_cluster <- scSTM_nC_P_cluster[match(rownames(dat), names(scSTM_nC_P_cluster))]
  } else {dat$scSTM_nC_P_cluster <- NA}
  
  # scSTM_C_nP
  scSTM_C_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/scSTM_Content_noPrevalance/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_P_name)){
    scSTM_C_nP <- readRDS(scSTM_C_P_name)
    scSTM_C_nP_cluster <- process_scSTM(scSTM_C_nP)
    dat$scSTM_C_nP_cluster <- scSTM_nC_P_cluster[match(rownames(dat), names(scSTM_C_nP_cluster))]
  } else {dat$scSTM_C_nP_cluster <- NA}
  
  # scSTM_C_P
  scSTM_C_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/scSTM_Content_Prevalance/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_P_name)){
    scSTM_C_P <- readRDS(scSTM_C_P_name)
    scSTM_C_P_cluster <- process_scSTM(scSTM_C_P)
    dat$scSTM_C_P_cluster <- scSTM_C_P_cluster[match(rownames(dat), names(scSTM_C_P_cluster))]
  } else {dat$scSTM_C_P_cluster <- NA}
  
  # Seurat
  seurat_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/seurat/seurat_", set_level, ".rds")
  if(file.exists(seurat_name)){
    seurat.sims <- readRDS(seurat_name)
    smeta <- seurat.sims@meta.data %>% as.data.frame()
    sub_sims <- sims[,rownames(smeta)] # filter by the rows
    seurat.adj <- sapply(smeta[,4:7], function(x) {
      adjustedRandIndex(x, sub_sims$Group)
    })
    # select the resolution that has the highest ARI. When there are multiple, select the first one
    best_res <- names(seurat.adj)[seurat.adj == max(seurat.adj)][1]
    seurat_cluster <- seurat.sims@meta.data %>% as.data.frame() %>% dplyr::select(all_of(best_res))
    dat$seurat_cluster <- seurat_cluster[match(rownames(dat), rownames(seurat_cluster)),]
    rm(sub_sims)
  } else{dat$seurat_cluster <- NA}

  
  # sctransform
  sctf_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sctransform/sctransform_", set_level, ".rds")
  if(file.exists(sctf_name)){
    sctf <- readRDS(sctf_name)
    sct_meta <- sctf@meta.data %>% as.data.frame()
    sub_sims <- sims[,rownames(sct_meta)] # filter by the rows
    sctf.adj <- sapply(sct_meta[,6:9], function(x) {
      adjustedRandIndex(x, sub_sims$Group)
    })
    # select the resolution that has the highest ARI. When there are multiple, select the first one
    best_res <- names(sctf.adj)[sctf.adj == max(sctf.adj)][1]
    sctf_cluster <- sctf@meta.data %>% as.data.frame() %>% dplyr::select(all_of(best_res))
    dat$sctransform_cluster <- sctf_cluster[match(rownames(dat), rownames(sctf_cluster)),]
    rm(sub_sims)
  }else{dat$sctransform_cluster <- NA}
  
  # fastTopics
  fasttopic_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/fastTopics/fastTopics_", set_level, ".rds")
  if(file.exists(fasttopic_name)){
    nmf.sims <- readRDS(fasttopic_name)
    max_indices <- apply(nmf.sims$L, 1, which.max)
    fastTopics_cluster <- colnames(nmf.sims$L)[max_indices]
    names(fastTopics_cluster) <- rownames(nmf.sims$L)
    dat$fastTopics_cluster <- fastTopics_cluster[match(rownames(dat), names(fastTopics_cluster))]
  } else{ dat$fastTopics_cluster <- NA }

  # monocle3
  monocle_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/monocle3/monocle3_", set_level, ".rds")
  if(file.exists(monocle_name)){
    monocle3 <- readRDS(monocle_name)
    dat$monocle3_cluster <- partitions(monocle3)[match(rownames(dat), names(partitions(monocle3)))]
  } else {dat$monocle3_cluster <- NA}
  
  adjusted_rand_indices <- sapply(dat[, 7:ncol(dat)], function(x) {
    adjustedRandIndex(x, sims$Group)
  })

  # Create a data frame to store results
  res.temp <- data.frame(
    sim = set_level,
    seed = unlist(strsplit(set_level, "_"))[1],
    control = unlist(strsplit(set_level, "_"))[2],
    level = unlist(strsplit(set_level, "_"))[3],
    nCellType = unlist(strsplit(set_level, "_"))[4],
    nSample = unlist(strsplit(set_level, "_"))[5],
    t(adjusted_rand_indices)  # Transpose to match the original data frame structure
  )
  
  res <- bind_rows(res, res.temp)
  cat(file_name, "\n")
  rm(sims)
}

write.csv(res, file = "res/adjRandIndex_single_sample_benchmark.csv")
