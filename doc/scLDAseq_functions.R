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

scLDAseq <- function(sims, verbose = TRUE) {
  
  ngroup <- length(unique(sims$Group))
  dat <- colData(sims) %>% 
    data.frame() 
  
  ###########################################################
  ################## scLDAseq ###############################
  ###########################################################
  r.file <- paste0("R/",list.files("R/"))
  sapply(r.file, source)
  sourceCpp("src/STMCfuns.cpp")
  
  t1 <- proc.time()
  
  K <- ngroup
  
  stm_dat <- prepsce(sims)
  prevalence <- as.formula(~stm_dat$meta$time)
  content <- NULL
  sce <- stm_dat$sce
  documents  <- stm_dat$documents
  vocab <- stm_dat$vocab
  data <- stm_dat$meta
  sample <- "Batch"
  
  res.stm <- multi_stm(documents = documents, vocab = vocab,
                       K = K, prevalence = prevalence, content = NULL,
                       data = data,
                       sce = sce,
                       sample = sample,
                       init.type= "Spectral",
                       gamma.prior= "Pooled",
                       kappa.prior= "L1",
                       control = list(gamma.maxits=3000))
  return(res.stm)
  # max_indices <- apply(res.stm$theta, 1, which.max)
  # colnames(res.stm$theta) <- paste0("topic_", 1:ncol(res.stm$theta))
  # rownames(res.stm$theta) <- colnames(res.stm$mu$mu)
  # res_cluster <- colnames(res.stm$theta)[max_indices]
  # names(res_cluster) <- rownames(res.stm$theta)
  # dat$scSTM_cluster <- res_cluster[match(names(res_cluster), dat$Cell)]
  
  msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
  if(verbose) cat(msg)
  
}

