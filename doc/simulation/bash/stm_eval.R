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
library(cluster)

set.seed(1)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

file <- list.files(path = "data/", pattern = "_sims.rds")

swap = TRUE

for (index in file) {
  
  samp <- sub("_sims.rds", "", index)
  sims <- readRDS(paste0("data/", index))
  
  if (swap) {
    sampled_data <- colData(sims) %>%
      data.frame() %>%
      mutate(new_time = ifelse(time == 1 & Batch %in% c("Batch1", "Batch2"), 2, 1)) %>%
      mutate(time = ifelse(Batch %in% c("Batch1", "Batch2"), new_time, time))
    sims$time <- sampled_data$time
    samp <- paste0(samp,"_oppo")
  }
  
  #### QC ######
  sims <- quickPerCellQC(sims)
  #### feature selection #####
  sims <- scuttle::logNormCounts(sims)
  dec.p2 <- modelGeneVar(sims)
  # feature selection
  p2.chosen <- getTopHVGs(dec.p2, n=2000)
  sims<- sims[p2.chosen,]
  
  ngroup <- length(unique(sims$Group))
  K <- ngroup
  
  stm_dat <- prepsce(sims)
  cat("prepsc completed for", samp, "\n")
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
  cat("multi_stm completed for", samp, "\n")
  
 saveRDS(res.stm, file = paste0("res/", samp, "_stmRes.rds"))
}





library(tidyr)
adjR_long <- gather(adj_dat, method, adjustedRandIndex, scSTM_adjR:CIDR_adjR, factor_key=TRUE)

ggplot(adjR_long, aes(x=method, y=adjustedRandIndex)) + 
  geom_boxplot()

sil_long <- gather(adj_dat, method, silhoutte_score, scSTM_sil:CIDR_sil, factor_key=TRUE)
ggplot(sil_long, aes(x=method, y=silhoutte_score)) + 
  geom_boxplot()
