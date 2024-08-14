# this script is to run scSTMseq with filtered gene
# without content
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(stm)
library(scater)
library(scran)
library(doParallel)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
gc <- as.integer(args[2])
cat(file_name, "\n")

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/de_gene_eval/sims/", file_name))
# quick qc
sims <- quickPerCellQC(sims, filter=TRUE)

### remove genes with count 0 
sims <- sims[rowSums(counts(sims)) != 0,]
#### feature selection #####
sims <- scuttle::logNormCounts(sims)

dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=800)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()
scSTMobj <- scSTMseq(sims, K = ngroup, sample = "Sample",
                                 prevalence=~Time, content=~Time, 
                     init.type = "Random", seed=1229671)
# scSTM.mod <- stm(sce = sims, K = ngroup, prevalence = ~time, init.type = "Spectral")
# scSTM.mod <- selectModel(sce = sims,
#                          K = ngroup, prevalence = ~Batch, content = ~Batch,
#                          N = 1, ts_runs = 1, random_run = 1,
#                          max.em.its = 100, net.max.em.its = 5) 
scSTM.mod <- selectModel_parallel(sce = sims, sample = "Sample",
                                  K = ngroup, prevalence = ~Time, content = ~Time,
                                  N = 10, ts_runs = 50, random_run = 50,
                                  max.em.its = 100, net.max.em.its = 15, gc = gc)

msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(scSTM.mod, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/de_gene_eval/scSTM_Content_Prevalance/", 
                                 "scSTM_", set_level, ".rds"))

