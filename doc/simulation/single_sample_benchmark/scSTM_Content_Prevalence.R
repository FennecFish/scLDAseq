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

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/sims/", file_name))
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
# scSTM.mod <- stm(sce = sims, K = ngroup, prevalence = ~time, init.type = "Spectral")
# scSTM.mod <- selectModel(sce = sims,
#                          K = ngroup, prevalence = ~Batch, content = ~Batch,
#                          N = 1, ts_runs = 1, random_run = 1,
#                          max.em.its = 100, net.max.em.its = 5) 
scSTM.mod <- selectModel(sce = sims,
                         K = ngroup, prevalence = ~Batch, content = ~Batch,
                         N = 6, ts_runs = 30, random_run = 30,
                         max.em.its = 100, net.max.em.its = 10) #, sample = "Batch",)

# all_values <- unlist(scSTM.mod$bound)
# max_value <- max(all_values)
# max_position_in_vector <- which(all_values == max_value)
# res <- scSTM.mod$runout[[max_position_in_vector]]
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(scSTM.mod, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/scSTM_Content_Prevalance/", 
                                 "scSTM_", set_level, ".rds"))

