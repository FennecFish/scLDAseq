setwd("/proj/milovelab/wu/scLDAseq")
# library(splatter)
library(Matrix)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(scater)
library(scran)
library(sva)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V7/sims/", file_name))
# quick qc
sims <- quickPerCellQC(sims, filter=TRUE)

# batch effect removal using combat seq
adjusted_counts <- ComBat_seq(counts(sims), batch=sims$Batch, group=NULL)
counts(sims) <- adjusted_counts
cat("Batch Effect Removed \n")

### remove genes with count 0
sims <- sims[rowSums(counts(sims)) != 0,]
#### feature selection #####
sims <- scuttle::logNormCounts(sims)

dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1000)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()

scSTM.mod <- selectModel(sce = sims,
                         K = ngroup, prevalence = ~time, content = NULL,
                         sample = "Batch", N = 3, ts_runs = 10, random_run =10)

all_values <- unlist(scSTM.mod$bound)
max_value <- max(all_values)
max_position_in_vector <- which(all_values == max_value)
res <- scSTM.mod$runout[[max_position_in_vector]]
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(res, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V7/scSTM/",
                           "scSTM_", set_level, ".rds"))
