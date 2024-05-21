setwd("/proj/milovelab/wu/scLDAseq")
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(scater)
library(scran)
library(sva)

dat <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_sce/anti_PD1_cohort1_sce.rds")
# dat$batch <- paste0(dat$patient_id, "_", dat$timepoint)
# quick qc
dat <- quickPerCellQC(dat, filter=TRUE)

# batch effect removal using combat seq
adjusted_counts <- ComBat_seq(counts(dat), batch=dat$patient_id, group=dat$timepoint)
counts(dat) <- adjusted_counts
saveRDS(dat, file = "/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_sce/PD1_cohort1_batchCorrected.rds")
cat("Batch Effect Removed \n")

### remove genes with count 0 
dat <- dat[rowSums(counts(dat)) != 0,]
#### feature selection #####
dat <- scuttle::logNormCounts(dat)

dec.p2 <- modelGeneVar(dat)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=2500)
dat <- dat[p2.chosen,]

nsample <- length(unique(dat$patient_id))
ngroup <- length(unique(dat$cellType))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()

scSTM.mod <- selectModel(sce = dat,
                         K = ngroup, prevalence = ~timepoint, content = NULL,
                         sample = "patient_id", N = 5, ts_runs = 20, random_run =20)

all_values <- unlist(scSTM.mod$bound)
max_value <- max(all_values)
max_position_in_vector <- which(all_values == max_value)
res <- scSTM.mod$runout[[max_position_in_vector]]
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(res, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/scSTM_combat_f_nc/", 
                           "scSTM_", set_level, ".rds"))