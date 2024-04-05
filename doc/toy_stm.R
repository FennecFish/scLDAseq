# setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
setwd("/proj/milovelab/wu/scLDAseq")
set.seed(1)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(mclust)
library(lme4)
library(lmerTest)
library(dplyr)
library(tibble)
library(stats)
library(scuttle)
library(scran)

sims <- readRDS("data/toydat.rds")

#### QC ######
sims <- quickPerCellQC(sims)
#### feature selection #####
sims <- scuttle::logNormCounts(sims)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1000)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

#### scLDAseq#############
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims)
K <- length(unique(sims$Group))

# STM for reference
library(stm)
r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
sourceCpp("../stm/src/STMCfuns.cpp")
res.stm <- stm(documents = dat$documents, vocab = dat$vocab,
               K = K, prevalence = ~time, content = NULL,
               data = dat$meta, 
               init.type= "Spectral",
               gamma.prior= "Pooled",
               kappa.prior= "L1",
               control = list(gamma.maxits=3000))
saveRDS(res.scLDAseq, file = "data/stm_toydat.rds")
