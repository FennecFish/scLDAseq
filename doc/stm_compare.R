setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
# setwd("/proj/milovelab/wu/scLDAseq")
set.seed(1)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(dplyr)
library(tibble)
library(stats)
library(splatter)
library(scater)
library(MASS)

sims <- readRDS("data/sims_1716946229_neg_L1_c5.rds")
##### QA ######
sims <- quickPerCellQC(sims, filter=TRUE)

### remove genes with count 0 
sims <- sims[rowSums(counts(sims)) != 0,]
#### feature selection #####
sims <- scuttle::logNormCounts(sims)
library(scater)
library(scran)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1000)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

args <- prepsce(sims)
documents <- args$documents
vocab <- args$vocab
data <- args$meta

r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
sourceCpp("../stm/src/STMCfuns.cpp")

res <- stm(documents = documents, vocab = vocab, K = ngroup,
           data = data,
                 prevalence = ~time, content = NULL,
                 init.type= "Spectral",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 control = list(gamma.maxits=3000),
                 emtol=1e-5,
                 seed = 9308641, max.em.its = 10)
plot(res$convergence$bound)
plot(res$convergence$bound1)
