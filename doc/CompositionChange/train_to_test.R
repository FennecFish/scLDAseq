setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(fastTopics)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

# file_name <- "sims_1716586115_neg_L1_c5.rds"

set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file_name)

scSTM <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/scSTM/scSTM_", set_level,".rds"))
test_dat <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/test/",
                     "test_", set_level, ".rds"))

beta <- lapply(scSTM$beta$logbeta, exp)[[1]]
X <- test_dat
X <- X[rownames(X) %in% scSTM$vocab,colnames(X) %in% scSTM$DocName] # subset test data in case any filtering
tX <- as.matrix(t(X))

fit0 <- init_poisson_nmf(tX, F = t(beta), init.method = "random")
res <- fit_poisson_nmf(tX, fit0 = fit0, numiter = 100, update.factors = NULL, method = "em",control = list(nc=2))
fit.multinom <- poisson2multinom(res)
saveRDS(fit.multinom, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/test_fastTopics/",
                           "fastTopics_", set_level, ".rds"))