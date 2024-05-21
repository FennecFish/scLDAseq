# this script is designed to assess if the pvalue is uniformly distributed.
# we are going to simulate only the null distribution
# with decreased cells and genes
# simulate 1000 samples each
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(tidyverse)
library(splatter)
library(SingleCellExperiment)
library(tidyr)


args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 0,1000))

##################################### cellType3 #####################################
# 10 batches, 2000 genes, 200 cells each batch, 3 cell tupes
# each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 10
nCellType <- 3
batchCells <- batchCells <- rep(200, nsample)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.01
out.prob <- 0.05 # outlier expr prob
out.facScale = 0.5 # count of outlier
out.facLoc = 4 # how far from the main
dropout.type = "experiment" #uses the same parameters for every cell in the same batch
dropout.shape = -1

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = 2000, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# first simulate negative control: no change
sims_neg <- sims
sc_data_neg <- colData(sims_neg) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))))
sims_neg$time <- sc_data_neg$time

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", 
                                "sims_", seed, "_neg_L1_c3.rds"))
cat("Generation Completed for c3 neg.\n")


##################################### cellType 5 ###################################
nsample <- 10
nCellType <- 5
batchCells <- batchCells <- rep(200, nsample)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.01
out.prob <- 0.05 # outlier expr prob
out.facScale = 0.5 # count of outlier
out.facLoc = 4 # how far from the main
dropout.type = "experiment" #uses the same parameters for every cell in the same batch
dropout.shape = -1

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = 2000, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# first simulate negative control: no change
sims_neg <- sims
sc_data_neg <- colData(sims_neg) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))))
sims_neg$time <- sc_data_neg$time

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", 
                                "sims_", seed, "_neg_L1_c5.rds"))
cat("Generation Completed for c5 neg.\n")

##################################### cellType 9 ##################################
nsample <- 10
nCellType <- 9
batchCells <- batchCells <- rep(200, nsample)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.01
out.prob <- 0.05 # outlier expr prob
out.facScale = 0.5 # count of outlier
out.facLoc = 4 # how far from the main
dropout.type = "experiment" #uses the same parameters for every cell in the same batch
dropout.shape = -1

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = 2000, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# first simulate negative control: no change
sims_neg <- sims
sc_data_neg <- colData(sims_neg) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))))
sims_neg$time <- sc_data_neg$time

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", 
                                "sims_", seed, "_neg_L1_c9.rds"))
cat("Generation Completed for c9 neg.\n")
