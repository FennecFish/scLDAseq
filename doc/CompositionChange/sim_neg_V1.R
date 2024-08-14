# this script is designed to assess if the pvalue is uniformly distributed.
# and test Type I error
# we are going to simulate only the null distribution
# using multivariate normal distribution
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(tidyverse)
library(splatter)
library(SingleCellExperiment)
library(tidyr)
library(MASS)

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 0,1000))

##################################### cellType5 #####################################
# 10 batches, 2000 genes, 200 cells each batch, 5 cell tupes
# each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 5
nCellType <- 5
nGenes <- 1000
batchCells <- rep(10000, nsample)
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.1

group.prob <- rep(1/nCellType, nCellType)

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    # out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    # dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

sims$time <- NA

# the end product has 3000 cells per batch
N = 4000/2

comp <- c(0.1, 0.2, 0.25 ,0.25, 0.2) * N

cond_A <- round(mvrnorm(n = nsample, mu = comp, Sigma = diag(x=5, nrow =5)))
cond_B <- round(mvrnorm(n = nsample, mu = comp, Sigma = diag(x=5, nrow =5)))

# Apply the function for each combination of Batch and Group
for (batch in unique(sims$Batch)) {
  
  n_A <- cond_A[match(batch, unique(sims$Batch)), ]
  n_B <- cond_B[match(batch, unique(sims$Batch)), ]
  
  idx <- which(sims$Batch == batch)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA)
    
    sims$time[match(idx$Cell, sims$Cell)] <- 1
    
    remaining_pool <- colData(sims) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB)
    
    sims$time[match(remaining_pool$Cell, sims$Cell)] <- 2
  }
}

# sims$time <- ifelse(is.na(sims$time), 2,sims$time)
sims <- sims[,!is.na(sims$time)]
# colData(sims) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V1/sims/",
                                "sims_", seed, "_neg_L1_c5.rds"))
cat("Created for Cell Type 5 \n")
##################################### cellType9 #####################################
# 5 batches, 1000 genes, 500 cells each batch, 9 cell tupes
# each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 5
nCellType <- 9
nGenes <- 1000
batchCells <- rep(10000, nsample)
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.1

group.prob <- rep(1/nCellType, nCellType)

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    # out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    # dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

sims$time <- NA

# the end product has 2000 cells per batch/timepoint
N = 4000/2

comp <- c(0.09, 0.11, 0.13 ,0.1, 0.17, 0.12, 0.1, 0.08, 0.1)  * N

cond_A <- round(mvrnorm(n = nsample, mu = comp, Sigma = diag(x=5, nrow =nCellType)))
cond_B <- round(mvrnorm(n = nsample, mu = comp, Sigma = diag(x=5, nrow =nCellType)))

# Apply the function for each combination of Batch and Group
for (batch in unique(sims$Batch)) {
  
  n_A <- cond_A[match(batch, unique(sims$Batch)), ]
  n_B <- cond_B[match(batch, unique(sims$Batch)), ]
  
  idx <- which(sims$Batch == batch)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA)
    
    sims$time[match(idx$Cell, sims$Cell)] <- 1
    
    remaining_pool <- colData(sims) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB)
    
    sims$time[match(remaining_pool$Cell, sims$Cell)] <- 2
  }
}

# sims$time <- ifelse(is.na(sims$time), 2,sims$time)
sims <- sims[,!is.na(sims$time)]
# colData(sims) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V1/sims/", 
                            "sims_", seed, "_neg_L1_c9.rds"))

cat("Created for Cell Type 9 \n")
