# this script is designed to assess if the pvalue is uniformly distributed.
# we are going to simulate only the null distribution
# with decreased cells and genes
# simulate 500 samples each
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

##################################### cellType5 #####################################
# 10 batches, 2000 genes, 200 cells each batch, 5 cell tupes
# each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 6
nCellType <- 5
batchCells <- rep(5000, nsample)
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.01
out.prob <- 0.05 # outlier expr prob
out.facScale = 0.5 # count of outlier
out.facLoc = 4 # how far from the main
dropout.type = "experiment" #uses the same parameters for every cell in the same batch
dropout.shape = -1

comp <- sample(1:200,nCellType)
sd <- 2

# simulate timepoint prob
cond_A <- do.call(cbind, lapply(comp, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp, function(mean) rnorm(nsample, mean=mean, sd = sd)))
group.prob <- colSums(cond_A + cond_B)
group.prob <- group.prob/sum(group.prob)

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = 1000, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

sims$time <- NA
cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)

# the end product has 800 cells per batch
N = 600/2
# Apply the function for each combination of Batch and Group
for (batch in unique(sims$Batch)) {

  cond_A_value <- cond_A[match(batch, unique(sims$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims$Batch)), ]

  idx <- which(sims$Batch == batch)

  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)

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

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V6/sims/",
                                "sims_", seed, "_neg_L1_c5.rds"))

##################################### cellType9 #####################################
# 6 batches, 1000 genes, 500 cells each batch, 9 cell tupes
# each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 6
nCellType <- 9
batchCells <- rep(5000, nsample)
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.01
out.prob <- 0.05 # outlier expr prob
out.facScale = 0.5 # count of outlier
out.facLoc = 4 # how far from the main
dropout.type = "experiment" #uses the same parameters for every cell in the same batch
dropout.shape = -1

comp <- sample(1:200,nCellType)
sd <- 2

# simulate timepoint prob
cond_A <- do.call(cbind, lapply(comp, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp, function(mean) rnorm(nsample, mean=mean, sd = sd)))
group.prob <- colSums(cond_A + cond_B)
group.prob <- group.prob/sum(group.prob)

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = 1000, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

sims$time <- NA
cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)

# the end product has 800 cells per batch
N = 600/2
# Apply the function for each combination of Batch and Group
for (batch in unique(sims$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims$Batch)), ]
  
  idx <- which(sims$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
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

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V6/sims/", 
                            "sims_", seed, "_neg_L1_c9.rds"))