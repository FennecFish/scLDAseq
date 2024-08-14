# this script is designed to create simulations under alternative hypothesis
# to evaluate composition changes
# with true comp follows a normal distribution, across each patient
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
# 6 batches, 2000 genes, 2000 cells each batch/timepoint, 5 cell tupes
# each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 6
nCellType <- 5
nGenes <- 3000
batchCells <- rep(10000, nsample)
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.1

# simulate timepoint prob
group.prob <- rep(1/nCellType, nCellType)
params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)
sims$time <- NA

# the end product has 600 cells per batch
N = 4000/2

# leve1, large change, small sd
sims_pos <- sims
comp_A <- c(0.1, 0.1, 0.3 ,0.3, 0.2) * N
comp_B <- c(0.3, 0.3, 0.1, 0.1, 0.2) * N
sd <- 5

cond_A <- do.call(cbind, lapply(comp_A, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp_B, function(mean) rnorm(nsample, mean=mean, sd = sd)))

cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)


# Apply the function for each combination of Batch and Group
for (batch in unique(sims_pos$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims_pos$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims_pos$Batch)), ]
  
  idx <- which(sims_pos$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA) 
    
    sims_pos$time[match(idx$Cell, sims_pos$Cell)] <- 1
    
    remaining_pool <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB) 
    
    sims_pos$time[match(remaining_pool$Cell, sims_pos$Cell)] <- 2
  }
}

sims_pos <- sims_pos[,!is.na(sims_pos$time)]

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/pos_V1/sims/", 
                                "sims_", seed, "_pos_L1_c5.rds"))
cat("nCellType 5 with Level 1 \n")
##### L2 #######
# leve1, large change, small sd
sims_pos <- sims
comp_A <- c(0.15, 0.15, 0.25 ,0.25, 0.2) * N
comp_B <- c(0.2, 0.2, 0.2, 0.2, 0.2) * N
sd <- 10

cond_A <- do.call(cbind, lapply(comp_A, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp_B, function(mean) rnorm(nsample, mean=mean, sd = sd)))

cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)


# Apply the function for each combination of Batch and Group
for (batch in unique(sims_pos$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims_pos$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims_pos$Batch)), ]
  
  idx <- which(sims_pos$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA) 
    
    sims_pos$time[match(idx$Cell, sims_pos$Cell)] <- 1
    
    remaining_pool <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB) 
    
    sims_pos$time[match(remaining_pool$Cell, sims_pos$Cell)] <- 2
  }
}

sims_pos <- sims_pos[,!is.na(sims_pos$time)]

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/pos_V1/sims/", 
                                "sims_", seed, "_pos_L2_c5.rds"))
cat("nCellType 5 with Level 2 \n")
##### L3 #######
# leve1, small change
sims_pos <- sims
comp_A <- c(0.18, 0.18, 0.21 ,0.21, 0.2) * N
comp_B <- c(0.22, 0.22, 0.19, 0.19, 0.2) * N
sd <- 5

cond_A <- do.call(cbind, lapply(comp_A, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp_B, function(mean) rnorm(nsample, mean=mean, sd = sd)))

cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)


# Apply the function for each combination of Batch and Group
for (batch in unique(sims_pos$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims_pos$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims_pos$Batch)), ]
  
  idx <- which(sims_pos$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA) 
    
    sims_pos$time[match(idx$Cell, sims_pos$Cell)] <- 1
    
    remaining_pool <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB) 
    
    sims_pos$time[match(remaining_pool$Cell, sims_pos$Cell)] <- 2
  }
}

sims_pos <- sims_pos[,!is.na(sims_pos$time)]

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/pos_V1/sims/", 
                                "sims_", seed, "_pos_L3_c5.rds"))
cat("nCellType 5 with Level 3 \n")
##### L4 #######
# leve1, minimal change
sims_pos <- sims
comp_A <- c(0.19, 0.21, 0.2 ,0.2, 0.2) * N
comp_B <- c(0.21, 0.19, 0.2, 0.2, 0.2) * N
sd <- 2

cond_A <- do.call(cbind, lapply(comp_A, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp_B, function(mean) rnorm(nsample, mean=mean, sd = sd)))

cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)


# Apply the function for each combination of Batch and Group
for (batch in unique(sims_pos$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims_pos$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims_pos$Batch)), ]
  
  idx <- which(sims_pos$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA) 
    
    sims_pos$time[match(idx$Cell, sims_pos$Cell)] <- 1
    
    remaining_pool <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB) 
    
    sims_pos$time[match(remaining_pool$Cell, sims_pos$Cell)] <- 2
  }
}

sims_pos <- sims_pos[,!is.na(sims_pos$time)]

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/pos_V1/sims/", 
                                "sims_", seed, "_pos_L4_c5.rds"))
cat("nCellType 5 with Level 4 \n")
##################################### cellType9 #####################################
# 6 batches, 2000 genes, 200 cells each batch, 9 cell tupes
# each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 6
nCellType <- 9
batchCells <- rep(10000, nsample)
nGenes <- 3000
de.prob <- runif(nCellType, 0.05, 0.4) # group de gene prob
de.facLoc <- 0.1

# simulate timepoint prob
group.prob <- rep(1/nCellType, nCellType)
params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)
sims$time <- NA

# the end product has 600 cells per batch
N = 4000/2

#### leve1 #####
# large change, small sd
sims_pos <- sims
comp_A <- c(0.05, 0.1, 0.03 ,0.02, 0.3, 0.25, 0.05, 0.1, 0.1) * N
comp_B <- c(0.2, 0.02, 0.15 ,0.06, 0.1, 0.2, 0.07, 0.1, 0.1) * N
sd <- 1

cond_A <- do.call(cbind, lapply(comp_A, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp_B, function(mean) rnorm(nsample, mean=mean, sd = sd)))

cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)


# Apply the function for each combination of Batch and Group
for (batch in unique(sims_pos$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims_pos$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims_pos$Batch)), ]
  
  idx <- which(sims_pos$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA) 
    
    sims_pos$time[match(idx$Cell, sims_pos$Cell)] <- 1
    
    remaining_pool <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB) 
    
    sims_pos$time[match(remaining_pool$Cell, sims_pos$Cell)] <- 2
  }
}

sims_pos <- sims_pos[,!is.na(sims_pos$time)]

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/pos_V1/sims/", 
                                "sims_", seed, "_pos_L1_c9.rds"))
cat("nCellType 9 with Level 1 \n")
##### L2 #######
# leve1, large change, small sd
sims_pos <- sims
comp_A <- c(0.15, 0.1, 0.05 ,0.13, 0.17, 0.15, 0.05, 0.1, 0.1) * N
comp_B <- c(0.1, 0.15, 0.05 ,0.14, 0.16, 0.13, 0.07, 0.1, 0.1) * N
sd <- 2

cond_A <- do.call(cbind, lapply(comp_A, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp_B, function(mean) rnorm(nsample, mean=mean, sd = sd)))

cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)


# Apply the function for each combination of Batch and Group
for (batch in unique(sims_pos$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims_pos$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims_pos$Batch)), ]
  
  idx <- which(sims_pos$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA) 
    
    sims_pos$time[match(idx$Cell, sims_pos$Cell)] <- 1
    
    remaining_pool <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB) 
    
    sims_pos$time[match(remaining_pool$Cell, sims_pos$Cell)] <- 2
  }
}

sims_pos <- sims_pos[,!is.na(sims_pos$time)]

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/pos_V1/sims/", 
                                "sims_", seed, "_pos_L2_c9.rds"))
cat("nCellType 9 with Level 2 \n")
##### L3 #######
# leve1, small change
sims_pos <- sims
comp_A <- c(0.09, 0.11, 0.1549 ,0.0951, 0.075, 0.1075, 0.1675, 0.1, 0.1) * N
comp_B <- c(0.11, 0.09, 0.1551 ,0.0949, 0.070, 0.1080, 0.172, 0.1, 0.1) * N
sd <- 1

cond_A <- do.call(cbind, lapply(comp_A, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp_B, function(mean) rnorm(nsample, mean=mean, sd = sd)))

cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)


# Apply the function for each combination of Batch and Group
for (batch in unique(sims_pos$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims_pos$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims_pos$Batch)), ]
  
  idx <- which(sims_pos$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA) 
    
    sims_pos$time[match(idx$Cell, sims_pos$Cell)] <- 1
    
    remaining_pool <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB) 
    
    sims_pos$time[match(remaining_pool$Cell, sims_pos$Cell)] <- 2
  }
}

sims_pos <- sims_pos[,!is.na(sims_pos$time)]

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/pos_V1/sims/", 
                                "sims_", seed, "_pos_L3_c9.rds"))
cat("nCellType 9 with Level 3 \n")
##### L4 #######
# leve1, minimal change
sims_pos <- sims
comp_A <- c(0.09, 0.11, 0.15 ,0.15, 0.1, 0.1, 0.1, 0.1, 0.1) * N
comp_B <- c(0.11, 0.09, 0.15 ,0.15, 0.1, 0.1, 0.1, 0.1, 0.1) * N
sd <- 1

cond_A <- do.call(cbind, lapply(comp_A, function(mean) rnorm(nsample, mean=mean, sd = sd)))
cond_B <- do.call(cbind, lapply(comp_B, function(mean) rnorm(nsample, mean=mean, sd = sd)))

cond_A <- cond_A/rowSums(cond_A)
cond_B <- cond_B/rowSums(cond_B)


# Apply the function for each combination of Batch and Group
for (batch in unique(sims_pos$Batch)) {
  
  cond_A_value <- cond_A[match(batch, unique(sims_pos$Batch)), ]
  cond_B_value <- cond_B[match(batch, unique(sims_pos$Batch)), ]
  
  idx <- which(sims_pos$Batch == batch)
  
  n_A <- round(N * cond_A_value)
  n_B <- round(N * cond_B_value)
  
  for (i in 1:length(n_A)) {
    group_name <- paste0("Group", i)
    nA <- n_A[i]
    nB <- n_B[i]
    
    idx <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      sample_n(nA) 
    
    sims_pos$time[match(idx$Cell, sims_pos$Cell)] <- 1
    
    remaining_pool <- colData(sims_pos) %>%
      as.data.frame() %>%
      filter(Batch == batch) %>%
      filter(Group == group_name) %>%
      filter(!(Cell %in% idx$Cell)) %>%
      sample_n(nB) 
    
    sims_pos$time[match(remaining_pool$Cell, sims_pos$Cell)] <- 2
  }
}

sims_pos <- sims_pos[,!is.na(sims_pos$time)]

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/pos_V1/sims/", 
                                "sims_", seed, "_pos_L4_c9.rds"))
cat("nCellType 9 with Level 4 \n")