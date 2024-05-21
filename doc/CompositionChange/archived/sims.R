# this script is to create simulation with true known cell type proportion
setwd("/proj/milovelab/wu/scLDAseq")
library(dplyr)
library(scuttle)
library(tidyverse)
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))

##################################### level 2 #####################################
# One batch, 3000 genes, 5000 cells, 5 cell
# types, each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 1
nCellType <- 3
batchCells <- 2000
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.5
out.prob <- 0.05 # outlier expr prob
out.facScale = 0.5 # count of outlier
out.facLoc = 4 # how far from the main
dropout.type = "batch" #uses the same parameters for every cell in the same batch
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
  group_by(Group) %>%
  mutate(time = sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))))
sims_neg$time <- sc_data_neg$time

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                            "sims_neg_", seed, "_L2.rds"))

# simulate negative control
sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group3" ~ sample(c(1,2),n(), replace = TRUE)
  ))
sims_pos$time <- sc_data_pos$time

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                            "sims_pos_", seed, "_L2.rds"))
cat("Generation Completed for L2.\n")
# prop.table(table(sc_data$Group, sc_data$time), margin = 2)

##################################### level 3 ###################################
# One batch, 3000 genes, 5000 cells, 
# 5 cell types, each with differential probability of 0.5, and de.facLoc = 0.5
nsample <- 1
nCellType <- 5
batchCells <- 2000
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.5
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
  group_by(Group) %>%
  mutate(time = sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))))
sims_neg$time <- sc_data_neg$time
# prop.table(table(sc_data_neg$Group, sc_data_neg$time), margin = 2)

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                                "sims_neg_", seed, "_L3.rds"))

# simulate negative control
sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group3" ~ sample(c(rep(1, 2* floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3* floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group5" ~ sample(c(1,2),n(), replace = TRUE)
  ))
#prop.table(table(sc_data_pos$Group, sc_data_pos$time), margin = 2)
sims_pos$time <- sc_data_pos$time

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                                "sims_pos_", seed, "_L3.rds"))

cat("Generation Completed for L3.\n")


##################################### level 4 ###################################
#  Six batches with batch effect. 3000 genes, 3000 cells each batch, 5 cell types, 
# each with differential probability of 0.5

nsample <- 6 
nCellType <- 5
batchCells <- rep(600, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
# rn_group <- runif(nCellType)
# group.prob <- rn_group / sum(rn_group) # cell type proportion
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.5
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
                    nGenes = 3000, 
                    batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale,
                    seed = seed)

sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# let batch1,2,3 be negative control, and 4,5,6 be positive control
sc_data <- colData(sims) %>%
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Batch %in% c("Batch1", "Batch2", "Batch3") ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))),
    Batch %in% c("Batch4", "Batch5", "Batch6") ~ case_when(
      Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group5" ~ sample(c(1, 2), n(), replace = TRUE)
    )
  ))
sims$time <- sc_data$time
# sc_data %>%
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)  # Adjust this line based on actual column names from time categories
saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                                "sims_", seed, "_L4.rds"))

cat("Generation Completed for L4.\n")

##################################### level 5 ###################################
# On top of experiment 4, decrease de.facLoc
nsample <- 6 
nCellType <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.1
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
                    nGenes = 3000, 
                    batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale,
                    seed = seed)

sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# let batch1,2,3 be negative control, and 4,5,6 be positive control
sc_data <- colData(sims) %>%
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Batch %in% c("Batch1", "Batch2", "Batch3") ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))),
    Batch %in% c("Batch4", "Batch5", "Batch6") ~ case_when(
      Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group5" ~ sample(c(1, 2), n(), replace = TRUE)
    )
  ))
sims$time <- sc_data$time
# sc_data %>%
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)  # Adjust this line based on actual column names from time categories
saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                            "sims_", seed, "_L5.rds"))
cat("Generation Completed for L5.\n")

##################################### level 6 #################################
nsample <- 6 
nCellType <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0.1, 0.3) # group de gene prob
de.facLoc <- 0.1
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
                    nGenes = 3000, 
                    batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale,
                    seed = seed)

sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# let batch1,2,3 be negative control, and 4,5,6 be positive control
sc_data <- colData(sims) %>%
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Batch %in% c("Batch1", "Batch2", "Batch3") ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))),
    Batch %in% c("Batch4", "Batch5", "Batch6") ~ case_when(
      Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group5" ~ sample(c(1, 2), n(), replace = TRUE)
    )
  ))
sims$time <- sc_data$time
# sc_data %>%
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)  # Adjust this line based on actual column names from time categories
saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                            "sims_", seed, "_L6.rds"))

cat("Generation Completed for L6.\n")

##################################### level 7 #################################
# further decrease de.facLoc
nsample <- 6 
nCellType <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0.1, 0.3) # group de gene prob
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
                    nGenes = 3000, 
                    batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale,
                    seed = seed)

sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# let batch1,2,3 be negative control, and 4,5,6 be positive control
sc_data <- colData(sims) %>%
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Batch %in% c("Batch1", "Batch2", "Batch3") ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))),
    Batch %in% c("Batch4", "Batch5", "Batch6") ~ case_when(
      Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group5" ~ sample(c(1, 2), n(), replace = TRUE)
    )
  ))
sims$time <- sc_data$time
# sc_data %>%
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)  # Adjust this line based on actual column names from time categories
saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                            "sims_", seed, "_L7.rds"))
cat("Generation Completed for L7.\n")

##################################### level 8 #################################
nsample <- 6 
nCellType <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0, 0.2) # group de gene prob
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
                    nGenes = 3000, 
                    batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale,
                    seed = seed)

sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)
# let batch1,2,3 be negative control, and 4,5,6 be positive control
sc_data <- colData(sims) %>%
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Batch %in% c("Batch1", "Batch2", "Batch3") ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))),
    Batch %in% c("Batch4", "Batch5", "Batch6") ~ case_when(
      Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group5" ~ sample(c(1, 2), n(), replace = TRUE)
    )
  ))
sims$time <- sc_data$time
# sc_data %>%
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)  # Adjust this line based on actual column names from time categories
saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                            "sims_", seed, "_8.rds"))


cat("Generation Completed for L8.\n")

##################################### level 9 ##################################
nsample <- 10
nCellType <- 8
batchCells <- rep(500, nsample)
batch.facLoc <- runif(nsample,0,0.3)
batch.facScale <- runif(nsample,0,0.3)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0, 0.2) # group de gene prob
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
                    nGenes = 3000, 
                    batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale,
                    seed = seed)

sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# let batch1,2,3 be negative control, and 4,5,6 be positive control
sc_data <- colData(sims) %>%
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Batch %in% c("Batch1", "Batch2", "Batch3", "Batch4", "Batch5") ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))),
    Batch %in% c("Batch6", "Batch7", "Batch8", "Batch9", "Batch10") ~ case_when(
      Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group5" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group6" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
      Group == "Group7" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
      Group == "Group8" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
    )
  ))
sims$time <- sc_data$time
# sc_data %>%
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)  # Adjust this line based on actual column names from time categories
saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", 
                            "sims_", seed, "_L9.rds"))

cat("Generation Completed for L9.\n")
