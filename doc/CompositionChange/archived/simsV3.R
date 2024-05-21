setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(tidyverse)
library(splatter)
library(SingleCellExperiment)
library(tidyr)


args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))

##################################### level 1 #####################################
# One batch, 5000 genes, 3000 cells, 3 cell
# types, each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 1
nCellType <- 3
batchCells <- 3000
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
                    nGenes = 5000, batchCells=batchCells,
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

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3/sims/", 
                                "sims_", seed, "_neg_c3.rds"))
cat("Generation Completed for c3 neg.\n")

# simulate positive control
sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group3" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3/sims/", 
                                "sims_", seed, "_pos_c3.rds"))
cat("Generation Completed for c3 pos.\n")
# prop.table(table(sc_data$Group, sc_data$time), margin = 2)

##################################### level 3 ###################################
nsample <- 1
nCellType <- 5
batchCells <- 3000
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
                    nGenes = 5000, batchCells=batchCells,
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

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3/sims/", 
                                "sims_", seed, "_neg_c5.rds"))
cat("Generation Completed for c5 neg.\n")
# simulate positive control
sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group5" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3/sims/", 
                                "sims_", seed, "_pos_c5.rds"))

cat("Generation Completed for c5 pos.\n")
# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

##################################### level 3 ##################################
nsample <- 1
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

# first simulate negative control: no change
sims_neg <- sims
sc_data_neg <- colData(sims_neg) %>% 
  as.data.frame() %>%
  group_by(Group) %>%
  mutate(time = sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2)))))
sims_neg$time <- sc_data_neg$time

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3/sims/", 
                                "sims_", seed, "_neg_c8.rds"))
cat("Generation Completed for c8 neg.\n")
# simulate positive control
sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group5" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group6" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group7" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group8" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3/sims/", 
                                "sims_", seed, "_pos_c8.rds"))
cat("Generation Completed for c8 pos.\n")
