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

##################################### cellType3 #####################################
# 10 batches, 5000 genes, 3000 cells each batch, 3 cell tupes
# each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 10
nCellType <- 3
batchCells <- batchCells <- rep(3000, nsample)
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
                    nGenes = 5000, batchCells=batchCells,
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

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_neg_L1_c3.rds"))
cat("Generation Completed for c3 neg.\n")

# simulate positive control
sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group3" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time
saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L1_c3.rds"))
cat("Generation Completed for c3 pos-L1.\n")

sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n() * 0.6)), rep(2, n() - floor(n() * 0.6)))),  # Adjust Group 1
    Group == "Group2" ~ sample(c(rep(1, ceiling(n() * 0.4)), rep(2, n() - ceiling(n() *0.4)))),  # Adjust Group 2
    Group == "Group3" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time
saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L2_c3.rds"))
cat("Generation Completed for c3 pos-L2.\n")

sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n() * 0.45)), rep(2, n() - floor(n() * 0.45)))),  # Adjust Group 1
    Group == "Group2" ~ sample(c(rep(1, ceiling(n() * 0.55)), rep(2, n() - ceiling(n() *0.55)))),  # Adjust Group 2
    Group == "Group3" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time
saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L3_c3.rds"))
cat("Generation Completed for c3 pos-L3.\n")

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n() * 0.495)), rep(2, n() - floor(n() * 0.495)))),  # Adjust Group 1
    Group == "Group2" ~ sample(c(rep(1, ceiling(n() * 0.505)), rep(2, n() - ceiling(n() *0.505)))),  # Adjust Group 2
    Group == "Group3" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time
saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L4_c3.rds"))
cat("Generation Completed for c3 pos-L4.\n")

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

# prop.table(table(sc_data$Group, sc_data$time), margin = 2)

##################################### cellType 5 ###################################
nsample <- 10
nCellType <- 5
batchCells <- batchCells <- rep(3000, nsample)
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
                    nGenes = 5000, batchCells=batchCells,
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

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_neg_L1_c5.rds"))
cat("Generation Completed for c5 neg.\n")

# simulate positive control wtih different effect size
sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group5" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time
# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)
saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L1_c5.rds"))
cat("Generation Completed for c5 pos-L1.\n")


sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group3" ~ sample(c(rep(1, floor(n() * 0.55)), rep(2, n() - floor(n() * 0.55)))),  
    Group == "Group4" ~ sample(c(rep(1, ceiling(n() * 0.45)), rep(2, n() - ceiling(n() *0.45)))), 
    Group == "Group5" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time
# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)
# 
saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L2_c5.rds"))
cat("Generation Completed for c5 pos-L2.\n")


sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n() * 0.45)), rep(2, n() - floor(n() * 0.45)))),  
    Group == "Group2" ~ sample(c(rep(1, ceiling(n() * 0.55)), rep(2, n() - ceiling(n() *0.55)))),
    Group == "Group3" ~ sample(c(rep(1, floor(n() * 0.55)), rep(2, n() - floor(n() * 0.55)))),  
    Group == "Group4" ~ sample(c(rep(1, ceiling(n() * 0.45)), rep(2, n() - ceiling(n() *0.45)))), 
    Group == "Group5" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L3_c5.rds"))
cat("Generation Completed for c5 pos-L3.\n")



sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n() * 0.495)), rep(2, n() - floor(n() * 0.495)))),  
    Group == "Group2" ~ sample(c(rep(1, ceiling(n() * 0.505)), rep(2, n() - ceiling(n() *0.505)))),  
    Group == "Group3" ~ sample(c(rep(1, floor(n() * 0.49)), rep(2, n() - floor(n() * 0.49)))),  
    Group == "Group4" ~ sample(c(rep(1, ceiling(n() * 0.50)), rep(2, n() - ceiling(n() *0.50)))), 
    Group == "Group5" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time

# colData(sims_neg) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L4_c5.rds"))
cat("Generation Completed for c5 pos-L4.\n")



##################################### cellType 8 ##################################
nsample <- 10
nCellType <- 9
batchCells <- batchCells <- rep(3000, nsample)
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
                    nGenes = 5000, batchCells=batchCells,
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

saveRDS(sims_neg, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_neg_L1_c9.rds"))
cat("Generation Completed for c9 neg.\n")

# simulate positive control wtih different effect size
sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group5" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group6" ~ sample(c(rep(1, 3 * floor(n()/4)), rep(2, floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group7" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group8" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group9" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L1_c9.rds"))
cat("Generation Completed for c9 pos-L1.\n")

sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n()/3)), rep(2, 2 * floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group2" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group3" ~ sample(c(rep(1, 2 * floor(n()/3)), rep(2, floor(n()/3)), rep(1, n() %% 3))),
    Group == "Group4" ~ sample(c(rep(1, floor(n()/4)), rep(2, 3 * floor(n()/4)), rep(1, n() %% 4))),
    Group == "Group5" ~ sample(c(rep(1, floor(n() * 0.55)), rep(2, n() - floor(n() * 0.55)))),  
    Group == "Group6" ~ sample(c(rep(1, ceiling(n() * 0.45)), rep(2, n() - ceiling(n() *0.45)))), 
    Group == "Group7" ~ sample(c(rep(1, floor(n() * 0.495)), rep(2, n() - floor(n() * 0.495)))),  
    Group == "Group8" ~ sample(c(rep(1, ceiling(n() * 0.505)), rep(2, n() - ceiling(n() *0.505)))), 
    Group == "Group9" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time
# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)
# 
saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L2_c9.rds"))
cat("Generation Completed for c0 pos-L2.\n")


sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n() * 0.4)), rep(2, n() - floor(n() * 0.4)))),  
    Group == "Group2" ~ sample(c(rep(1, ceiling(n() * 0.6)), rep(2, n() - ceiling(n() *0.6)))),
    Group == "Group3" ~ sample(c(rep(1, floor(n() * 0.55)), rep(2, n() - floor(n() * 0.55)))),  
    Group == "Group4" ~ sample(c(rep(1, ceiling(n() * 0.45)), rep(2, n() - ceiling(n() *0.45)))), 
    Group == "Group5" ~ sample(c(rep(1, floor(n() * 0.55)), rep(2, n() - floor(n() * 0.55)))),  
    Group == "Group6" ~ sample(c(rep(1, ceiling(n() * 0.45)), rep(2, n() - ceiling(n() *0.45)))), 
    Group == "Group7" ~ sample(c(rep(1, floor(n() * 0.495)), rep(2, n() - floor(n() * 0.495)))),  
    Group == "Group8" ~ sample(c(rep(1, ceiling(n() * 0.505)), rep(2, n() - ceiling(n() *0.505)))), 
    Group == "Group9" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L3_c9.rds"))
cat("Generation Completed for c9 pos-L3.\n")



sims_pos <- sims
sc_data_pos <- colData(sims_pos) %>% 
  as.data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = case_when(
    Group == "Group1" ~ sample(c(rep(1, floor(n() * 0.495)), rep(2, n() - floor(n() * 0.495)))),  
    Group == "Group2" ~ sample(c(rep(1, ceiling(n() * 0.505)), rep(2, n() - ceiling(n() *0.505)))),  
    Group == "Group3" ~ sample(c(rep(1, floor(n() * 0.49)), rep(2, n() - floor(n() * 0.49)))),  
    Group == "Group4" ~ sample(c(rep(1, ceiling(n() * 0.50)), rep(2, n() - ceiling(n() *0.50)))), 
    Group == "Group5" ~ sample(c(rep(1, floor(n() * 0.48)), rep(2, n() - floor(n() * 0.48)))),  
    Group == "Group6" ~ sample(c(rep(1, ceiling(n() * 0.52)), rep(2, n() - ceiling(n() *0.52)))), 
    Group == "Group7" ~ sample(c(rep(1, floor(n() * 0.495)), rep(2, n() - floor(n() * 0.495)))),  
    Group == "Group8" ~ sample(c(rep(1, ceiling(n() * 0.505)), rep(2, n() - ceiling(n() *0.505)))), 
    Group == "Group9" ~ sample(c(rep(1, ceiling(n()/2)), rep(2, floor(n()/2))))
  ))
sims_pos$time <- sc_data_pos$time

# colData(sims_pos) %>% as.data.frame() %>% count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
#   mutate(ratio = `1` / `2`)

saveRDS(sims_pos, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V4/sims/", 
                                "sims_", seed, "_pos_L4_c9.rds"))
cat("Generation Completed for c9 pos-L4.\n")
