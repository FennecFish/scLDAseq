setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(splatter)
# library(Seurat)
set.seed(1)

################ Single Sample Simulation ##################
# simulate five cell types
params <- newSplatParams()
params <- setParams(params, group.prob = c(0.3,0.05,0.15,0.1,0.4),
                    de.prob = c(0.1, 0.2, 0.3, 0.2, 0.05), 
                    nGenes = 500, batchCells=c(400,600),  
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = TRUE)

# create change of proportions
sim1 <- sims
table(sims$Group)
g1_cell <- sample(sim1$Cell[sim1$Group=="Group1"],300,replace=FALSE)
g2_cell <- sample(sim1$Cell[sim1$Group=="Group2"],5,replace=FALSE)
g3_cell <- sample(sim1$Cell[sim1$Group=="Group3"],50,replace=FALSE)
g4_cell <- sample(sim1$Cell[sim1$Group=="Group4"],30,replace=FALSE)
g5_cell <- sample(sim1$Cell[sim1$Group=="Group5"],200,replace=FALSE)
sim1_cell <- c(g1_cell,g2_cell, g3_cell, g4_cell, g5_cell)
sim1 <- sim1[,colnames(sim1) %in% sim1_cell]
sim1$time <- 1
print(prop.table(table(sim1$Group)))

sim2 <- sims
sim2 <- sim2[,!colnames(sim2) %in% sim1_cell]
sim2$time <- 2
print(prop.table(table(sim2$Group)))

# combine the two timepoints
sim_new_p1 <- cbind(sim1,sim2)  
sim_new_p1$time <- factor(sim_new_p1$time)

# save the data
saveRDS(sim_new_p1, "data/sim_single_sample.rds")
