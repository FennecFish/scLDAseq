setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(splatter)
library(scMerge)
library(scran)
# library(Seurat)
set.seed(1)

###### read estimated parameters ### 
params <- readRDS("data/est_params.rds")


##### estimating parameters #####
# estimate from real data
p2 <- readRDS("data/p10_p15_sub.rds")
dec.p2 <- modelGeneVar(p2)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, prop=0.3)
p2.sub <- p2[p2.chosen,]
p2.sub <- p2.sub[,colSums(assay(p2.sub)) > 500 & colSums(assay(p2.sub)) < 5000]
assay(p2.sub) <- as.matrix(assay(p2.sub))

params <- splatEstimate(p2.sub)
saveRDS(params, file = "data/est_params.rds")

p10 <- p2.sub[,colData(p2.sub)$patient_id == "BIOKEY_10"]
p15 <- p2.sub[,colData(p2.sub)$patient_id == "BIOKEY_15"]

params_p10 <- splatEstimate(p10)
saveRDS(params_p10, file = "data/est_params_p10.rds")

params_p15 <- splatEstimate(p15)
saveRDS(params_p15, file = "data/est_params_p15.rds")

######### simulating batches as sampels ##########
params <- setParams(params, group.prob = c(0.3,0.05,0.15,0.1,0.4),
                    de.prob = c(0.1, 0.2, 0.3, 0.2, 0.05), 
                    nGenes = 1000, batchCells=c(5500,4500))
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = TRUE)

# # create change of proportions
# sim1 <- sims
# table(sims$Group)
# g1_cell <- sample(sim1$Cell[sim1$Group=="Group1"],2000,replace=FALSE)
# g2_cell <- sample(sim1$Cell[sim1$Group=="Group2"],200,replace=FALSE)
# g3_cell <- sample(sim1$Cell[sim1$Group=="Group3"],800,replace=FALSE)
# g4_cell <- sample(sim1$Cell[sim1$Group=="Group4"],730,replace=FALSE)
# g5_cell <- sample(sim1$Cell[sim1$Group=="Group5"],1000,replace=FALSE)
# sim1_cell <- c(g1_cell,g2_cell, g3_cell, g4_cell, g5_cell)
# sim1 <- sim1[,colnames(sim1) %in% sim1_cell]
# sim1$time <- 1
# print(prop.table(table(sim1$Group)))
# 
# sim2 <- sims
# sim2 <- sim2[,!colnames(sim2) %in% sim1_cell]
# sim2$time <- 2
# print(prop.table(table(sim2$Group)))
# 
# # combine the two timepoints
# sim_new <- cbind(sim1,sim2)  
# sim_new$time <- factor(sim_new$time)
# 
# # save the data
# saveRDS(sim_new, "data/sim_2samples.rds")

cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.1,0.05,0.15,0.1, 0.6))
pre_count <- cell_count %*% t(pre_prp)
colnames(pre_count) <- paste0("Group",1:5)
rownames(pre_count) <- paste0("Batch", 1:2)

sampled_data <- colData(sims) %>%
    data.frame() %>%
    group_by(Group, Batch) %>%
    mutate(time = 2) %>% 
    ungroup() 

for (i in 1:nrow(pre_count)) {
    batch_name <- rownames(pre_count)[i]
    for (j in 1:ncol(pre_count)) {
        group_name <- colnames(pre_count)[j]
        sampled_data <- sampled_data %>%
            group_by(Group, Batch) %>%
            mutate(time = ifelse(
                Group == group_name & Batch == batch_name & 
                    row_number() %in% sample(row_number(), min(pre_count[i,j], n())), 1, time)) %>%
            ungroup() 
    }
}
sims$time <- sampled_data$time
##### simulating patient 10 #####
params_p10 <- setParams(params_p10, group.prob = c(0.3,0.05,0.15,0.1,0.4),
                    de.prob = c(0.1, 0.2, 0.3, 0.2, 0.05))
sims_p10 <- splatSimulate(params_p10, method = "groups",
                      verbose = FALSE)

sim1 <- sims_p10
table(sim1$Group)
g1_cell <- sample(sim1$Cell[sim1$Group=="Group1"],1300,replace=FALSE)
g2_cell <- sample(sim1$Cell[sim1$Group=="Group2"],100,replace=FALSE)
g3_cell <- sample(sim1$Cell[sim1$Group=="Group3"],400,replace=FALSE)
g4_cell <- sample(sim1$Cell[sim1$Group=="Group4"],250,replace=FALSE)
g5_cell <- sample(sim1$Cell[sim1$Group=="Group5"],300,replace=FALSE)
sim1_cell <- c(g1_cell,g2_cell, g3_cell, g4_cell, g5_cell)
sim1 <- sim1[,colnames(sim1) %in% sim1_cell]
sim1$time <- 1
print(prop.table(table(sim1$Group)))

sim2 <- sims_p10
sim2 <- sim2[,!colnames(sim2) %in% sim1_cell]
sim2$time <- 2
print(prop.table(table(sim2$Group)))

# combine the two timepoints
sim_new_p10 <- cbind(sim1,sim2)  
sim_new_p10$time <- factor(sim_new_p10$time)

colnames(sim_new_p10) <- paste0("Sample10_",colnames(sim_new_p10))
colData(sim_new_p10)$sample <- "sample10"

# save the data
saveRDS(sim_new_p10, "data/sim_p10.rds")

##### simulating patient 15 #####
params_p15 <- setParams(params_p15, group.prob = c(0.3,0.05,0.15,0.1,0.4),
                        de.prob = c(0.1, 0.2, 0.3, 0.2, 0.05))
sims_p15 <- splatSimulate(params_p15, method = "groups",
                          verbose = FALSE)

sim1 <- sims_p15
table(sim1$Group)
g1_cell <- sample(sim1$Cell[sim1$Group=="Group1"],1000,replace=FALSE)
g2_cell <- sample(sim1$Cell[sim1$Group=="Group2"],100,replace=FALSE)
g3_cell <- sample(sim1$Cell[sim1$Group=="Group3"],400,replace=FALSE)
g4_cell <- sample(sim1$Cell[sim1$Group=="Group4"],350,replace=FALSE)
g5_cell <- sample(sim1$Cell[sim1$Group=="Group5"],800,replace=FALSE)
sim1_cell <- c(g1_cell,g2_cell, g3_cell, g4_cell, g5_cell)
sim1 <- sim1[,colnames(sim1) %in% sim1_cell]
sim1$time <- 1
print(prop.table(table(sim1$Group)))

sim2 <- sims_p15
sim2 <- sim2[,!colnames(sim2) %in% sim1_cell]
sim2$time <- 2
print(prop.table(table(sim2$Group)))

# combine the two timepoints
sim_new_p15 <- cbind(sim1,sim2)  
sim_new_p15$time <- factor(sim_new_p15$time)

colnames(sim_new_p15) <- paste0("Sample15_",colnames(sim_new_p15))
colData(sim_new_p15)$sample <- "sample2"
# save the data
saveRDS(sim_new_p15, "data/sim_p15.rds")


#### combined data ####
im_new_p10
cbind(sim_new_p10, sim_new_p15)

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


################ Multi-sample Simulation ##################
params <- newSplatParams()
params <- setParams(params, group.prob = c(0.1,0.4,0.05,0.25,0.2),
                    de.prob = c(0.1, 0.2, 0.3, 0.2, 0.2), 
                    nGenes = 1000, batchCells=c(1500,2500),  
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = TRUE)

# create change of proportions
sim1 <- sims
table(sim1$Group)
g1_cell <- sample(sim1$Cell[sim1$Group=="Group1"],10,replace=FALSE)
g2_cell <- sample(sim1$Cell[sim1$Group=="Group2"],400,replace=FALSE)
g3_cell <- sample(sim1$Cell[sim1$Group=="Group3"],40,replace=FALSE)
g4_cell <- sample(sim1$Cell[sim1$Group=="Group4"],150,replace=FALSE)
g5_cell <- sample(sim1$Cell[sim1$Group=="Group5"],100,replace=FALSE)
sim1_cell <- c(g1_cell,g2_cell, g3_cell, g4_cell, g5_cell)
sim1 <- sim1[,colnames(sim1) %in% sim1_cell]
sim1$time <- 1
print(prop.table(table(sim1$Group)))

sim2 <- sims
sim2 <- sim2[,!colnames(sim2) %in% sim1_cell]
sim2$time <- 2
print(prop.table(table(sim2$Group)))

# combine the two timepoints
sim_new_p2 <- cbind(sim1,sim2)  
sim_new_p2$time <- factor(sim_new_p2$time)
saveRDS(sim_new_p2, file = "data/sim_p2_sample.rds")

p2 <- readRDS("data/sim_p2_sample.rds")
p1 <- readRDS("data/sim_single_sample.rds")

colnames(p2) <- paste0("Sample2_",colnames(p2))
colData(p2)$sample <- "sample2"
colnames(p1) <- paste0("Sample1_",colnames(p1))
colData(p1)$sample <- "sample2"

cbind(p1,p2)
p1sce_list <- list(p1,p2)
sce_combo <- sce_cbind(sce_list, method = "intersect", 
          exprs = c("BatchCellMeans", "BaseCellMeans",
                    "counts"),
          colData_names =c("Group,time,sample"))

saveRDS()