### This script serves as a positive control for the scLDAseq 

setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(mclust)
# library(Seurat)
set.seed(1)
##################################################
################ Three Patients ##################
##################################################
params <- newSplatParams()
params <- setParams(params, group.prob = c(0.5,0.5),
                    de.prob = c(0.3, 0.3), 
                    nGenes = 1000, batchCells=c(1000,1500, 800),  
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)# create change of proportions

# we assume that the cells pre and post treatment are equal
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.5,0.5))
pre_count <- cell_count %*% t(pre_prp)


sampled_data <- colData(sims) %>%
    data.frame() %>%
    group_by(Group, Batch) %>%
    mutate(time = case_when(
        Group == "Group1" & Batch == "Batch1" & 
            row_number() %in% sample(row_number(), min(pre_count[1,1], n())) ~ 1,
        Group == "Group2" & Batch == "Batch1" & 
            row_number() %in% sample(row_number(), min(pre_count[1,2], n())) ~ 1,
        Group == "Group1" & Batch == "Batch2" & 
            row_number() %in% sample(row_number(), min(pre_count[2,1], n())) ~ 1,
        Group == "Group2" & Batch == "Batch2" & 
            row_number() %in% sample(row_number(), min(pre_count[2,2], n())) ~ 1,
        Group == "Group1" & Batch == "Batch3" & 
            row_number() %in% sample(row_number(), min(pre_count[3,1], n())) ~ 1,
        Group == "Group2" & Batch == "Batch3" & 
            row_number() %in% sample(row_number(), min(pre_count[3,2], n())) ~ 1,
        TRUE ~ 2 
    )) %>%
    ungroup() 


time_prop <- sampled_data %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))


sims$time <- sampled_data$time

##### eval 
r.file <- paste0("R/",list.files("R/"))
# r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
# sourceCpp("../stm/src/STMCfuns.cpp")
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims)
K <- 2 
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
data <- dat$meta
sample <- "Batch"

res <- multi_stm(documents = documents, vocab = vocab,
                 K = K, prevalence = prevalence, content = NULL,
                 data = data, 
                 sce = sce,
                 sample = sample,
                 init.type= "Spectral",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 control = list(gamma.maxits=3000))

# access the true proportion distribution
time_prop <- colData(sims) %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))

png("../res/neg_control_3sample_true_prop.png", height = 2000, width = 2000, res = 300)
ggplot(time_prop, aes(x = Group, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type True Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") + 
    facet_grid(~Batch)
dev.off()

max_indices <- apply(res$theta, 1, which.max)
colnames(res$theta) <- paste0("topic_", 1:ncol(res$theta))
res_cluster <- colnames(res$theta)[max_indices]
adjustedRandIndex(res_cluster,sims$Group) # 1

res_dat <- colData(sims) %>% 
    data.frame() %>% 
    mutate(assigned_cluster = res_cluster)

res_prop <- res_dat %>% 
    group_by(time, Batch) %>%
    count(assigned_cluster) %>%
    mutate(Proportion = n / sum(n))

png("../res/neg_control_3sample_assigned_prop.png", height = 2000, width = 2000, res = 300)
ggplot(res_prop, aes(x = assigned_cluster, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "assigned_cluster", y = "Proportion", fill = "time", 
         title = "Cell Type Assigned Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    facet_grid(~Batch)
dev.off()

##################################################
################### Two Patients #################
##################################################
params <- newSplatParams()
params <- setParams(params, group.prob = c(0.5,0.5),
                    de.prob = c(0.3, 0.3), 
                    nGenes = 1000, batchCells=c(1000,1500),  
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)# create change of proportions

# we assume that the cells pre and post treatment are equal
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.5,0.5))
pre_count <- cell_count %*% t(pre_prp)


sampled_data <- colData(sims) %>%
    data.frame() %>%
    group_by(Group, Batch) %>%
    mutate(time = case_when(
        Group == "Group1" & Batch == "Batch1" & 
            row_number() %in% sample(row_number(), min(pre_count[1,1], n())) ~ 1,
        Group == "Group2" & Batch == "Batch1" & 
            row_number() %in% sample(row_number(), min(pre_count[1,2], n())) ~ 1,
        Group == "Group1" & Batch == "Batch2" & 
            row_number() %in% sample(row_number(), min(pre_count[2,2], n())) ~ 1,
        Group == "Group2" & Batch == "Batch2" & 
            row_number() %in% sample(row_number(), min(pre_count[2,1], n())) ~ 1,
        TRUE ~ 2 
    )) %>%
    ungroup() 

time_prop <- sampled_data %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))


sims$time <- sampled_data$time

##### eval 
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims)
K <- 2 
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
data <- dat$meta
sample <- "Batch"

res <- multi_stm(documents = documents, vocab = vocab,
                 K = K, prevalence = prevalence, content = NULL,
                 data = data, 
                 sce = sce,
                 sample = sample,
                 init.type= "Spectral",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 control = list(gamma.maxits=3000))

# access the true proportion distribution
time_prop <- colData(sims) %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))

png("../res/neg_control_2sample_true_prop.png", height = 2000, width = 2000, res = 300)
ggplot(time_prop, aes(x = Group, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type True Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") + 
    facet_grid(~Batch)
dev.off()

max_indices <- apply(res$theta, 1, which.max)
colnames(res$theta) <- paste0("topic_", 1:ncol(res$theta))
res_cluster <- colnames(res$theta)[max_indices]
adjustedRandIndex(res_cluster,sims$Group) # 1

res_dat <- colData(sims) %>% 
    data.frame() %>% 
    mutate(assigned_cluster = res_cluster)

res_prop <- res_dat %>% 
    group_by(time, Batch) %>%
    count(assigned_cluster, time, Batch) %>%
    mutate(Proportion = n / sum(n))

png("../res/neg_control_2sample_assigned_prop.png", height = 2000, width = 2000, res = 300)
ggplot(res_prop, aes(x = assigned_cluster, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "assigned_cluster", y = "Proportion", fill = "time", 
         title = "Cell Type Assigned Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    facet_grid(~Batch)
dev.off()

##################################################
########## Single Patient Simulation ######## ####
##################################################
params <- newSplatParams()
params <- setParams(params, group.prob = c(0.5,0.5),
                    de.prob = c(0.3, 0.3), 
                    nGenes = 1000, batchCells=c(1000,1500),  
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = TRUE)# create change of proportions

# we assume that the cells pre and post treatment are equal
cell_count <- sum(table(sims$Group))/length(unique(sims$Group))
pre_prp <- c(0.5,0.5)
pre_count <- cell_count* pre_prp

sampled_data <- colData(sims) %>%
    data.frame() %>%
    group_by(Group) %>%
    mutate(time = case_when(
        Group == "Group1" & row_number() %in% sample(row_number(), min(pre_count[1], n())) ~ 1,
        Group == "Group2" & row_number() %in% sample(row_number(), min(pre_count[2], n())) ~ 1,
        TRUE ~ 2  # Assuming 'time' already exists and is numeric; otherwise, adjust as needed
    )) %>%
    ungroup() 

prop.table(table(sampled_data$time,sampled_data$Group), margin = 2)
sims$time <- sampled_data$time

##### eval 
r.file <- paste0("R/",list.files("R/"))
# r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
# sourceCpp("../stm/src/STMCfuns.cpp")
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims)
K <- 2 
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
data <- dat$meta

res <- multi_stm(documents = documents, vocab = vocab,
                 K = K, prevalence = prevalence, content = NULL,
                 data = data, 
                 sce = sce,
                 sample = NULL,
                 init.type= "Spectral",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 control = list(gamma.maxits=3000))

# access the true proportion distribution
time_prop <- colData(sims) %>% 
    data.frame() %>%
    group_by(time) %>%
    count(Group, time) %>%
    mutate(Proportion = n / sum(n))

png("../res/neg_control_1sample_true_prop.png", height = 2000, width = 2000, res = 300)
ggplot(time_prop, aes(x = Group, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type True Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") 
dev.off()

max_indices <- apply(res$theta, 1, which.max)
colnames(res$theta) <- paste0("topic_", 1:ncol(res$theta))
res_cluster <- colnames(res$theta)[max_indices]
adjustedRandIndex(res_cluster,sims$Group) # 1

res_dat <- colData(sims) %>% data.frame() %>% mutate(assigned_cluster = res_cluster)
res_prop <- res_dat %>% 
    group_by(time) %>%
    count(assigned_cluster, time) %>%
    mutate(Proportion = n / sum(n))

png("../res/neg_control_1sample_assigned_prop.png", height = 2000, width = 2000, res = 300)
ggplot(res_prop, aes(x = assigned_cluster, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "assigned_cluster", y = "Proportion", fill = "time", 
         title = "Cell Type Assigned Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
dev.off()
