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

params <- newSplatParams()
params <- setParams(params, group.prob = c(0.3,0.05,0.15,0.1,0.4),
                    de.prob = c(0.1, 0.2, 0.3, 0.2, 0.05), 
                    nGenes = 5000, batchCells=c(5500,4500, 6000))
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)# create change of proportions

# we assume that the cells pre and post treatment are equal
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.1,0.05,0.15,0.1, 0.6))
pre_count <- cell_count %*% t(pre_prp)
colnames(pre_count) <- paste0("Group",1:5)
rownames(pre_count) <- paste0("Batch", 1:3)

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
    
time_prop <- sampled_data %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))


sims$time <- sampled_data$time
save(sims, file = "data/sims_3samples_5groups_same_direction.rds")

#### feature selection #####
sims <- scuttle::logNormCounts(sims)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1000)
sims.sub <- sims[p2.chosen,]
saveRDS(sims.sub, file = "data/sims_3samples_5groups_same_direction_1000genes.rds")
# p2.sub <- p2.sub[,colSums(assay(p2.sub)) > 500 & colSums(assay(p2.sub)) < 5000]
# assay(p2.sub) <- as.matrix(assay(p2.sub))

##### eval 
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims.sub)
K <- 5
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

# png("../res/2sample_5type_true_prop.png", height = 2000, width = 2000, res = 300)
ggplot(time_prop, aes(x = Group, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type True Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") + 
    facet_grid(~Batch)
# dev.off()

max_indices <- apply(res$theta, 1, which.max)
colnames(res$theta) <- paste0("topic_", 1:ncol(res$theta))
res_cluster <- colnames(res$theta)[max_indices]
adjustedRandIndex(res_cluster,sims$Group) # 0.7667998

res_dat <- colData(sims) %>% 
    data.frame() %>% 
    mutate(assigned_cluster = res_cluster)

res_prop <- res_dat %>% 
    group_by(time, Batch) %>%
    count(assigned_cluster) %>%
    mutate(Proportion = n / sum(n))

# png("../res/2sample_5type_assigned_prop.png", height = 2000, width = 2000, res = 300)
ggplot(res_prop, aes(x = assigned_cluster, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "assigned_cluster", y = "Proportion", fill = "time", 
         title = "Cell Type Assigned Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    facet_grid(~Batch)
# dev.off()

############### Using Seurat ###########################
library(Seurat)
# counts <- assay(sims, "counts")
# libsizes <- colSums(counts)
# size.factors <- libsizes/mean(libsizes)
# logcounts(sims) <- log2(t(t(counts)/size.factors) + 1)
seurat.sims <- as.Seurat(sims, counts = "counts", data = "logcounts")
# 
# seurat.sims <- NormalizeData(seurat.sims, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.sims)
seurat.sims <- ScaleData(seurat.sims, features = all.genes)
seurat.sims <- RunPCA(seurat.sims)

seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)
adjustedRandIndex(Idents(seurat.sims), sims$Group) # 0.4218298

########################## ########################## ########################## 
########################## opposite diretions ###########################################
########################## ########################## ########################## 
sims.sub <- readRDS("data/sims_3samples_5groups_same_direction_1000genes.rds")
sampled_data_sub <- colData(sims.sub) %>%
    data.frame() %>%
    mutate(new_time = ifelse(time == 1 & Batch == "Batch2", 2, 1)) %>%
    mutate(time = ifelse(Batch == "Batch2", new_time, time))
sims.sub$time <- sampled_data_sub$time

##### eval 
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims.sub)
K <- 5
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
time_prop <- colData(sims.sub) %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))

# png("../res/2sample_5type_true_prop.png", height = 2000, width = 2000, res = 300)
ggplot(time_prop, aes(x = Group, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type True Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") + 
    facet_grid(~Batch)
# dev.off()

max_indices <- apply(res$theta, 1, which.max)
colnames(res$theta) <- paste0("topic_", 1:ncol(res$theta))
res_cluster <- colnames(res$theta)[max_indices]
adjustedRandIndex(res_cluster,sims.sub$Group) # 0.7665575

res_dat <- colData(sims.sub) %>% 
    data.frame() %>% 
    mutate(assigned_cluster = res_cluster)

res_prop <- res_dat %>% 
    group_by(time, Batch) %>%
    count(assigned_cluster) %>%
    mutate(Proportion = n / sum(n))

# png("../res/2sample_5type_assigned_prop.png", height = 2000, width = 2000, res = 300)
ggplot(res_prop, aes(x = assigned_cluster, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "assigned_cluster", y = "Proportion", fill = "time", 
         title = "Cell Type Assigned Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    facet_grid(~Batch)
# dev.off()

############### Using Seurat ###########################
library(Seurat)
# counts <- assay(sims, "counts")
# libsizes <- colSums(counts)
# size.factors <- libsizes/mean(libsizes)
# logcounts(sims) <- log2(t(t(counts)/size.factors) + 1)
seurat.sims <- as.Seurat(sims.sub, counts = "counts", data = "logcounts")
# 
# seurat.sims <- NormalizeData(seurat.sims, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.sims)
seurat.sims <- ScaleData(seurat.sims, features = all.genes)
seurat.sims <- RunPCA(seurat.sims)

seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)
adjustedRandIndex(Idents(seurat.sims), sims$Group) # 0.4218298
