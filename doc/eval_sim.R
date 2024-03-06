setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
set.seed(1)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(mclust)
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

# two samples, use batches as samples
sce <- readRDS("data/sim_2samples.rds")
dat <- prepsce(sce)
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

saveRDS(res, file = "data/res_batch_as_sample.rds")
png("../res/batch_as_sample_convergence.png", height = 1000, width = 1000, res = 300)
plot(res$convergence$bound, type = "l", ylab = "Approximate Objective", main = "Convergence")
dev.off()

# access the true proportion distribution
sim <- readRDS("data/sim_2samples.rds")
time_prop <- colData(sim) %>% 
    data.frame() %>%
    group_by(Batch, time) %>%
    count(Group, time, Batch) %>%
    mutate(Proportion = n / sum(n))
# prop.table(table(sim$time,sim$Group, sim$Batch),margin=1)
# time_prop <- prop.table(table(sim$time,sim$Group, sim$Batch),margin=1) %>% 
#     data.frame() %>%
#     mutate(time = paste0("Time",Var1), sample = paste0("sample_",gsub("[^0-9]", "", Var3))) 


png("../res/batch_as_sample_true_prop.png", height = 2000, width = 2000, res = 300)
ggplot(time_prop, aes(x = Group, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type True Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") + 
    facet_grid(~Batch)
dev.off()

# library(topicmodels)
# library(cluster)
# dist <- distHellinger(res$theta)
# dist <- as.dist(dist)
# kc <- pam(dist, K)
# kc_res <- kc$clustering
# names(kc_res) <- rownames(res$settings$covariates$X)

max_indices <- apply(res$theta, 1, which.max)
colnames(res$theta) <- paste0("topic_", 1:ncol(res$theta))
res_cluster <- colnames(res$theta)[max_indices]
adjustedRandIndex(res_cluster,sims$Group) # -0.0009481646



res_dat <- colData(sim) %>% data.frame() %>% mutate(assigned_cluster = kc_res)
res_prop <- res_dat %>% 
    group_by(Batch, time) %>%
    count(assigned_cluster, time) %>%
    mutate(Proportion = n / sum(n))

png("../res/batch_as_sample_assigned_prop.png", height = 2000, width = 2000, res = 300)
ggplot(res_prop, aes(x = assigned_cluster, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "assigned_cluster", y = "Proportion", fill = "time", 
         title = "Cell Type Assigned Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") + 
    facet_grid(~Batch)
dev.off()

### Using Seurat
library(Seurat)
counts <- assay(sce, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)
seurat.sims <- as.Seurat(sce, counts = "counts", data = "logcounts")

seurat.sims <- NormalizeData(seurat.sims, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.sims)
seurat.sims <- ScaleData(seurat.sims, features = all.genes)
seurat.sims <- RunPCA(seurat.sims)

seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)
adjustedRandIndex(Idents(seurat.sims), sims$Group)



########### single sample ############
sce <- readRDS("data/sim_single_sample.rds")
prop.table(table(sce$Group, sce$time))
dat <- prepsce(sce)
K <- 5
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
data <- dat$meta
# sample <- "patient_id"


res <- multi_stm(documents = documents, vocab = vocab,
           K = K, prevalence = prevalence, content = NULL,
           data = data, 
           sce = sce,
           sample = NULL,
           init.type= "Spectral",
           gamma.prior= "Pooled",
           kappa.prior= "L1",
           control = list(gamma.maxits=3000))
plot(res$convergence$bound, type = "l", ylab = "Approximate Objective", main = "Convergence")

time_prop <- prop.table(table(sce$time,sce$Group),margin=1) %>% 
    data.frame() %>%
    mutate(time = paste0("Time",Var1))
# png("time proportion.png", width = 2500, height = 1500, res = 300)
ggplot(time_prop, aes(x = Var2, y = Freq, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type Proportion Distribution For Each Timepoint") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
#dev.off()

# png("batch proportion.png", width = 2500, height = 1500, res = 300)
# ggplot(batch_prop, aes(x = Var2, y = Freq, fill = as.factor(batch))) +
#     geom_bar(stat = "identity", position = "dodge") +
#     labs(x = "Group", y = "Proportion", fill = "batch", 
#          title = "Cell Type Proportion Distribution For Each Batch") +
#     theme_minimal() +
#     scale_fill_brewer(palette = "Set2")
# dev.off()

# genes for each topic
labelTopics(res,c(1:K))
# labelTopics(b)

# covariance effect
prep <-estimateEffect(1:K ~ time, stmobj = res, meta= data,uncertainty = "Global")

summary(prep, topics=1:K)

summary.estimateEffect(prep)
# png("estimate effect.png", width = 3500, height = 2500, res = 300)
plot(prep, "time", xlim = c(-0.2,0.2), 
     main = "Estimated Effect of Time as a Covariate to Cell-Topic Proportion",
     method="difference",cov.value1=1,cov.value2=2,
     model = res,
     labeltype = "frex")
# dev.off()