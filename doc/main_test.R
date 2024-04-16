setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
# setwd("/proj/milovelab/wu/scLDAseq")
set.seed(1)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(mclust)
library(lme4)
library(lmerTest)
library(dplyr)
library(tibble)
library(stats)
library(splatter)
library(scater)


###### simple simulation #########################
params <- newSplatParams()
params <- setParams(params, group.prob = c(0.3,0.3, 0.4),
                    de.prob = c(0.3, 0.3, 0.3),
                    nGenes = 500, batchCells=c(200,200,100,200),
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)# create change of proportions

cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.3,0.3, 0.4))
pre_count <- cell_count %*% t(pre_prp)
colnames(pre_count) <- paste0("Group",1:length(unique(sims$Group)))
rownames(pre_count) <- paste0("Batch", 1:length(unique(sims$Batch)))

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
saveRDS(sims, file = "data/toydat.rds")

sims <- readRDS("data/toydat.rds")

colData(sims) %>% as.data.frame() %>%
  count(Batch, time, Group) %>%  # Count combinations of Batch, time, and Group
  group_by(Batch, time) %>%  # Group by Batch to calculate proportions within each Batch
  mutate(Proportion = n / sum(n)) %>% 
  arrange(Batch, Group) %>%  # Calculate proportion of each combination within each Batch
  ungroup()  

#### QC ######
sims <- quickPerCellQC(sims)

sims <- scuttle::logNormCounts(sims)
library(scran)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1000)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

################### difficult simluation ###################
nsample <- 8
nCellType <- 5
batchCells <- rep(200, nsample)
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
                    nGenes = 2000, 
                    batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale)

sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# assume number of cells are equal across time
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/2)
# generate cell type proportion at time 1
all_rn_pre <- matrix(runif(nsample * nCellType), ncol = nsample)
pre_prp <- sweep(all_rn_pre, 2, colSums(all_rn_pre), FUN="/")
#calculate count for each group/batch at time 1
pre_count <- ceiling(sweep(pre_prp, 2, t(cell_count)[1, ], "*")) 
rownames(pre_count) <- paste0("Group",1:length(unique(sims$Group)))
colnames(pre_count) <- paste0("Batch", 1:length(unique(sims$Batch)))

sampled_data <- colData(sims) %>%
    data.frame() %>%
    group_by(Group, Batch) %>%
    mutate(time = 2) %>% 
    ungroup() 

# randomly assign time = 1 based on the count in pre_count
for (i in 1:ncol(pre_count)) {
    batch_name <- colnames(pre_count)[i]
    for (j in 1:nrow(pre_count)) {
        group_name <- rownames(pre_count)[j]
        sampled_data <- sampled_data %>%
            group_by(Group, Batch) %>%
            mutate(time = ifelse(
                Group == group_name & Batch == batch_name & 
                    row_number() %in% sample(row_number(), min(pre_count[j,i], n())), 1, time)) %>%
            ungroup() 
    }
}
# assign the new time
sims$time <- sampled_data$time

#### scLDAseq#############

# sims <-readRDS("data/sims_1712865809_L7.rds")
# sims <-readRDS("data/toydat.rds")


# library(sctransform)
# library(Seurat)
# sims <- scuttle::logNormCounts(sims)
# seurat <- CreateSeuratObject(counts = counts(sims))
# seurat <- SCTransform(seurat, verbose = FALSE)
sims <-readRDS("data/toydat.rds")
sims <- sims[rowSums(counts(sims))!=0,]

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

# sims <- readRDS("data/sims_1712873779_L3.rds")
# scSTM <- readRDS("data/scSTM_allgenes_noContent_1712873779_L3.rds")

K <- length(unique(sims$Group))
res <- multi_stm(sce = sims,
                 K = K, prevalence = ~time, content = NULL,
                 sample = "Batch",
                 init.type= "Spectral",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 seed = 123,
                 control = list(gamma.maxits=3000),
                 emtol=1e-5)

res <- selectModel(sce = sims,
                    K = K, prevalence = ~time, content = NULL,
                    sample = "Batch", N = 5, ts_runs = 5, random_run = 5)
all_values <- unlist(res$bound)
max_value <- max(all_values)
max_position_in_vector <- which(all_values == max_value)
res.final <- res$runout[[max_position_in_vector]]

### clustering
scSTMobj <- res.final
max_indices <- apply(scSTMobj$theta, 1, which.max)
colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
rownames(scSTMobj$theta) <- colnames(scSTMobj$mu$mu)
res_cluster <- colnames(scSTMobj$theta)[max_indices]
names(res_cluster) <- rownames(scSTMobj$theta)
adjustedRandIndex(res_cluster,sims$Group)

res <- multi_stm(sce = sims,
                 K = K, prevalence = ~time, content = NULL,
                 sample = "Batch",
                 init.type= "Random",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 seed = 9248714,
                 control = list(gamma.maxits=3000),
                 emtol=1e-6)

res <- multi_stm(sce = sims,
                 K = K, prevalence = ~time, content = NULL,
                 sample = "Batch",
                 init.type= "Random",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 seed = 1,
                 control = list(gamma.maxits=3000),
                 emtol=1e-4)


# dat <- prepsce(sims)
K <- length(unique(sims$Group))

res <- multi_stm(sce = sims,
                 K = K, prevalence = ~time, content = NULL,
                 sample = "Batch",
                 init.type= "Random",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 seed = 100,
                 control = list(gamma.maxits=3000),
                 emtol=1e-3)


test <- selectModel(sce = sims,
                    K = K, prevalence = ~time, content = NULL,
                    sample = "Batch", N = 10, ts_runs = 10, random_run = 10)


# verify with STM
dat <- prepsce(sims)
r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
sourceCpp("../stm/src/STMCfuns.cpp")

res.stm <- stm(documents = dat$documents, vocab = dat$vocab,
               K = K, prevalence = ~time, content = NULL,
               data = dat$meta, 
               init.type= "Random",
               gamma.prior= "Pooled",
               kappa.prior= "L1",
               emtol=1e-3,
               control = list(gamma.maxits=3000))

all_values <- unlist(test$bound)
max_value <- max(all_values)
max_position_in_vector <- which(all_values == max_value)

res <- test$runout[[max_position_in_vector]]

res <- multi_stm(sce = sims,
                          K = K, prevalence = ~time, content = NULL,
                          sample = "Batch",
                          init.type= "NMF",
                          gamma.prior= "Pooled",
                          kappa.prior= "L1",
                          control = list(gamma.maxits=3000))

saveRDS(res.scLDAseq, file = "data/scLDAseq_toydat_neg.rds")

max_indices <- apply(res$theta, 1, which.max)
colnames(res$theta) <- paste0("topic_", 1:ncol(res$theta))
rownames(res$theta) <- colnames(res$mu$mu)
res_cluster <- colnames(res$theta)[max_indices]
names(res_cluster) <- rownames(res$theta)
scSTM_cluster <- res_cluster[match(names(res_cluster), sims$Cell)]
adjustedRandIndex(scSTM_cluster, sims$Group)


# STM for reference
library(stm)


eff <- estimateEffect(1:3 ~ time, 
                      stmobj = res.stm, meta = dat$meta, uncertainty = "Global")
summary(eff)

# calculate adjRDI
max_indices <- apply(res.stm$theta, 1, which.max)
colnames(res.stm$theta) <- paste0("topic_", 1:ncol(res.stm$theta))
rownames(res.stm$theta) <- colnames(res.stm$mu$mu)
res_cluster <- colnames(res.stm$theta)[max_indices]
names(res_cluster) <- rownames(res.stm$theta)
scSTM_cluster <- res_cluster[match(names(res_cluster), sims$Cell)]
adjustedRandIndex(scSTM_cluster, sims$Group)

######### seurat##########

seurat.sims <- as.Seurat(sims, counts = "counts", data = "logcounts")
seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 500)
all.genes <- rownames(seurat.sims)
seurat.sims <- ScaleData(seurat.sims, features = all.genes)
seurat.sims <- RunPCA(seurat.sims)

seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)

seurat_cluster <- Idents(seurat.sims)[match(names(Idents(seurat.sims)), sims$Cell)]
adjustedRandIndex(seurat_cluster, sims$Group)

####### scProportion ##############
# sc_utils_obj = seurat.sims
seurat.sims$Time <- ifelse(seurat.sims$time == 1, "pre","on")
#seurat.sims$on <- ifelse(seurat.sims$time == 2, 1, 0)
library(scProportionTest)
# metadata <- as.data.table(seurat.sims[[]],keep.rownames = "cell_id")
sc_utils_obj <- sc_utils(seurat.sims)
cluster_identity = "Group"
sample_1 = "pre"
sample_2 = "on"
sample_identity = "Time"
n_permutations = 1000

## Prepare data.
meta_data <- copy(sc_utils_obj@meta_data)

meta_data <- meta_data[
  get(sample_identity) %in% c(sample_1, sample_2),
  c(..sample_identity, ..cluster_identity)
]

setnames(
  meta_data,
  old = c(sample_identity, cluster_identity),
  new = c("samples", "clusters")
)

meta_data[, clusters := as.character(clusters)]
cluster_cases <- unique(meta_data[["clusters"]])

## Get observed differences in fraction.
obs_diff <- meta_data[, .(count = .N), by = .(samples, clusters)]
obs_diff <- obs_diff[
  CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
  on = .(samples, clusters)
][
  is.na(count), count := 0
][]
obs_diff[, fraction := count / sum(count), by = samples]
obs_diff <- dcast(obs_diff, clusters ~ samples, value.var = "fraction")
obs_diff[, obs_log2FD := log2(get(sample_2)) - log2(get(sample_1))]
# log2(obs_diff$2) - log2(obs_diff$1)
## Permutation test.
perm_results <- matrix(NA, nrow(obs_diff), n_permutations)
rownames(perm_results) <- sort(cluster_cases)

for (i in seq_len(n_permutations)) {
  permuted <- copy(meta_data)
  permuted[["samples"]] <- sample(permuted[["samples"]])
  permuted <- permuted[, .(count = .N), by = .(samples, clusters)]
  permuted <- permuted[
    CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
    on = .(samples, clusters)
  ][
    is.na(count), count := 0
  ][]
  permuted[, fraction := count / sum(count), by = samples]
  permuted <- dcast(permuted, clusters ~ samples, value.var = "fraction")
  permuted[, perm_log2FD := log2(get(sample_2)) - log2(get(sample_1))]
  
  perm_results[, i] <- permuted[["perm_log2FD"]]
}

increased <- rowSums(apply(perm_results, 2, function(x) obs_diff[["obs_log2FD"]] <= x))
increased <- (increased + 1) / (n_permutations + 1)

decreased <- rowSums(apply(perm_results, 2, function(x) obs_diff[["obs_log2FD"]] >= x))
decreased <- (decreased + 1) / (n_permutations + 1)

obs_diff[, pval := ifelse(obs_log2FD > 0, increased[.I], decreased[.I])]
obs_diff[, FDR := p.adjust(pval, "fdr")]

## Boostrap log2FD CI.
boot_results <- matrix(NA, nrow(obs_diff), n_permutations)
rownames(boot_results) <- sort(cluster_cases)

for (i in seq_len(n_permutations)) {
  booted <- copy(meta_data)
  booted[, clusters := sample(clusters, replace = TRUE), by = samples]
  booted <- booted[, .(count = .N), by = .(samples, clusters)]
  booted <- booted[
    CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
    on = .(samples, clusters)
  ][
    is.na(count), count := 0
  ][]
  booted[, fraction := count / sum(count), by = samples]
  booted <- dcast(booted, clusters ~ samples, value.var = "fraction")
  booted[, boot_log2FD := log2(get(sample_2)) - log2(get(sample_1))]
  
  boot_results[, i] <- booted[["boot_log2FD"]]
}

boot_results[!is.finite(boot_results)] <- NA
boot_mean <- rowMeans(boot_results, na.rm = TRUE)
boot_ci <- t(apply(boot_results, 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))
boot_ci <- as.data.table(boot_ci)
setnames(boot_ci, old = c(1, 2), new = c("boot_CI_2.5", "boot_CI_97.5"))

obs_diff[, boot_mean_log2FD := boot_mean]
obs_diff <- cbind(obs_diff, boot_ci)




########## fishpond ######################

########### permutation test #################
# library(scProportionTest)
n_permutations = 1000
meta_data <- 
sc_utils_obj <- sc_utils(seurat.sims)
cluster_identity = "Group"
sample_1 = "pre"
sample_2 = "on"
sample_identity = "Time"

metadata <- data.frame(time = res.scLDAseq$settings$covariates$X[,"time"], 
                       sampleID = res.scLDAseq$sampleID)
prop <- data.frame(res.scLDAseq$theta )
colnames(prop) <- paste0("topic", 1:K)
metadata <- cbind(metadata, prop)
metadata <- data.frame(metadata)
## Get observed differences in fraction.
result <- metadata %>%
  group_by(time, sampleID) %>%
  summarise(
    sum_topic1 = sum(topic1),
    sum_topic2 = sum(topic2),
    sum_topic3 = sum(topic3)
  )
# Calculate the proportion for each topic
proportions <- metadata %>%
  group_by(time, sampleID) %>%
  summarise(
    topic1 = sum(topic1) / n(),
    topic2 = sum(topic2) / n(),
    topic3 = sum(topic3) / n()
  )

library(tidyr)
log2_diff <- proportions %>% gather(topic, topic_proportion, topic1:topic3, factor_key=TRUE)  %>%
  arrange(sampleID, topic, time) %>% # Ensure the data is ordered
  group_by(sampleID, topic) %>%
  summarise(
    log2_diff = log2(topic_proportion[time == 2] / topic_proportion[time == 1])
  ) %>%
  ungroup()

## Permutation test.
perm_res <- vector(mode = "list", length = nsample)
# perm_results <- matrix(NA, K, n_permutations)
# rownames(perm_results) <- sort(cluster_cases)
perm_res[[1]] <- matrix(, nrow = K, ncol = n_permutations)
perm_res[[2]] <- matrix(, nrow = K, ncol = n_permutations)

for (i in seq_len(n_permutations)) {
  
  metadata_permuted <- metadata %>%
    group_by(sampleID) %>%
    mutate(time_permuted = sample(time)) %>%
    ungroup()
  
  perm_log2_diff <- metadata_permuted %>%
    group_by(time_permuted, sampleID) %>%
    summarise(
      topic1 = sum(topic1) / n(),
      topic2 = sum(topic2) / n(),
      topic3 = sum(topic3) / n()
    ) %>% gather(topic, topic_proportion, topic1:topic3, factor_key=TRUE)  %>%
    arrange(sampleID, topic, time_permuted) %>% # Ensure the data is ordered
    group_by(sampleID, topic) %>%
    summarise(
      log2_diff = log2(topic_proportion[time_permuted == 2] / topic_proportion[time_permuted == 1])
    ) %>%
    ungroup()
  sample1 <- perm_log2_diff %>% filter(sampleID == "Batch1")
  perm_res[[1]][,i] <- sample1[["log2_diff"]]
  
  sample2 <- perm_log2_diff %>% filter(sampleID == "Batch2")
  perm_res[[2]][,i] <- sample2[["log2_diff"]]
}

for (s in 1:nsample) {
  sample_dat_log2_diff <- log2_diff %>% filter(sampleID == paste0("Batch", s))
  
  increased <- rowSums(apply(perm_res[[s]], 2, function(x) sample_dat_log2_diff[["log2_diff"]] <= x))
  increased <- (increased + 1) / (n_permutations + 1)
  
  decreased <- rowSums(apply(perm_res[[s]], 2, function(x) sample_dat_log2_diff[["log2_diff"]] >= x))
  decreased <- (decreased + 1) / (n_permutations + 1)
  
  sample_dat_log2_diff <- data.table(sample_dat_log2_diff)
  sample_dat_log2_diff[, pval := ifelse(log2_diff > 0, increased[.I], decreased[.I])]
  sample_dat_log2_diff[, FDR := p.adjust(pval, "fdr")]
}


## Boostrap log2FD CI.
boot_results <- matrix(NA, nrow(obs_diff), n_permutations)
rownames(boot_results) <- sort(cluster_cases)

########### linear mixed model ################
K <- 3
stmobj <- res.scLDAseq
metadata <- dat$meta
uncertainty = "Global"
formula = as.formula(1:K ~ time)
# library(parallel)
effM <- estimateMixedEffect(formula = as.formula(1:K ~ time), 
                           stmobj = res, meta= meta, uncertainty = "Global")
summary.estimateMixedEffect(effM)

response <- as.character(formula)[2] #second object is the response in this cases
K <- eval(parse(text=response))
formula <- formula(paste(as.character(formula)[c(1,3)], collapse = " "))
termobj <- terms(formula, data=metadata)

mf <- model.frame(termobj, data=metadata)
xmat <- model.matrix(termobj,data=metadata)
varlist <- all.vars(termobj)
metadata <- metadata[, varlist, drop=FALSE]
xmat <- as.data.frame(xmat)
xmat$sample <- stmobj$sampleID

thetasims <- thetaPosteriorSample(stmobj, nsims=1) # draw from posterior distribution
thetasims <- do.call(rbind, thetasims)

##### limma test ######
library(limma)
y.prop <- t(thetasims)
rownames(y.prop) <- paste0("topic",K)

Time <- factor(xmat$time)
levels(Time) <- c("pre", "on")
tdesign <- model.matrix(~0+Time)
colnames(tdesign) <- levels(Time)
corfit <- duplicateCorrelation(y.prop,tdesign,block=stmobj$sampleID)
fit <- lmFit(y.prop,design = tdesign,block=xmat$sampleID,correlation=corfit$consensus)
cm <- makeContrasts(Time1v2 = on-pre, levels=tdesign)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
topTable(fit2, coef="Time1v2")

# linear mixed model

thetaLogit <- log(thetasims/(1-thetasims))
output <- vector(mode="list", length=length(K))
for(k in K) {
  y <- thetaLogit[,k]
  # lm.mod = lmer(y ~ time + (1|sample), REML = F, data = xmat)
  
  form <- as.formula(paste0("y ~ ", as.character(formula)[-1], " + (1|sample)"))
  lm.mod = lmer(form, REML = F, data = xmat)
  # est <- summary(lm.mod)$coefficients[,1]
  lm_0 <-lmer(y ~ (1|sample), REML = F, data = xmat)
  anova(lm.mod,lm_0)
  PBmodcomp(lm.mod,lm_0,nsim=200)
  output[[which(k==K)]] <- est
}

# target <- data.frame(Subject =  rep(1:6, each =2), 
#                      +                      Condition = rep(c("Diseased", "Normal"), each = 6),
#                      +                      Tissue = rep(c("A", "B"), rep = 6))
# colnames(thetasims) <- pas
# compute 
# thetaLogit <- log(thetasims/(1-thetasims))
# output <- vector(mode="list", length=length(K))


##  
#Step 3: Calculate Coefficients
##

# first simulate theta 
storage <- vector(mode="list", length=nsims)

# setup parallel running if numCores is specified
if(numCores|is.numeric(numCores)) {
  if(is.numeric(numCores)) {
    cl <- numCores 
    if(cl > detectCores()) stop(paste0("Number of cores is more than available. The max cores available is ", detectCores()))
  } else{
    numCores <- detectCores()
  }
  cl <- makeCluster(numCores)
  clusterEvalQ(cl, {
    library(lme4)
  })
  clusterExport(cl, c("stmobj","xmat","K", "mixed.lm", "thetaPosteriorSample", "rmvnorm",
                      "row.lse", "formula"))
  # start_time <- Sys.time()
  storage <- parLapply(cl, 1:nsims,function(i) mixed.lm(stmobj, xmat, K, formula))
  stopCluster(cl)
  # end_time <- Sys.time()
  # elapsed_time <- end_time - start_time
} else{
  
  for(i in 1:nsims) {
    storage[[i]] <- mixed.lm(stmobj, xmat, K)
    if (i %% 10 == 0) {  # Print progress every 10 iterations
      cat(sprintf("...%d%% ", i /nsims * 100))
      flush.console()
    }
  }
  cat("\n")
}



# # using swish?
# # library(MASS)
# lambda <- res$eta
# nu <- res$nu
# nsims <- 10
# storage <- vector(mode="list", length=nsims)
# for (j in 1:nsims) {
#     out <- vector(mode="list",length=nrow(lambda)) 
#     for (i in 1:length(out)) {
#         sigma <- nu[[i]]
#         choleskydecomp <- chol(sigma)
#         mat <- rmvnorm(1, lambda[i,],nu[[i]],choleskydecomp)
#         mat <- cbind(mat, 0)
#         out[[i]] <- exp(mat - row.lse(mat))
#     }
#     out <- do.call(rbind, out)
#     rownames(out) <- paste0("cell", 1:nrow(lambda))
#     colnames(out) <- paste0("topic", 1:(ncol(lambda)+1))
#     storage[[j]] <- out
# }
# names(storage) <- paste0("infRep",1:nsims)
# library(fishpond)
# 
# infRepsArray <- abind::abind(storage,along=3)
