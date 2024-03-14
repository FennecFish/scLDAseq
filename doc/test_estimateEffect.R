# setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
setwd("/proj/milovelab/wu/scLDAseq/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(mclust)

##### estimating effect size for real data
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

res.multi <- readRDS("res/Anti-PD1_mix_K40_stmRes.rds")
K = res.multi$settings$dim$K
sims <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_sub_mix_response_samples_1000g.rds")
meta <- colData(sims) %>% data.frame()

# multi_eff <- estimateEffect(1:K ~ timepoint,
#                             stmobj = res.multi, 
#                             meta= meta, uncertainty = "Global")


multi_eff_26 <- estimateEffect(36:39 ~ timepoint, 
                            stmobj = res.multi, 
                            sampleNames = "patient_id", 
                            sampleIDs = "BIOKEY_26",
                            meta= meta, uncertainty = "Global")

summary(multi_eff_26)
plot(multi_eff_26, "timepoint", xlim = c(-0.3,0.3),
     main = "Estimated Effect of Time as a Covariate to Cell-Topic Proportion for Sample 2",
     method="difference",cov.value1="Pre",cov.value2="On",
     model = res.multi,
     labeltype = "frex",
     xlab = "Proportion Difference (95% Confidence Interval) ")

# plot()
##################################################
##### estimating effectsize for simulation #####
##################################################
sims <- readRDS("data/PositiveControl_2sample_3group_flip.rds")

time_prop <- colData(sims) %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))

# run multiSTM
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims)
K <- length(unique(sims$Group))
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
meta <- dat$meta
sample <- "Batch"

res.multi <- multi_stm(documents = documents, vocab = vocab,
                 K = K, prevalence = prevalence, content = NULL,
                 data = meta, 
                 sce = sce,
                 sample = sample,
                 init.type= "Spectral",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 control = list(gamma.maxits=3000))

multi_eff <- estimateEffect(1:K ~ time, stmobj = res.multi, meta= meta, uncertainty = "Global")
summary(multi_eff)

multi_eff_p1 <- estimateEffect(1:K ~ time, stmobj = res.multi, meta= meta, 
                               sampleNames = "Batch", sampleIDs = "Batch1",uncertainty = "Global")
summary(multi_eff_p1)
plot(multi_eff_p1, "time", model=res.multi,
     method="difference",cov.value1=1,cov.value2=2)
png("res/estimate_effect_TwoGroupChange_patient1.png", width = 3500, height = 2500, res = 300)
plot(multi_eff_p1, "time", xlim = c(-0.3,0.3),
     main = "Estimated Effect of Time as a Covariate to Cell-Topic Proportion for Sample 1",
     method="difference",cov.value1=1,cov.value2=2,
     model = res.multi,
     labeltype = "frex", 
     xlab = "Proportion Difference (95% Confidence Interval) ")
dev.off()

multi_eff_p2 <- estimateEffect(1:K ~ time, stmobj = res.multi, meta= meta, 
                               sampleNames = "Batch", sampleIDs = "Batch2",uncertainty = "Global")
summary(multi_eff_p2)
png("res/estimate_effect_TwoGroupChange_patient2.png", width = 3500, height = 2500, res = 300)
plot(multi_eff_p2, "time", xlim = c(-0.3,0.3),
     main = "Estimated Effect of Time as a Covariate to Cell-Topic Proportion for Sample 2",
     method="difference",cov.value1=1,cov.value2=2,
     model = res.multi,
     labeltype = "frex",
     xlab = "Proportion Difference (95% Confidence Interval) ")
dev.off()


# plot the true proportion
time_prop <- colData(sims) %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n)) %>%
    mutate(time = ifelse(time==1, "pre-treatment", "post-treatment")) %>%
    mutate(Batch = gsub("Batch", "Sample", Batch)) %>%
    mutate(Group = gsub("Group", "CellType_", Group))
time_prop$time <- factor(time_prop$time, levels=c("pre-treatment", "post-treatment"))
png("res/true_prop.png", height = 2500, width = 2500, res = 300)
ggplot(time_prop, aes(x = Group, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type True Proportion Change") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") + 
    facet_grid(~Batch) + 
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=0.2)) +
    xlab(NULL)
dev.off()


# run starndard stm
r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
sourceCpp("../stm/src/STMCfuns.cpp")

stm.res <- stm(documents = documents, vocab = vocab,
                       K = K, prevalence = prevalence, content = NULL,
                       data = meta, 
                       sce = sce,
                       init.type= "Spectral",
                       gamma.prior= "Pooled",
                       kappa.prior= "L1",
                       control = list(gamma.maxits=3000))




stm_eff <-estimateEffect(1:K ~ time, stmobj = stm.res, meta= meta,uncertainty = "Global")
summary(stm_eff)

#################################################
################ positive control ################
#################################################

params <- newSplatParams()
params <- setParams(params, group.prob = c(0.1,0.45,0.45),
                    de.prob = c(0.2, 0.2, 0.2),
                    nGenes = 500, batchCells=c(1000,1000))
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# we assume that the cells pre and post treatment are equal
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.1,0.65,0.25))
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

#### QC ######
sims <- quickPerCellQC(sims)
#### feature selection #####
sims <- scuttle::logNormCounts(sims)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=2000)
sims.sub <- sims[p2.chosen,]

nsample <- length(unique(sims.sub$Batch))
ngroup <- length(unique(sims.sub$Group))

file_name <- paste0("data/PositiveControl_", nsample, "sample_",
                    ngroup, "group.rds")

saveRDS(sims.sub, file = file_name)


# swap patient prop change
rand.samp <- sample(unique(sims.sub$Batch), size = 1, replace = FALSE)
sampled_data <- colData(sims.sub) %>%
    data.frame() %>%
    mutate(new_time = ifelse(time == 1 & Batch %in% rand.samp, 2, 1)) %>%
    mutate(time = ifelse(Batch %in% rand.samp, new_time, time))
sims.sub$time <- sampled_data$time
file_name <- paste0("data/PositiveControl_", nsample, "sample_",
                    ngroup, "group_flip.rds")
saveRDS(sims.sub, file = file_name)
