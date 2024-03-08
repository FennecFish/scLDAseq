setwd("/proj/milovelab/wu/scLDAseq/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("No arguments supplied. Usage: Rscript test.R <filename>", call. = FALSE)
}

seed <- as.numeric(args[1])
batch <- as.logical(args[2])


cat(seed, "\n")
cat(batch, "\n")

# source("doc/eval.R")
# source("doc/methods.R")
#################################################
################ positive control ################
#################################################

params <- newSplatParams()
# params <- setParams(params, group.prob = c(0.45,0.45,0.1),
#                     de.prob = c(0.2, 0.2, 0.2), 
#                     nGenes = 5000, batchCells=c(2000,2000,2000), seed = seed)
params <- setParams(params, group.prob = c(0.2,0.3,0.5),
                    de.prob = c(0.1, 0.1, 0.05),
                    nGenes = 5000, batchCells=c(2000,2000,2000), seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = TRUE)

# we assume that the cells pre and post treatment are equal
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.1,0.8,0.1))
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
                    ngroup, "group_seed_", seed, ".rds")

saveRDS(sims.sub, file = file_name)

