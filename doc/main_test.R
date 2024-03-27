setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
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
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

params <- newSplatParams()
params <- setParams(params, group.prob = c(0.5,0.2, 0.3),
                    de.prob = c(0.3, 0.3, 0.3), 
                    nGenes = 200, batchCells=c(100,150),  
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)# create change of proportions

cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.2,0.5, 0.3))
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
prop.table(table(sims$time,sims$Group), margin = 1)



dat <- prepsce(sims)
K <- length(unique(sims$Group))
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
meta <- dat$meta
sample <- "Batch"

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

res <- multi_stm(documents = documents, vocab = vocab,
                       K = K, prevalence = prevalence, content = NULL,
                       data = meta, 
                       sce = sce,
                       sample = sample,
                       init.type= "Spectral",
                       gamma.prior= "Pooled",
                       kappa.prior= "L1",
                       control = list(gamma.maxits=3000))

K <- 3
stmobj <- res
metadata <- meta
effM <- estimateMixedEffect(alt.formula = 1:K ~ time, 
                           null.formula = NULL,
                           stmobj = res, meta= meta, uncertainty = "Global")
plot.fixedEffect(effM, "time")

eff <- estimateEffect(1:K ~ time, 
                            stmobj = res, meta= meta, uncertainty = "Global")
summary(multi_eff)

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
