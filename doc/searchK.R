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
params <- setParams(params, group.prob = c(0.5,0.5),
                    de.prob = c(0.3, 0.3), 
                    nGenes = 5000, batchCells=c(1000,1500),  
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)# create change of proportions

# we assume that the cells pre and post treatment are equal
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.2,0.8))
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

#### feature selection #####
sims <- scuttle::logNormCounts(sims)
sub.feature <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(sub.feature, n=1000)
sims.sub <- sims[p2.chosen,]

##### eval 
r.file <- paste0("R/",list.files("R/"))
#r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
#sourceCpp("../stm/src/STMCfuns.cpp")
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims.sub)
K <- 2 
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
data <- dat$meta
sample <- "Batch"

storage <- searchK(documents, vocab, K = 2:8, prevalence =prevalence, data = data)
storage.high <- searchK(documents, vocab, K = 9:12, prevalence =prevalence, data = data)
x<-searchK(documents, vocab, K = 30, prevalence =prevalence, data = data)
plot(storage.high)
