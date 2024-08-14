setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)
library(fastTopics)

# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/sims/", pattern = "sims*")
res.adj <- data.frame()
level <- "neg_L1_c5"

for(l in level){
  files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/scSTM/",
                         pattern = l)[-1]
  res <- data.frame()
  for(file in files){
    set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file)
  
    sim_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/sims/sims_", set_level, ".rds")
    sims <- readRDS(sim_name)
    dat <- colData(sims) %>% data.frame() %>% dplyr::select(Cell:Group,time)
    
    fasttopic_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/test_fastTopics/fastTopics_", set_level, ".rds")
    ft_res <- readRDS(fasttopic_name)
    max_indices <- apply(ft_res$L, 1, which.max)
    # colnames(nmf.sims$L) <- paste0("topic_", 1:ncol(nmf.sims$L))
    # rownames(nmf.sims$L) <- colnames(scSTMobj$mu$mu)
    fastTopics_cluster <- colnames(ft_res$L)[max_indices]
    names(fastTopics_cluster) <- rownames(ft_res$L)
    dat$cluster <- fastTopics_cluster[match(dat$Cell, names(fastTopics_cluster))]
    
    adjusted_rand_indices <- adjustedRandIndex(dat$Group, dat$cluster)
    
    # Create a data frame to store results
    res.temp <- data.frame(
      sim = set_level,
      level = l,
      t(adjusted_rand_indices)  # Transpose to match the original data frame structure
    )
    
    res <- bind_rows(res, res.temp)
    cat(file, "\n")
  }
  res.adj <- bind_rows(res.adj, res)
}


write.csv(res.adj, file = "res/res_CompositionChange_adjustedRandIndex_V2_test_dat.csv")

acc <- read.csv("res/res_CompositionChange_adjustedRandIndex_V2_test_dat.csv")
colnames(acc) <- c("row_num", "sim", "level", "adjR")

ggplot(acc, aes(x=row_num, y=adjR)) +
  geom_point() +
  labs(title = "Adjusted Rand Index on Train Data",
       x = "replicate",
       y = "adjustedRandIndex")


ggplot(acc, aes(x = level, y = adjR)) +
  geom_boxplot() +
  labs(title = "Boxplot of Adjusted Rand Indices by Level",
       x = "Level",
       y = "Adjusted Rand Indices") +
  theme_minimal()
