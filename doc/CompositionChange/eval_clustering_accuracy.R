setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)


process_scSTM <- function(scSTMobj) {
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- colnames(scSTMobj$mu$mu)
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/sims/", pattern = "sims*")
dat <- data.frame()

# # control <- c("neg_L1", paste0("pos_L", 1:4))
# control <- c("neg_L1")
# # nCellType <- paste0("c",c(5,9))
# nCellType <- "c5"
# level <- expand.grid(control, nCellType) %>% mutate(level = paste0(Var1,"_",Var2)) %>% select(level)
# level <- level$level
level <- "neg_L1_c5"

calc_cluster <- function(stmFiles, l){
  res <- data.frame()
  for(file_name in stmFiles){
    set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file_name)
    # seed <- sub("(.*)_.*$", "\\1", set_level)
    sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/sims/sims_", set_level, ".rds"))
    dat <- colData(sims) %>% data.frame() %>% select(Cell:Group,time)
    
    scSTM_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/test_fastTopics/",file_name)
    scSTMobj <- readRDS(scSTM_name)
    cluster <- process_scSTM(scSTMobj)
    dat$cluster <- cluster[match(dat$Cell, names(cluster))]
    
    adjusted_rand_indices <- adjustedRandIndex(dat$Group, dat$cluster)
    
    # Create a data frame to store results
    res.temp <- data.frame(
      sim = set_level,
      level = l,
      t(adjusted_rand_indices)  # Transpose to match the original data frame structure
    )
    
    res <- bind_rows(res, res.temp)
    cat(stmFiles, "\n")
  }
  return(res)
}

res.adj <- data.frame()
for(l in level){
  stmFiles <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V2/sims/",
                         pattern = l)
  res <- calc_cluster(stmFiles, l)
  res.adj <- bind_rows(res.adj, res)
}

write.csv(res.adj, file = "res/res_CompositionChange_adjustedRandIndex_V2_train_dat.csv")

# acc <- read.csv("res/res_CompositionChange_adjustedRandIndex_V2_train_dat.csv")
# colnames(acc) <- c("row_num", "sim", "level", "adjR")
# 
# ggplot(acc, aes(x=row_num, y=adjR)) +
#   geom_point() +
#   labs(title = "Adjusted Rand Index on Train Data",
#        x = "replicate",
#        y = "adjustedRandIndex")
# 
# ggplot(acc, aes(x = level, y = adjR)) +
#   geom_boxplot() +
#   labs(title = "Boxplot of Adjusted Rand Indices by Level",
#        x = "Level",
#        y = "Adjusted Rand Indices") +
#   theme_minimal()
