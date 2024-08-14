setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_test4/")
dat <- data.frame()

all.Y <- vector(mode = "list")
for (file_name in files){
  scSTMobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_test4/", file_name))
  set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file_name)

  time <- scSTMobj$settings$covariates$X[,-1]
  names(time) <- scSTMobj$DocName

  theta <- scSTMobj$theta
  rownames(theta) <- scSTMobj$DocName

  t1 <- theta[match(names(time)[time==1], rownames(theta)),]
  t2 <- theta[match(names(time)[time==2], rownames(theta)),]

  names(scSTMobj$sampleID) <- scSTMobj$DocName
  sampleID <- unique(scSTMobj$sampleID)

  # the following code is to calculate the composition mean for each patient
  res <- data.frame()
  for (sample in sampleID){
    x1 <- t1[rownames(t1) %in% names(which(scSTMobj$sampleID==sample)),]
    x2 <- t2[rownames(t2) %in% names(which(scSTMobj$sampleID==sample)),]
    Y <- rbind(colMeans(x1),colMeans(x2))
    cat(file_name, "\n")
    res <- rbind(res, Y)
  }
  all.Y[[set_level]] <- res
}
saveRDS(all.Y, file = "res/allY_singleV3_prop.rds")


# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", pattern = "c3")
# all.Y <- vector(mode = "list")
# for (file_name in files){
#   sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", file_name))
#   set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
#   
#   dat <- colData(sims) %>%
#     as.data.frame() %>%
#     dplyr::count(Batch, Group, time) %>%
#     dplyr::group_by(Batch, time) %>%
#     dplyr::mutate(Proportion = n / sum(n)) %>%
#     dplyr::select(Batch, time, Proportion) %>%
#     dplyr::group_by(time, Batch) %>%
#     dplyr::mutate(id = row_number()) %>%
#     tidyr::pivot_wider(
#       names_from = id,
#       values_from = Proportion,
#       names_prefix = "K_")
#   all.Y[[set_level]] <- dat
#   cat(file_name, "\n")
# }
# saveRDS(all.Y, file = "allY_V5_c3_sims.rds")
# 
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", pattern = "c9")
# all.Y <- vector(mode = "list")
# for (file_name in files){
#   sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V5/sims/", file_name))
#   set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
#   
#   dat <- colData(sims) %>%
#     as.data.frame() %>%
#     dplyr::count(Batch, Group, time) %>%
#     dplyr::group_by(Batch, time) %>%
#     dplyr::mutate(Proportion = n / sum(n)) %>%
#     dplyr::select(Batch, time, Proportion) %>%
#     dplyr::group_by(time, Batch) %>%
#     dplyr::mutate(id = row_number()) %>%
#     tidyr::pivot_wider(
#       names_from = id,
#       values_from = Proportion,
#       names_prefix = "K_")
#   all.Y[[set_level]] <- dat
#   cat(file_name, "\n")
# }
# saveRDS(all.Y, file = "allY_V5_c9_sims.rds")
# 
# 
# # #### construct all.Y for scSTM V3 #####
# # level <- paste0("L", 4:9)
# # 
# # all.Y <- vector(mode = "list")
# # for (l in level){
# #   files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/scSTM/", pattern = l)
# #   res <- data.frame()
# #   for (file_name in files){
# #     scSTMobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/scSTM/", file_name))
# #     theta <- scSTMobj$theta
# #     rownames(theta) <- scSTMobj$DocName
# #     time <- scSTMobj$settings$covariates$X[,2]
# #     names(time) <- scSTMobj$DocName
# #     t1 <- theta[match(names(time)[time==1], rownames(theta)),]
# #     t2 <- theta[match(names(time)[time==2], rownames(theta)),]
# #     Y <- rbind(colMeans(t1),colMeans(t2))
# #     cat(file_name, "\n")
# #     res <- rbind(res, Y)
# #   }
# #   all.Y[[l]] <- res
# # }
# # saveRDS(all.Y, file = "allY_V1_scSTM.rds")
# 
# # #### construct all.Y for scSTM V3 #####
# # control <- c("neg", "pos")
# # nCellType <- paste0("c",c(3,5,8))
# # level <- expand.grid(control, nCellType) %>% mutate(level = paste0(Var1,"_",Var2)) %>% select(level)
# # level <- level$level
# # 
# # all.Y <- vector(mode = "list")
# # for (l in level){
# #   files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3/scSTM/", pattern = l)
# #   res <- data.frame()
# #   for (file_name in files){
# #     scSTMobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3/scSTM/", file_name))
# #     theta <- scSTMobj$theta
# #     rownames(theta) <- scSTMobj$DocName
# #     time <- scSTMobj$settings$covariates$X[,2]
# #     names(time) <- scSTMobj$DocName
# #     t1 <- theta[match(names(time)[time==1], rownames(theta)),]
# #     t2 <- theta[match(names(time)[time==2], rownames(theta)),]
# #     Y <- rbind(colMeans(t1),colMeans(t2))
# #     cat(file_name, "\n")
# #     res <- rbind(res, Y)
# #   }
# #   all.Y[[l]] <- res
# # }
# # saveRDS(all.Y, file = "all_Y_scSTM.rds")
# 
# ######### construct dataset for sims V1 ###########
# # control <- c("neg", "pos")
# # nCellType <- paste0("c",c(3,5,8))
# # level <- expand.grid(control, nCellType) %>% mutate(level = paste0(Var1,"_",Var2)) %>% select(level)
# # level <- level$level
# # level <- paste0("L", 4:9)
# # 
# # all.Y <- vector(mode = "list")
# # for (l in level){
# #   files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", pattern = l)
# #   res <- data.frame()
# #   for (file_name in files){
# #     sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", file_name))
# #     dat <- colData(sims) %>%
# #       as.data.frame() %>%
# #       dplyr::count(Batch, Group, time) %>%
# #       dplyr::group_by(time) %>%  
# #       dplyr::mutate(Proportion = n / sum(n)) %>%
# #       dplyr::select(Batch, time, Proportion) %>%
# #       dplyr::group_by(time, Batch) %>%
# #       dplyr::mutate(id = row_number()) %>%
# #       tidyr::pivot_wider(
# #         names_from = id,
# #         values_from = Proportion,
# #         names_prefix = "K_")
# #     
# #     cat(file_name, "\n")
# #     res <- rbind(res, dat)
# #   }
# #   all.Y[[l]] <- res
# # }
# # saveRDS(all.Y, file = "allY_V1_sims.rds")
# 
# 
# ######## construct the dataset for sims V3 ########
# # control <- c("neg", "pos")
# # nCellType <- paste0("c",c(3,5,8))
# # level <- expand.grid(control, nCellType) %>% mutate(level = paste0(Var1,"_",Var2)) %>% select(level)
# # level <- level$level
# # 
# # all.Y <- vector(mode = "list")
# # for (l in level){
# #   files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", pattern = l)
# #   res <- data.frame()
# #   for (file_name in files){
# #     sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", file_name))
# #     dat <- colData(sims) %>%
# #       as.data.frame() %>%
# #       dplyr::count(Batch, Group, time) %>%
# #       dplyr::group_by(time) %>%  
# #       dplyr::mutate(Proportion = n / sum(n)) %>%
# #       dplyr::select(time, Proportion) %>%
# #       dplyr::group_by(time) %>%
# #       dplyr::mutate(id = row_number()) %>%
# #       tidyr::pivot_wider(
# #         names_from = id,
# #         values_from = Proportion,
# #         names_prefix = "K_")
# #     
# #     cat(file_name, "\n")
# #     res <- rbind(res, dat)
# #   }
# #   all.Y[[l]] <- res
# # }
# # saveRDS(all.Y, file = "all_Y.rds")
# # 
