setwd("/proj/milovelab/wu/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(mclust)
library(Seurat)
library(cluster)
library(monocle3)
library(sctransform)

process_scSTM <- function(scSTMobj) {
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- colnames(scSTMobj$mu$mu)
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", pattern = "sims*")
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", pattern = "scSTM*")
dat <- data.frame()
level <- paste0("L", 1:9)
res.adj <- data.frame()

for(l in level){
  scSTM_file_name <- grep(paste0(l,".rds$"), files, value = TRUE)
  sim_name <- unique(sub(".*_([0-9]+_L[0-9]+).*", "\\1", scSTM_file_name))
  res <- data.frame()
  
  for(file in sim_name){
    set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file)
    # set_level <- file
    sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/sims_", file, ".rds"))
    dat <- colData(sims) %>% data.frame() %>% select(Cell:Group,time)
    
    # scSTM
    # scSTM_path <- grep(file, scSTM_file_name, value = TRUE)
    # for(path in scSTM_path){
    #   scSTM <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", path))
    #   scSTM_cluster <- process_scSTM(scSTM)
    #   pattern <- "(_[0-9]+_L[0-9]+\\.rds)$"
    #   scSTM_setup <- sub(pattern, "", path)
    # 
    #   dat[,scSTM_setup] <- scSTM_cluster[match(dat$Cell, names(scSTM_cluster))]
    #   rm(scSTM)
    # }
    
    # fileter gene with no content model
    f_nc_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_f_nc/scSTM_filtergenes_noContent_",set_level,".rds")
    # scSTM_f_nc <- readRDS(f_nc_name)
    # scSTM_f_nc_cluster <- process_scSTM(scSTM_f_nc)
    # dat$scSTM_f_nc_cluster <- scSTM_f_nc_cluster[match(dat$Cell, names(scSTM_f_nc_cluster))]
    if(file.exists(f_nc_name)){
      scSTM_f_nc <- readRDS(f_nc_name)
      scSTM_f_nc_cluster <- process_scSTM(scSTM_f_nc)
      dat$scSTM_f_nc_cluster <- scSTM_f_nc_cluster[match(dat$Cell, names(scSTM_f_nc_cluster))]
    } else {dat$scSTM_f_nc_cluster = NA}

    
    # all gene with no content model
    a_nc_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_a_nc/scSTM_allGenes_noContent_",set_level,".rds")
    if(file.exists(a_nc_name)){
      scSTM_a_nc <- readRDS(a_nc_name)
      scSTM_a_nc_cluster <- process_scSTM(scSTM_a_nc)
      dat$scSTM_a_nc_cluster <- scSTM_a_nc_cluster[match(dat$Cell, names(scSTM_a_nc_cluster))]
    }else {dat$scSTM_a_nc_cluster = NA}
    
    # filter genes with content model
    f_c_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_f_c/scSTM_filterGenes_Content_",set_level,".rds")
    if(file.exists(f_c_name)){
      scSTM_f_c <- readRDS(f_c_name)
      scSTM_f_c_cluster <- process_scSTM(scSTM_f_c)
      dat$scSTM_f_c_cluster <- scSTM_f_c_cluster[match(dat$Cell, names(scSTM_f_c_cluster))]
    } else {dat$scSTM_f_c_cluster = NA}
    
    # all genes wtih content model
    a_c_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_a_c/scSTM_allGenes_Content_",set_level,".rds")
    if(file.exists(a_c_name)){
      scSTM_a_c <- readRDS(a_c_name)
      scSTM_a_c_cluster <- process_scSTM(scSTM_a_c)
      dat$scSTM_a_c_cluster <- scSTM_a_c_cluster[match(dat$Cell, names(scSTM_a_c_cluster))]
    } else {dat$scSTM_a_c_cluster = NA}
    
    # filter gene with combat no content 
    if (l %in% c("L1","L2","L3")){
      combat_f_nc_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_f_nc/scSTM_filtergenes_noContent_",set_level,".rds")
    } else{
      combat_f_nc_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_combat_f_nc/scSTM_combat_filtergenes_noContent_",set_level,".rds")
    }
    if(file.exists(combat_f_nc_name)){
      scSTM_combat_f_nc <- readRDS(combat_f_nc_name)
      scSTM_combat_f_nc_cluster <- process_scSTM(scSTM_combat_f_nc)
      dat$scSTM_combat_f_nc_cluster <- scSTM_combat_f_nc_cluster[match(dat$Cell, names(scSTM_combat_f_nc_cluster))]
    } else {dat$scSTM_combat_f_nc_cluster = NA}
    
    # filter gene with content with combat
    if (l %in% c("L1","L2","L3")){
      combat_f_c_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_f_c/scSTM_filterGenes_Content_",set_level,".rds")
    } else{
      combat_f_c_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_combat_f_c/scSTM_combat_filterGenes_Content_",set_level,".rds")
    }
    if(file.exists(combat_f_c_name)){
      scSTM_combat_f_c <- readRDS(combat_f_c_name)
      scSTM_combat_f_c_cluster <- process_scSTM(scSTM_combat_f_c)
      dat$scSTM_combat_f_c_cluster <- scSTM_combat_f_c_cluster[match(dat$Cell, names(scSTM_combat_f_c_cluster))]
    } else {dat$scSTM_combat_f_c_cluster = NA}

    # Seurat
    seurat.sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/seurat/seurat_", set_level, ".rds"))
    smeta <- seurat.sims@meta.data %>% as.data.frame()
    sub_sims <- sims[,rownames(smeta)] # filter by the rows
    seurat.adj <- sapply(smeta[,4:7], function(x) {
      adjustedRandIndex(x, sub_sims$Group)
    })
    # select the resolution that has the highest ARI. When there are multiple, select the first one
    best_res <- names(seurat.adj)[seurat.adj == max(seurat.adj)][1]
    seurat_cluster <- seurat.sims@meta.data %>% as.data.frame() %>% select(all_of(best_res))
    dat$seurat_cluster <- seurat_cluster[match(dat$Cell, rownames(seurat_cluster)),]
    rm(sub_sims)

    # sctransform
    sctf <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/sctransform_multi/sctransform_", set_level, ".rds"))
    sct_meta <- sctf@meta.data %>% as.data.frame()
    sub_sims <- sims[,rownames(sct_meta)] # filter by the rows
    sctf.adj <- sapply(sct_meta[,6:9], function(x) {
      adjustedRandIndex(x, sub_sims$Group)
    })
    # select the resolution that has the highest ARI. When there are multiple, select the first one
    best_res <- names(sctf.adj)[sctf.adj == max(sctf.adj)][1]
    sctf_cluster <- sctf@meta.data %>% as.data.frame() %>% select(all_of(best_res))
    dat$sctransform_cluster <- sctf_cluster[match(dat$Cell, rownames(sctf_cluster)),]
    rm(sub_sims)
    
    # fastTopics
    fasttopic_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/fastTopics/fastTopics_", set_level, ".rds")
    nmf.sims <- readRDS(fasttopic_name)
    max_indices <- apply(nmf.sims$L, 1, which.max)
    # colnames(nmf.sims$L) <- paste0("topic_", 1:ncol(nmf.sims$L))
    # rownames(nmf.sims$L) <- colnames(scSTMobj$mu$mu)
    fastTopics_cluster <- colnames(nmf.sims$L)[max_indices]
    names(fastTopics_cluster) <- rownames(nmf.sims$L)
    dat$fastTopics_cluster <- fastTopics_cluster[match(dat$Cell, names(fastTopics_cluster))]
    # if(file.exists(fasttopic_name)){
    #   nmf.sims <- readRDS(fasttopic_name)
    #   max_indices <- apply(nmf.sims$L, 1, which.max)
    #   # colnames(nmf.sims$L) <- paste0("topic_", 1:ncol(nmf.sims$L))
    #   # rownames(nmf.sims$L) <- colnames(scSTMobj$mu$mu)
    #   fastTopics_cluster <- colnames(nmf.sims$L)[max_indices]
    #   names(fastTopics_cluster) <- rownames(nmf.sims$L)
    #   dat$fastTopics_cluster <- fastTopics_cluster[match(dat$Cell, names(fastTopics_cluster))]
    # } else {dat$fastTopics_cluster = NA}
    
    # monocle3
    monocle3 <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/monocle3/monocle3_", set_level, ".rds"))
    dat$monocle3_cluster <- partitions(monocle3)[match(dat$Cell, names(partitions(monocle3)))]
    
    # SC3
    sc3 <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/sc3/sc3_", set_level, ".rds"))
    sc3_cluster <- sc3@colData %>% as.data.frame() %>% select(ends_with("clusters"))
    dat$sc3_cluster <- sc3_cluster[match(dat$Cell, rownames(sc3_cluster)),]
      
    adjusted_rand_indices <- sapply(dat[, 5:ncol(dat)], function(x) {
      adjustedRandIndex(x, sims$Group)
    })
    
    # Create a data frame to store results
    res.temp <- data.frame(
      sim = set_level,
      level = l,
      t(adjusted_rand_indices)  # Transpose to match the original data frame structure
    )
    # res.temp <- data.frame(
    #   sim = set_level,
    #   level = l,
    #   scSTM_a_nc_adjR = adjustedRandIndex(dat$scSTM_a_nc_cluster,sims$Group),
    #   scSTM_f_nc_adjR = adjustedRandIndex(dat$scSTM_f_nc_cluster,sims$Group),
    #   scSTM_f_c_adjR = adjustedRandIndex(dat$scSTM_f_c_cluster,sims$Group),
    #   Seurat_adjR = adjustedRandIndex(dat$seurat_cluster, sims$Group), 
    #   fastTopics_adjR = adjustedRandIndex(dat$fastTopics_cluster,sims$Group))

    res <- bind_rows(res, res.temp)
    rm(sims)
  }
  
  res.adj <- bind_rows(res.adj, res)
}

write.csv(res.adj, file = "res/res_fig1_adj.csv")



# read in scSTM
# 
# 
# msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
# if(verbose) cat(msg)
# ###########################################################
# #################### Seurat ###############################
# ###########################################################
# t1 <- proc.time()
# 
# seurat.sims <- as.Seurat(sims, counts = "counts", data = "logcounts")
# seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 500)
# all.genes <- rownames(seurat.sims)
# seurat.sims <- ScaleData(seurat.sims, features = all.genes)
# seurat.sims <- RunPCA(seurat.sims)
# 
# seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
# seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)
# 
# 
# msg <- sprintf("Completed Seurat (%d seconds). \n", floor((proc.time()-t1)[3]))
# if(verbose) cat(msg)
# saveRDS(seurat.sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/seurat_", sim_name, ".rds"))
# ###########################################################
# ###################### CIDR ###############################
# ###########################################################
# # t1 <- proc.time()
# # 
# # res.cidr <- scDataConstructor(as.matrix(counts(sims)))
# # res.cidr <- determineDropoutCandidates(res.cidr)
# # res.cidr <- wThreshold(res.cidr)
# # res.cidr <- scDissim(res.cidr)
# # res.cidr <- scPCA(res.cidr)
# # res.cidr <- nPC(res.cidr)
# # nCluster(res.cidr)
# # res.cidr <- scCluster(res.cidr)
# # dat$cidr_cluster <- res.cidr@clusters[match(colnames(res.cidr@tags), dat$Cell)]
# # 
# # msg <- sprintf("Completed CIDR (%d seconds). \n", floor((proc.time()-t1)[3]))
# # if(verbose) cat(msg)
# ###########################################################
# ##################### RACEID ###############################
# ###########################################################
# t1 <- proc.time()
# 
# # tutorial
# # https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
# sc <- SCseq(counts(sims))
# # Cells with a relatively low total number of transcripts are discarded.
# sc <- filterdata(sc,mintotal=2000) 
# # retrieve filtered and normalized expression matrix 
# # (normalized to the minimum total transcript count across all cells retained after filtering) 
# fdata <- getfdata(sc)
# # If all genes should be used, then the parameter FSelect needs to be set to FALSE. 
# sc <- compdist(sc, metric="pearson", FSelect = FALSE)
# # sc <- clustexp(sc)
# sc <- clustexp(sc,cln=ngroup,sat=FALSE) # FUNcluster for other methods
# sc <- findoutliers(sc)
# dat$raceID_cluster <- NA
# raceID_cluster <- sc@cluster$kpart[match(names(sc@cluster$kpart), dat$Cell)]
# dat$raceID_cluster <- sc@cluster$kpart[match(dat$Cell,names(sc@cluster$kpart))]
# msg <- sprintf("Completed RaceID (%d seconds). \n", floor((proc.time()-t1)[3]))
# if(verbose) cat(msg)
# 
# saveRDS(sc, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/raceID_", sim_name, ".rds"))
# 
# set.seed(1)
# 
# #### evaluating all methods #####
# sc_eval <- function(sims, dat) {
#   
#   res <- data.frame()
#   # compute silhouette score
#   # dist.matrix <- dist(t(counts(sims)))
#   # scSTM.sil <- silhouette(as.numeric(as.factor(dat$scSTM_cluster)), dist.matrix)
#   # seurat.sil <- silhouette(as.numeric(as.factor(dat$seurat_cluster)), dist.matrix)
#   # raceid.sil <- silhouette(as.numeric(as.factor(dat$raceID_cluster)), dist.matrix)
#   # cidr.sil <- silhouette(as.numeric(as.factor(dat$cidr_cluster)), dist.matrix)
#   
#   res <- data.frame(
#     scSTM_adjR = adjustedRandIndex(dat$scSTM_cluster,sims$Group),
#     Seurat_adjR = adjustedRandIndex(dat$seurat_cluster, sims$Group), 
#     raceID_adjR = adjustedRandIndex(dat$raceID_cluster,sims$Group),#,
#     CIDR_adjR = adjustedRandIndex(dat$cidr_cluster,sims$Group))#,
#   # scSTM_sil = mean(scSTM.sil[,3]),
#   # seurat_sil = mean(seurat.sil[,3]),
#   # raceID_sil = mean(raceid.sil[,3]))
#   
#   return(res)
# }
# 
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/", pattern = "3cellTypes_sims.rds")
# sim_name <- sub("\\_sims.rds$", "", files)
# res <- data.frame()
# 
# for (i in sim_name){
#   dat <- read.csv(paste0("res/clustering_benchmark/colData_",i,".csv"))
#   sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",i,"_sims.rds"))
#   temp <- sc_eval(sims, dat = dat)
#   res <- rbind(res, temp)
#   rm(dat)
#   rm(sims)
#   cat(i, "\n")
# }
# 
# rownames(res) <- sim_name
# res <- tibble::rownames_to_column(res, "sim")
# res$NumCellType <- sub(".*samples_(\\d+)cellTypes", "\\1", res$sim)
# res$NumSample <- sub(".*_(\\d+samples).*", "\\1", res$sim)
# write.csv(res, file = "res/res_eval_clustering_4methods.csv")
# 
# # change the format from wide to long
# library(tidyr)
# res <- read.csv("res/res_eval_methods.csv")
# wcidr <- read.csv("res/res_eval_clustering_4methods.csv")
# 
# wcidr_long <- wcidr %>% gather(method, adjR, scSTM_adjR:CIDR_adjR, factor_key=TRUE)
# # gather(res, method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
# ggplot(wcidr_long, aes(x=method, y=adjR)) + 
#   geom_boxplot() + 
#   ggtitle("Comparison of Clustering Accuracy among Methods with 3 cell Types")
# 
# res <- res %>%
#   mutate(num_samples = as.numeric(gsub(".*_(\\d+)samples.*", "\\1", sim)),
#          num_cell_types = as.numeric(gsub(".*_(\\d+)cellTypes", "\\1", sim)))
# data_long <- res %>% filter(num_cell_types == 5) %>% gather(method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
# # gather(res, method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
# ggplot(data_long, aes(x=method, y=adjR)) + 
#   geom_boxplot() + 
#   ggtitle("Comparison of Clustering Accuracy among Methods with 5 cell Types")
# 
# 
# # looking at why scSTM is low
# 
# sim_name <- "BIOKEY_11_10samples_5cellTypes"
# sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",sim_name,"_sims.rds"))
