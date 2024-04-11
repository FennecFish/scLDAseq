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
library(RaceID)
library(cidr)
library(cluster)

sc_methods <- function(sims, sim_name, verbose = TRUE) {
    
    ngroup <- length(unique(sims$Group))
    dat <- colData(sims) %>% 
        data.frame() 

    ###########################################################
    ################## scLDAseq ###############################
    ###########################################################
    r.file <- paste0("R/",list.files("R/"))
    sapply(r.file, source)
    sourceCpp("src/STMCfuns.cpp")

    t1 <- proc.time()

    K <- ngroup

    stm_dat <- prepsce(sims)
    prevalence <- as.formula(~stm_dat$meta$time)
    content <- NULL
    sce <- stm_dat$sce
    documents  <- stm_dat$documents
    vocab <- stm_dat$vocab
    data <- stm_dat$meta
    sample <- "Batch"

    res.stm <- multi_stm(documents = documents, vocab = vocab,
                         K = K, prevalence = prevalence, content = NULL,
                         data = data,
                         sce = sce,
                         sample = sample,
                         init.type= "Spectral",
                         gamma.prior= "Pooled",
                         kappa.prior= "L1",
                         control = list(gamma.maxits=3000))
    saveRDS(res.stm, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/scLDAseq_", sim_name, ".rds"))
    max_indices <- apply(res.stm$theta, 1, which.max)
    colnames(res.stm$theta) <- paste0("topic_", 1:ncol(res.stm$theta))
    rownames(res.stm$theta) <- colnames(res.stm$mu$mu)
    res_cluster <- colnames(res.stm$theta)[max_indices]
    names(res_cluster) <- rownames(res.stm$theta)
    dat$scSTM_cluster <- res_cluster[match(names(res_cluster), dat$Cell)]

    msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
    if(verbose) cat(msg)
    ###########################################################
    #################### Seurat ###############################
    ###########################################################
    t1 <- proc.time()
    
    seurat.sims <- as.Seurat(sims, counts = "counts", data = "logcounts")
    seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 500)
    all.genes <- rownames(seurat.sims)
    seurat.sims <- ScaleData(seurat.sims, features = all.genes)
    seurat.sims <- RunPCA(seurat.sims)
    
    seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
    seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)
    
    dat$seurat_cluster <- Idents(seurat.sims)[match(names(Idents(seurat.sims)), dat$Cell)]
    msg <- sprintf("Completed Seurat (%d seconds). \n", floor((proc.time()-t1)[3]))
    if(verbose) cat(msg)
    saveRDS(seurat.sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/seurat_", sim_name, ".rds"))
    ###########################################################
    ###################### CIDR ###############################
    ###########################################################
    # t1 <- proc.time()
    # 
    # res.cidr <- scDataConstructor(as.matrix(counts(sims)))
    # res.cidr <- determineDropoutCandidates(res.cidr)
    # res.cidr <- wThreshold(res.cidr)
    # res.cidr <- scDissim(res.cidr)
    # res.cidr <- scPCA(res.cidr)
    # res.cidr <- nPC(res.cidr)
    # nCluster(res.cidr)
    # res.cidr <- scCluster(res.cidr)
    # dat$cidr_cluster <- res.cidr@clusters[match(colnames(res.cidr@tags), dat$Cell)]
    # 
    # msg <- sprintf("Completed CIDR (%d seconds). \n", floor((proc.time()-t1)[3]))
    # if(verbose) cat(msg)
    ###########################################################
    ##################### RACEID ###############################
    ###########################################################
    t1 <- proc.time()
    
    # tutorial
    # https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
    sc <- SCseq(counts(sims))
    # Cells with a relatively low total number of transcripts are discarded.
    sc <- filterdata(sc,mintotal=2000) 
    # retrieve filtered and normalized expression matrix 
    # (normalized to the minimum total transcript count across all cells retained after filtering) 
    fdata <- getfdata(sc)
    # If all genes should be used, then the parameter FSelect needs to be set to FALSE. 
    sc <- compdist(sc, metric="pearson", FSelect = FALSE)
    # sc <- clustexp(sc)
    sc <- clustexp(sc,cln=ngroup,sat=FALSE) # FUNcluster for other methods
    sc <- findoutliers(sc)
    dat$raceID_cluster <- NA
    raceID_cluster <- sc@cluster$kpart[match(names(sc@cluster$kpart), dat$Cell)]
    dat$raceID_cluster <- sc@cluster$kpart[match(dat$Cell,names(sc@cluster$kpart))]
    msg <- sprintf("Completed RaceID (%d seconds). \n", floor((proc.time()-t1)[3]))
    if(verbose) cat(msg)
    
    saveRDS(sc, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/raceID_", sim_name, ".rds"))
    
    return(dat)
}

###########################################################
######## adjusted rand index & silhouette  ################
###########################################################


sc_eval <- function(sims, dat) {
    
  res <- data.frame()
    # compute silhouette score
    # dist.matrix <- dist(t(counts(sims)))
    # scSTM.sil <- silhouette(as.numeric(as.factor(dat$scSTM_cluster)), dist.matrix)
    # seurat.sil <- silhouette(as.numeric(as.factor(dat$seurat_cluster)), dist.matrix)
    # raceid.sil <- silhouette(as.numeric(as.factor(dat$raceID_cluster)), dist.matrix)
    # cidr.sil <- silhouette(as.numeric(as.factor(dat$cidr_cluster)), dist.matrix)
    
    res <- data.frame(
        scSTM_adjR = adjustedRandIndex(dat$scSTM_cluster,sims$Group),
        Seurat_adjR = adjustedRandIndex(dat$seurat_cluster, sims$Group), 
        raceID_adjR = adjustedRandIndex(dat$raceID_cluster,sims$Group),
        CIDR_adjR = adjustedRandIndex(dat$cidr_cluster,sims$Group)#,
    #     scSTM_sil = mean(scSTM.sil[,3]),
    #     seurat_sil = mean(seurat.sil[,3]),
    #     raceID_sil = mean(raceid.sil[,3]),
    #     CIDR_sil = mean(cidr.sil[,3])
    )
    return(res)
}

