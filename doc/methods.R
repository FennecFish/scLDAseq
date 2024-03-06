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
library(Seurat)
library(RaceID)
library(cidr)
library(cluster)

set.seed(1)

sc_methods <- function(sims) {
    
    ngroup <- length(unique(sims$Group))
    dat <- colData(sims) %>% 
        data.frame() 

    ###########################################################
    ################## scLDAseq ###############################
    ###########################################################
    r.file <- paste0("R/",list.files("R/"))
    sapply(r.file, source)
    sourceCpp("src/STMCfuns.cpp")
    
    K <- ngroup
    
    dat <- prepsce(sims)
    prevalence <- as.formula(~dat$meta$time)
    content <- NULL
    sce <- dat$sce
    documents  <- dat$documents
    vocab <- dat$vocab
    data <- dat$meta
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
    
    max_indices <- apply(res.stm$theta, 1, which.max)
    colnames(res.stm$theta) <- paste0("topic_", 1:ncol(res.stm$theta))
    rownames(res.stm$theta) <- colnames(res.stm$mu$mu)
    res_cluster <- colnames(res.stm$theta)[max_indices]
    names(res_cluster) <- rownames(res.stm$theta)
    dat$scSTM_cluster <- res_cluster[match(names(res_cluster), dat$Cell)]
    
    ###########################################################
    #################### Seurat ###############################
    ###########################################################
    
    seurat.sims <- as.Seurat(sims, counts = "counts", data = "logcounts")
    seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 500)
    all.genes <- rownames(seurat.sims)
    seurat.sims <- ScaleData(seurat.sims, features = all.genes)
    seurat.sims <- RunPCA(seurat.sims)
    
    seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
    seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)
    dat$seurat_cluster <- Idents(seurat.sims)[match(names(Idents(seurat.sims)), dat$Cell)]

    ###########################################################
    ###################### CIDR ###############################
    ###########################################################
    res.cidr <- scDataConstructor(counts(sims))
    res.cidr <- determineDropoutCandidates(res.cidr)
    res.cidr <- wThreshold(res.cidr)
    res.cidr <- scDissim(res.cidr)
    res.cidr <- scPCA(res.cidr)
    res.cidr <- nPC(res.cidr)
    nCluster(res.cidr)
    res.cidr <- scCluster(res.cidr)
    dat$cidr_cluster <- res.cidr@clusters[match(colnames(res.cidr@tags), dat$Cell)]

    ###########################################################
    ##################### RACEID ###############################
    ###########################################################
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
    dat$raceID_cluster <- sc@cluster$kpart[match(names(sc@cluster$kpart), dat$Cell)]
    
    return(dat)
}
