# this script is to generate simulated data for multiple patients
# using splatPop
setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library("scater")
library(SingleCellExperiment)
library(MASS)
library(VariantAnnotation)
library(checkmate)

required_packages <- c(
    "BiocGenerics", "BiocParallel", "checkmate", "crayon", "edgeR", 
    "fitdistrplus", "grDevices", "locfit", "matrixStats", "methods", 
    "rlang", "S4Vectors", "scuttle", "stats", "SummarizedExperiment", 
    "utils", "withr"
)

# Function to install and load required packages
install_and_load <- function(packages) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg)
        }
        library(pkg, character.only = TRUE)
    }
}

bioc_packages <- c("BiocGenerics", "BiocParallel", "edgeR", "S4Vectors", "scuttle", "SummarizedExperiment")
for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg)
    }
}

# Install CRAN packages
cran_packages <- setdiff(required_packages, bioc_packages)
install_and_load(cran_packages)

# Load Bioconductor packages
sapply(bioc_packages, function(pkg) library(pkg, character.only = TRUE))


r.file <- paste0("doc/rev_splatter/",list.files("doc/rev_splatter/"))
sapply(r.file, source)

vcf <- mockVCF(n.samples = 3)
gff <- mockGFF(n.genes = 50)

params.batches <- newSplatPopParams(
    similarity.scale = 10,
    
    group.prob = list(c(0.2, 0.2, 0.2, 0.2, 0.2),c(0.1, 0.3, 0.2, 0.2, 0.2)),
    de.prob = 0.6,
    de.facLoc = 0.6,
    de.facScale = 0.1,
    
    batchCells = c(200, 200),
    batch.size = 3,
    batch.facLoc = 0.5,
    batch.facScale = 0.1,
    
    seed = 1,
    nGenes = 500
)

# condition.prob needs to be the same length of cell type
sim.pop.batches <- splatPopSimulate(
    vcf = vcf,
    gff = gff,
    params = params.batches,
    sparsify = FALSE
)

prop.table(table(sim.pop.batches$Group,sim.pop.batches$Sample, sim.pop.batches$Batch),
           margin = c(2, 3))
sim.pop.batches <- logNormCounts(sim.pop.batches)
sim.pop.batches <- runPCA(sim.pop.batches,ncomponents = 10)
plotPCA(sim.pop.batches, colour_by = "Sample", shape_by = "Group")
plotPCA(sim.pop.batches[,sim.pop.batches$Sample == "sample_1"], colour_by = "Group", shape_by = "Batch")

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

sim.pop.batches <- sim.pop.batches[rowSums(counts(sim.pop.batches)) != 0,]
sim.pop.batches <- sim.pop.batches[, colSums(counts(sim.pop.batches)) != 0]

ngroup <- length(unique(sim.pop.batches$Group))

scSTM.mod <- selectModel(sce = sim.pop.batches,
                         K = ngroup, prevalence = ~Batch, content = ~Batch,
                         N = 1, ts_runs = 1, random_run = 1,
                         max.em.its = 100, net.max.em.its = 5)

