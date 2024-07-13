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


r.file <- paste0("doc/splatter/R/",list.files("doc/splatter/R/"))
sapply(r.file, source)

vcf <- mockVCF(n.samples = 3)
gff <- mockGFF()

params.batches <- newSplatPopParams(
    similarity.scale = 10,
    
    group.prob = list(c(0.2, 0.2, 0.2, 0.2, 0.2),c(0.1, 0.3, 0.2, 0.2, 0.2)),
    de.prob = 0.6,
    de.facLoc = 0.6,
    de.facScale = 0.1,
    
    batchCells = c(100, 100),
    batch.size = 3,
    batch.facLoc = 0.5,
    batch.facScale = 0.1,
    
    seed = 1
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


r.file <- paste0("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/splatter/R/",list.files("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/splatter/R/"))
sapply(r.file, source)
params <- newSplatParams()
params <- setParams(params, group.prob = c(0.1,0.2,0.3,0,0.4),
                    nGenes = 200,
                    batchCells=batchCells <- rep(30, 2))

sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)

params.batches <- newSplatPopParams(
    similarity.scale = 5,
    
    group.prob = c(0.2, 0.01, 0.1, 0.6, 0.09),
    de.prob = 0.5,
    de.facLoc = 0.5,
    de.facScale = 0.5,
    
    condition.prob = c(0.5, 0.5),
    cde.downProb = 0.2,
    cde.facLoc = 0.2,
    cde.facScale = 0.2,
    
    seed = 1
)


sim.pop.batches2 <- splatPopSimulate(
    vcf = vcf,
    gff = gff,
    params = params.batches,
    sparsify = FALSE
)
table(sim.pop.batches2$Condition, sim.pop.batches2$Sample)

sims <- cbind(sim.pop.batches, sim.pop.batches2)
sim.pop.batches <- logNormCounts(sim.pop.batches)
sim.pop.batches <- runPCA(sim.pop.batches, ncomponents = 10)
plotPCA(sim.pop.batches, colour_by = "Condition", shape_by = "Group")

df <- colData(sims) %>%
    as.data.frame() %>%
    count(Sample, Group, Condition) %>%
    pivot_wider(names_from = Condition, values_from = n, values_fill = list(n = 0))
total_1 <- sum(df$Condition1)
total_2 <- sum(df$Condition2)
result <- df %>%
    mutate(
        Proportion_1 = Condition1 / total_1,
        Proportion_2 = Condition2 / total_2,
        Ratio = Proportion_1 / Proportion_2
    )
