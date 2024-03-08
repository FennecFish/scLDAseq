# This script is design to estimate parameters from Anti-PD1
setwd("/proj/milovelab/wu/scLDAseq/scLDAseq")

library("Seurat")
library(dplyr)
library(tidyverse)
library(splatter)

set.seed(1)

# subset Anti-PD1 data into individuals and estimate parameteres
dat <- readRDS("/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/raw/cohort1_filtered.rds")
# rand.index <- sample(unique(dat$patient_id), size = 5, replace = FALSE)
rand.index <- unique(dat$patient_id)

for (i in rand.index){
  subdat <- subset(dat, subset = patient_id==i)
  saveRDS(subdat, file = paste0("data/", i, ".rds"))
  
  count <- as.matrix(subdat[["RNA"]]$counts)
  params <- splatEstimate(count)
  saveRDS(params, file = paste0("data/", i , "_params.rds"))
  cat("Parameter generated for ", i, "\n")
  rm(subdat)
}

rm(dat)
