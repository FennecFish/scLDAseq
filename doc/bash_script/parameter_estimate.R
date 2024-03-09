# This script is design to estimate parameters from Anti-PD1
setwd("/proj/milovelab/wu/scLDAseq/scLDAseq")

library("Seurat")
library(dplyr)
library(tidyverse)
library(splatter)

set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
index <- as.character(args[1])

cat(index, "\n")
number_part <- as.numeric(gsub("[^0-9]", "", index))
cat(number_part, "\n")
patientID <- paste0("BIOKEY_",number_part)
cat(patientID, "\n")
# subset Anti-PD1 data into individuals and estimate parameteres
dat <- readRDS("/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/raw/cohort1_filtered.rds")
subdat <- subset(dat, subset = patient_id == patientID)
subcount <- as.matrix(subdat[["RNA"]]$counts)
params <- splatEstimate(subcount)
saveRDS(params, file = paste0("data/", patientID , "_params.rds"))
cat("Parameter generated for ", patientID, "\n")

# # rm(subdat)
# # 
# # # rand.index <- sample(unique(dat$patient_id), size = 5, replace = FALSE)
# # # rand.index <- unique(dat$patient_id)
# # # write.csv(rand.index, file= "doc/bash_script/patientID.csv", row.names = FALSE)
# # for (i in rand.index){
# #   subdat <- subset(dat, subset = patient_id==i)
# #   saveRDS(subdat, file = paste0("data/", i, ".rds"))
# # 
# #   count <- as.matrix(subdat[["RNA"]]$counts)
# #   params <- splatEstimate(count)
# #   saveRDS(params, file = paste0("data/", i , "_params.rds"))
# #   cat("Parameter generated for ", i, "\n")
# #   rm(subdat)
# # }
# # 
# # rm(dat)
