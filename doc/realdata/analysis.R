# Goal: This script is used to analyze scSTM from real data PD1
# 1) compare assigned cluster from paper cluster
# 2) changes between groups (time/expansion)
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

select_top_scSTM <- function(scSTMobj) {
  if(class(scSTMobj) == "selectModel") {
    all_values <- unlist(scSTMobj$bound)
    max_value <- max(all_values)
    if(length(which(all_values == max_value)) > 1){
      max_position_in_vector <- which(all_values == max_value)[1]
    } else{
      max_position_in_vector <- which(all_values == max_value)
    }
    scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
  }
  return(scSTMobj)
}

scSTMobj <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_TimeandResponse_Prevalence_TimeandSample.rds")

dat <- select_top_scSTM(scSTMobj)
sims <- dat$settings$sce
K <- length(unique(sims$cellType))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
time_effect <- estimateEffect(1:K ~ timepoint*expansion, 
                              stmobj = dat, 
                              sampleNames = "patient_id",
                              # sampleIDs = "BIOKEY_10",
                              uncertainty = "Global")

summary(time_effect)
plot.estimateEffect(time_effect, "Time", model=scSTM_nC_P,
                    method="difference",cov.value1="Time1",cov.value2="Time2")
labelTopics(scSTM_nC_P, 4)
