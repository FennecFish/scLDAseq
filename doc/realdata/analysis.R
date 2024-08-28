# Goal: This script is used to analyze scSTM from real data PD1
# 1) compare assigned cluster from paper cluster
# 2) changes between groups (time/expansion)
# setwd("/proj/milovelab/wu/scLDAseq")
setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
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

# scSTMobj <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_TimeandResponse_Prevalence_TimeandSample.rds")
scSTMobj <- readRDS("data/scSTM_Content_TimeandResponse_Prevalence_TimeandSample.rds")
scSTMobj <- select_top_scSTM(scSTMobj)
scSTMobj <- readRDS("data/test_interaction.rds")
sims <- scSTMobj$settings$sce
K <- length(unique(sims$cellType))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
structure_plot(scSTMobj, topics = 1:K, grouping = sims$cellType, n = 5000, gap = 5)
time_effect <- estimateEffect(1:K ~ timepoint*expansion, 
                              stmobj = scSTMobj, ref <- c("Pre", "NE"),
                              sampleNames = "patient_id",
                              # sampleIDs = "BIOKEY_10",
                              uncertainty = "Global")

summary(time_effect)
plot.estimateEffect(time_effect, covariate = "timepoint", model=scSTMobj,
                    method="difference",cov.value1="Pre",cov.value2="On")

plot(time_effect, covariate = "timepoint", model = scSTMobj, 
     #method = "difference",
     method = "difference", cov.value1="Pre",cov.value2="On",
     xlab = "timepoint", moderator = "expansion", 
     moderator.value = "E", linecol = "blue", printlegend = FALSE) 

plot(time_effect, covariate = "timepoint", model = scSTMobj, 
     #method = "difference",
     method = "difference", cov.value1="Pre",cov.value2="On",
     xlab = "timepoint", moderator = "expansion", 
     moderator.value = "NE", linecol = "red", add = T, 
     printlegend = F) 
legend(0, .08, c("Liberal", "Conservative"), + lwd = 2, col = c("blue", "red"))

labelTopics(scSTM_nC_P, 4)
