# Check to see if signature genes in topics
# align with DE genes in sims
setwd("/proj/milovelab/wu/scLDAseq")
library(ggplot2)
library(dplyr)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/sims/", pattern = "sims*")
file_name <- files[1]
set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/sims/sims_", set_level, ".rds"))
scSTMobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark_V2/scSTM_noContent_Prevalance/scSTM_",set_level,".rds"))

if(class(scSTMobj) == "selectModel") {
  all_values <- unlist(scSTMobj$bound)
  max_value <- max(all_values)
  max_position_in_vector <- which(all_values == max_value)
  scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
}

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

# identify the true DE genes 
sim_gene <- rowData(sims)
groupde_columns <- grep("^GroupDE", colnames(sim_gene), value = TRUE)

sim_gene[learned_label$frex[3,],]

learned_label <- labelTopics(scSTMobj)
prob_label <- learned_label$prob
