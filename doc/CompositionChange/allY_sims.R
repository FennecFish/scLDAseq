setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/scSTM/")
all.Y <- vector(mode = "list")
for (file_name in files){
  set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file_name)
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/neg_V3/sims/sims_", set_level, ".rds"))
  

  dat <- colData(sims) %>%
    as.data.frame() %>%
    dplyr::count(Batch, Group, time) %>%
    dplyr::group_by(Batch, time) %>%
    dplyr::mutate(Proportion = n / sum(n)) %>%
    dplyr::select(Batch, time, Proportion) %>%
    dplyr::group_by(time, Batch) %>%
    dplyr::mutate(id = row_number()) %>%
    tidyr::pivot_wider(
      names_from = id,
      values_from = Proportion,
      names_prefix = "K_")
  all.Y[[set_level]] <- dat
  cat(file_name, "\n")
}
saveRDS(all.Y, file = "allY_V2_sims.rds")