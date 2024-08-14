setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(Rcpp)
library(mclust)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

select_scSTM <- function(scSTMobj) {
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

cluster_scSTM <- function(scSTMobj, sims) {
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- scSTMobj$DocName
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  res_cluster <- res_cluster[match(colnames(sims), names(res_cluster))]
  return(adjustedRandIndex(res_cluster, sims$Group))
}

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/de_gene_eval/sims/", pattern = "sims*")
n <- 10
res <- vector(mode = "list")

for (file_name in files){
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/de_gene_eval/sims/sims_", set_level, ".rds"))
  truth <- rowData(sims)
  # truth <- truth[rownames(truth) %in% scSTM_C_P$vocab,]
  groupDE_cols <- grep("GroupDE", colnames(truth))
  # subset to only upregulated genes
  rows_to_keep <- apply(truth[, groupDE_cols], 1, function(row) any(row > 1))
  up_gene <- truth[rows_to_keep, ]
  K <- length(unique(sims$Group))
  # for batch
  batchDE_cols <- grep("BatchFac", colnames(truth))
  batch_up_gene <- truth[apply(truth[, batchDE_cols], 1, function(row) any(row > 1)), ]
  
  scSTM_C_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/de_gene_eval/scSTM_Content_Prevalance/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_P_name)){
    scSTM_C_P <- readRDS(scSTM_C_P_name)
    scSTM_C_P <- select_scSTM(scSTM_C_P)
    est_C_P <- labelTopics(scSTM_C_P, n = n)

    C_P <- data.frame(
      scSTM = "C_P",
      sim = set_level,
      seed = unlist(strsplit(set_level, "_"))[1],
      control = unlist(strsplit(set_level, "_"))[2],
      level = unlist(strsplit(set_level, "_"))[3],
      nCellType = unlist(strsplit(set_level, "_"))[4],
      nSample = unlist(strsplit(set_level, "_"))[5],
      num_of_top_DEgene = paste0(
        ifelse(is.na(table(est_C_P$topics %in% rownames(up_gene))["TRUE"]), 0, table(est_C_P$topics %in% rownames(up_gene))["TRUE"]), 
        "/", n*K ),
      num_of_top_Timegene = paste0(
        ifelse(is.na(table(est_C_P$covariate %in% rownames(batch_up_gene))["TRUE"]), 0, table(est_C_P$covariate %in% rownames(batch_up_gene))["TRUE"]),
        "/", n*K ),
      adjustedRandIndex = cluster_scSTM(scSTM_C_P, sims)
    )
  } else {
    C_P <- data.frame(
      scSTM = "C_P",
      sim = set_level,
      seed = unlist(strsplit(set_level, "_"))[1],
      control = unlist(strsplit(set_level, "_"))[2],
      level = unlist(strsplit(set_level, "_"))[3],
      nCellType = unlist(strsplit(set_level, "_"))[4],
      nSample = unlist(strsplit(set_level, "_"))[5],
      num_of_top_DEgene = NA,
      num_of_top_Timegene = NA,
      adjustedRandIndex = NA
    )
  }
  
 
  scSTM_nC_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/de_gene_eval/scSTM_noContent_Prevalance/scSTM_",set_level,".rds")
  if(file.exists(scSTM_nC_P_name)){
    scSTM_nC_P <- readRDS(scSTM_nC_P_name)
    scSTM_nC_P <- select_scSTM(scSTM_nC_P)
    est_nC_P <- labelTopics(scSTM_nC_P, n = n)
    
    nC_P <- data.frame(
      scSTM = "nC_P",
      sim = set_level,
      seed = unlist(strsplit(set_level, "_"))[1],
      control = unlist(strsplit(set_level, "_"))[2],
      level = unlist(strsplit(set_level, "_"))[3],
      nCellType = unlist(strsplit(set_level, "_"))[4],
      nSample = unlist(strsplit(set_level, "_"))[5],
      num_of_DEgene_prob = paste0(
        ifelse(is.na(table(est_nC_P$prob %in% rownames(up_gene))["TRUE"]), 0 , table(est_nC_P$prob %in% rownames(up_gene))["TRUE"]), 
        "/", n*K),
      num_of_DEgene_frex = paste0(
        ifelse(is.na(table(est_nC_P$frex %in% rownames(up_gene))["TRUE"]), 0 , table(est_nC_P$frex %in% rownames(up_gene))["TRUE"]), 
        "/", n*K),
      num_of_DEgene_lift = paste0(
        ifelse(is.na(table(est_nC_P$lift %in% rownames(up_gene))["TRUE"]), 0 , table(est_nC_P$lift %in% rownames(up_gene))["TRUE"]), 
        "/", n*K),
      num_of_DEgene_score = paste0(
        ifelse(is.na(table(est_nC_P$score %in% rownames(up_gene))["TRUE"]), 0 , table(est_nC_P$score %in% rownames(up_gene))["TRUE"]), 
        "/", n*K),
      adjustedRandIndex = cluster_scSTM(scSTM_nC_P, sims)
    )
    
  } else {
    nC_P <- data.frame(
      scSTM = "nC_P",
      sim = set_level,
      seed = unlist(strsplit(set_level, "_"))[1],
      control = unlist(strsplit(set_level, "_"))[2],
      level = unlist(strsplit(set_level, "_"))[3],
      nCellType = unlist(strsplit(set_level, "_"))[4],
      nSample = unlist(strsplit(set_level, "_"))[5],
      num_of_DEgene_prob = NA,
      num_of_DEgene_frex = NA,
      num_of_DEgene_lift = NA,
      num_of_DEgene_score = NA,
      adjustedRandIndex = NA
    )
  }
  
  scSTM_C_nP_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/de_gene_eval/scSTM_Content_noPrevalance/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_nP_name)){
    scSTM_C_nP <- readRDS(scSTM_C_nP_name)
    scSTM_C_nP <- select_scSTM(scSTM_C_nP)
    est_C_nP <- labelTopics(scSTM_C_nP, n = n)
    
    C_nP <- data.frame(
      scSTM = "C_nP",
      sim = set_level,
      seed = unlist(strsplit(set_level, "_"))[1],
      control = unlist(strsplit(set_level, "_"))[2],
      level = unlist(strsplit(set_level, "_"))[3],
      nCellType = unlist(strsplit(set_level, "_"))[4],
      nSample = unlist(strsplit(set_level, "_"))[5],
      num_of_top_DEgene = paste0(
        ifelse(is.na(table(est_C_nP$topics %in% rownames(up_gene))["TRUE"]), 0, table(est_C_nP$topics %in% rownames(up_gene))["TRUE"]), 
        "/", n*K ),
      num_of_top_Timegene = paste0(
        ifelse(is.na(table(est_C_nP$covariate %in% rownames(batch_up_gene))["TRUE"]), 0, table(est_C_nP$covariate %in% rownames(batch_up_gene))["TRUE"]),
        "/", n*K ),
      adjustedRandIndex = cluster_scSTM(scSTM_C_nP, sims)
    )
  } else {
    C_nP <- data.frame(
      scSTM = "C_nP",
      sim = set_level,
      seed = unlist(strsplit(set_level, "_"))[1],
      control = unlist(strsplit(set_level, "_"))[2],
      level = unlist(strsplit(set_level, "_"))[3],
      nCellType = unlist(strsplit(set_level, "_"))[4],
      nSample = unlist(strsplit(set_level, "_"))[5],
      num_of_top_DEgene = NA,
      num_of_top_Timegene = NA,
      adjustedRandIndex = NA
    )
  }
  res[[file_name]][["nC_P"]] <- nC_P
  res[[file_name]][["C_P"]] <- C_P
  res[[file_name]][["C_nP"]] <- C_nP
  
  cat(file_name, "\n")
  rm(sims)
}

saveRDS(res, file = "res/de_genes_eval/topic_definition_de.rds")

################################################################################
########################### generate plot ######################################
################################################################################

res <- readRDS("res/de_genes_eval/topic_definition_de.rds")
extract_and_combine <- function(res, element_name) {
  do.call(rbind, lapply(res, function(x) x[[element_name]]))
}

nC_P_df <- extract_and_combine(res, "nC_P")
C_nP_df <- extract_and_combine(res, "C_nP")
C_P_df <- extract_and_combine(res, "C_P")

# Function to convert proportion to numeric
convert_proportion <- function(x) {
  sapply(strsplit(x, "/"), function(y) as.numeric(y[1]) / as.numeric(y[2]))
}
################################################################################
########## plot nC_P different proportion, separate by nSample, nCellType ######
################################################################################
# Convert proportions to numeric
nC_P_df <- nC_P_df %>%
  mutate(across(starts_with("num_of_DEgene_"), convert_proportion))

# Transform data to long format
nC_P_long <- nC_P_df %>%
  rename_with(~ gsub("num_of_DEgene_", "", .), starts_with("num_of_DEgene_")) %>%
  mutate(nCellType = recode(nCellType, "c2" = "Number_of_CellType = 2", "c4" = "Number_of_CellType = 4"),
         nSample = recode(nSample, "nsample1" = "Number_of_Sample = 1", "nsample5" = "Number_of_Sample = 5")) %>%
  pivot_longer(cols = c(prob, frex, lift, score),
               names_to = "metric",
               values_to = "value")

# Create the box plot
png("res/de_genes_eval/nC_P_DEGs_proportion_fig.png", width = 2500, height = 2000, res =300)
ggplot(nC_P_long, aes(x = metric, y = value)) +
  geom_boxplot() +
  facet_wrap(~ nCellType + nSample) +
  labs(title = "Comparison of Metrics Defining Top Genes in the Topics",
       subtitle = "Without the Content Model",
       x = "Metrics",
       y = "Proportion of Overlapping Genes") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 13)
  )
dev.off()

################################################################################
########## plot C_nP different proportion, separate by nSample, nCellType ######
################################################################################

# there are a few C_nP had very bad performance. Not sure why but remove them
C_nP_df <- C_nP_df %>% filter(adjustedRandIndex > 0)

# Convert proportions to numeric
C_nP_df <- C_nP_df %>%
  mutate(across(starts_with("num_of_top_"), convert_proportion))

# Transform data to long format
C_nP_df <- C_nP_df %>%
  rename_with(~ gsub("num_of_", "", .), starts_with("num_of_")) %>%
  mutate(nCellType = recode(nCellType, "c2" = "Number_of_CellType = 2", "c4" = "Number_of_CellType = 4"),
         nSample = recode(nSample, "nsample1" = "Number_of_Sample = 1", "nsample5" = "Number_of_Sample = 5")) %>%
  rename("top_DEgene" = "CellType_Genes",
         "top_Timegene" = "Content_Genes") %>%
  pivot_longer(cols = c(CellType_Genes, Content_Genes),
               names_to = "Top_Genes",
               values_to = "value")

# Create the box plot
png("res/de_genes_eval/C_nP_DEGs_proportion_fig.png", width = 2500, height = 2000, res =300)
ggplot(C_nP_df, aes(x = Top_Genes, y = value)) +
  geom_boxplot() +
  facet_wrap(~ nCellType + nSample) +
  labs(title = "Proportion of Top Genes Defined by Different Parts of the Model",
       subtitle = "With the Content, Without the Prevalence Model",
       x = "Models",
       y = "Proportion of Overlapping Genes") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 13)
  )
dev.off()


################################################################################
########## plot C_P different proportion, separate by nSample, nCellType ######
###############################################################################
# there are a few C_nP had very bad performance. Not sure why but remove them
C_P_df <- C_P_df %>% filter(adjustedRandIndex > 0)

# Convert proportions to numeric
C_P_df <- C_P_df %>%
  mutate(across(starts_with("num_of_top_"), convert_proportion))

# Transform data to long format
C_P_df <- C_P_df %>%
  rename_with(~ gsub("num_of_", "", .), starts_with("num_of_")) %>%
  mutate(nCellType = recode(nCellType, "c2" = "Number_of_CellType = 2", "c4" = "Number_of_CellType = 4"),
         nSample = recode(nSample, "nsample1" = "Number_of_Sample = 1", "nsample5" = "Number_of_Sample = 5")) %>%
  rename("top_DEgene" = "CellType_Genes",
         "top_Timegene" = "Content_Genes") %>%
  pivot_longer(cols = c(CellType_Genes, Content_Genes),
               names_to = "Top_Genes",
               values_to = "value")

# Create the box plot
png("res/de_genes_eval/C_P_DEGs_proportion_fig.png", width = 2500, height = 2000, res =300)
ggplot(C_P_df, aes(x = Top_Genes, y = value)) +
  geom_boxplot() +
  facet_wrap(~ nCellType + nSample) +
  labs(title = "Proportion of Top Genes Defined by Different Parts of the Model",
       subtitle = "With the Content, With the Prevalence Model",
       x = "Models",
       y = "Proportion of Overlapping Genes") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 13)
  )
dev.off()
