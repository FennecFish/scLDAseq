setwd("/proj/milovelab/wu/scLDAseq/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)

set.seed(1)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

################# evaluation #############################

stm.file <- list.files(path = "res/", pattern = "stmRes")
res.dat <- data.frame()
plot_list1 = list()
plot_list2 = list()

for(index in stm.file){

  res.stm <- readRDS(paste0("res/",index))
  samp <- sub("(BIOKEY_[0-9]+).*", "\\1", index)
  sims <- readRDS(paste0("data/", samp, "_sims.rds"))
  swap <- grepl("oppo", index)
  
  max_indices <- apply(res.stm$theta, 1, which.max)
  colnames(res.stm$theta) <- paste0("topic_", 1:ncol(res.stm$theta))
  rownames(res.stm$theta) <- colnames(res.stm$mu$mu)
  res_cluster <- colnames(res.stm$theta)[max_indices]
  names(res_cluster) <- rownames(res.stm$theta)
  
  if (swap) {
    sampled_data <- colData(sims) %>%
      data.frame() %>%
      mutate(new_time = ifelse(time == 1 & Batch %in% c("Batch1", "Batch2"), 2, 1)) %>%
      mutate(time = ifelse(Batch %in% c("Batch1", "Batch2"), new_time, time))
    sims$time <- sampled_data$time
    samp <- paste0(samp, "_oppo")
  }
  
  sims <- sims[,res.stm$DocName]
  scSTM_cluster <- res_cluster[match(names(res_cluster), colnames(sims))]
  scSTM_adjR = adjustedRandIndex(scSTM_cluster,sims$Group)
  
  true_prop <- colData(sims) %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n)) %>%
    mutate(time = ifelse(time ==1, "pre-trt","post-trt"))
  true_prop$Batch <- gsub("Batch", "Sample", true_prop$Batch)
  true_prop$Group <- gsub("Group", "Cluster", true_prop$Group)
  
  # png(paste0("res/true_prop_",samp,".png"), width = 2000, height = 3000, res = 300)
  p <- ggplot(true_prop, aes(x = Group, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Group", y = "Proportion", fill = "time", 
         title = "Cell Type True Proportion",
         subtitle = samp) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") + 
    facet_grid(~Batch) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
    xlab(NULL)
  
  plot_list1[[samp]] = p

  
  res_dat <- colData(sims) %>% 
    data.frame() %>% 
    mutate(assigned_cluster = scSTM_cluster)
  
  res_prop <- res_dat %>% 
    group_by(time, Batch) %>%
    count(assigned_cluster) %>%
    mutate(Proportion = n / sum(n)) %>%
    mutate(time = ifelse(time ==1, "pre-trt","post-trt"))
  
  res_prop$Batch <- gsub("Batch", "Sample", res_prop$Batch)
  
  # png(paste0("res/stm_prop_",samp,".png"), height = 2000, width = 3000, res = 300)
  p2<-ggplot(res_prop, aes(x = assigned_cluster, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "assigned_cluster", y = "Proportion", fill = "time", 
         title = paste0("Cell Type Assigned Proportion with adjusted Rand Index=", round(scSTM_adjR, digits=3)),
         subtitle = samp) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    facet_grid(~Batch) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
    xlab(NULL)
  # dev.off()
  plot_list2[[samp]] = p2
  
  cat("Finished", samp, "\n")
  # meta <- colData(sims) %>% 
  #   data.frame()
  # library(stm)
  # est_b1 <- estimateEffect(1:res.stm$settings$dim$K ~ time, stmobj = res.stm, meta= meta,uncertainty = "Global")
  
  # res.dat <- rbind(res.dat, c(scSTM_adjR,scSTM.sil))
  
}

# Another option: create pdf where each page is a separate plot.
pdf("res/true_prop_plots.pdf")
for (index in stm.file) {
  samp <- sub("(BIOKEY_[0-9]+).*", "\\1", index)
  swap <- grepl("oppo", index)
  if(swap){
    samp <- paste0(samp, "_oppo")
  }
  print(plot_list1[[samp]])
}
dev.off()

pdf("res/stm_prop_plots.pdf")
for (index in stm.file) {
  samp <- sub("(BIOKEY_[0-9]+).*", "\\1", index)
  swap <- grepl("oppo", index)
  if(swap){
    samp <- paste0(samp, "_oppo")
  }
  print(plot_list2[[samp]])
}
dev.off()

files <- list.files(path = "res/", pattern = "^adjR_.*\\.csv$")
adj_dat <- data.frame()
for(i in files){
  csv <- read.csv(paste0("res/",i))
  sim_name <- sub(".*silhoutte_(.*?)\\.csv$", "\\1", i)
  rownames(csv) <- sim_name
  adj_dat <- rbind(adj_dat, csv)
}