setwd("/proj/milovelab/wu/scLDAseq")
library(compositions)
library(Matrix)
library(dplyr)
# library(scuttle)
library(tidyverse)
# library(splatter)
# library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(CompDTUReg)
library(stats)

power_error_plot <- function(dat, threshold = 0.05,
                             title1, title2){
  dat <- dat %>%
    mutate(padj = p.adjust(pval_CompDTU, method = "fdr")) %>% # not using adjusted pvalue
    mutate(sig_change = ifelse(padj < 0.05, 1, 0))
  
  power <- dat %>%
    filter(truth == 1) %>%  # Subset where the null hypothesis is false
    group_by(nCellType) %>%
    summarise(
      total_cases = n(),  # Total cases where truth = 1
      successful_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
      power = successful_rejections / total_cases  # Proportion of successful rejections
    ) 
  
  typeI <- dat %>%
    filter(truth == 0) %>%  # Subset where the null hypothesis is false
    group_by(nCellType) %>%
    summarise(
      total_cases = n(),  # Total cases where truth = 1
      false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
      alpha = false_rejections / total_cases  # Proportion of successful rejections
    ) 
  
  par(mfrow = c(1, 2))
  p <- ggplot(power, aes(x = nCellType, y = power, group = 1)) +
    geom_line(color = "blue") +  # Connect points with lines
    geom_point(size = 3, color = "red") +  # Highlight each point
    labs(title = title1,
         x = "nCellType",
         y = "Power") +
    theme_minimal()
  
  t1 <- ggplot(typeI, aes(x = nCellType, y = alpha, group = 1)) +
    geom_line(color = "blue") +  # Connect points with lines
    geom_point(size = 3, color = "red") +  # Highlight each point
    labs(title = title2,
         x = "nCellType",
         y = "Type I Error") +
    theme_minimal()
  
  return(grid.arrange(p, t1))
}
all.Y <- readRDS("allY_V3_scSTM.rds")
all.Y <- lapply(all.Y, function(x) x[1:100, ])
##### CompDTU ########
compDTU <- function(Y){
  # Y <- Y[,2:ncol(Y)]
  Y <- compositions::ilr(Y)
  Group <- rep(c(1,2), times = nrow(Y)/2)
  Group <- factor(Group)
  rep <- nrow(Y)/10
  res <- data.frame()
  for (n in 1:rep){
    start <- n*10-9
    end <- n*10
    sub_Y <- Y[start:end,]
    sub_G <- Group[start:end]
    temp <- CompDTUReg(genename = paste0("rep-",n), Y = sub_Y, Group = sub_G, runWithME = FALSE, YInfRep = NULL)
    res <- rbind(res, temp)
  }
  return(res)
}

res <- lapply(all.Y, FUN = compDTU)
res <- Map(function(df, name) {
  df$truth <- str_split(name, "_", simplify = TRUE)[1]
  df$nCellType <- str_split(name, "_", simplify = TRUE)[2]
  df
}, res, names(res))
res <- do.call(rbind, res) %>%
  as.data.frame() %>%
  mutate(truth = ifelse(truth == "neg", 0, 1))

power_error_plot(res, title1 = "Power Plot for CompDTU", title2 = "Type I error plot for CompDTU")



######################### using V1 ###########################################
# 
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V3", pattern = "sims")
# 
# calc_irl <- function(theta, stmobj, sample){
#   rownames(theta) <- stmobj$DocName
#   time <- stmobj$settings$covariates$X[,2]
#   t1 <- theta[time==1 & stmobj$sampleID == sample,]
#   t2 <- theta[time==2 & stmobj$sampleID == sample,]
#   single_Y <- rbind(colMeans(t1), colMeans(t2))
#   single_Y <- compositions::ilr(single_Y)
#   return(single_Y)
# }
# 
# theta_posterior <- function(stmobj, nsims){
#   thetasims <- thetaPosterior(stmobj, nsims=nsims, type="Global")
#   new_lists <- vector("list", nsims)
#   for (i in 1:nsims) {
#     new_lists[[i]] <- lapply(thetasims, function(x) x[i, ])
#   }
#   new_matrices <- lapply(new_lists, function(lst) {
#     do.call(rbind, lst)
#   })
#   thetasims <- new_matrices
#   return(thetasims)
# }
# 
# ## treat each sample as a replicate, and data generated under the same level are biological replicates ##
# 
# 
# compDTU_rep <- function(level, neg=NULL){
#   
#   #if(level %in% c("L2", "L3")){sample = gsub("_\\d+_", "_", set_level)} else{sample = level}
#   if(level %in% c("L2", "L3")){
#     if (neg){opt = "neg"} else{opt = "pos"}
#     stmFiles <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/",
#                            pattern = paste0(".*", opt, ".*", level, ".*"))
#   } else{
#     stmFiles <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/",
#                            pattern = level)
#   }
#   
#   all.Y <- vector(mode = "list")
#   for (file_name in stmFiles){
#     set_level <- sub("scSTM_combat_filterGenes_noContent_([^.]*)\\.rds", "\\1",  file_name)
#     seed <- sub("(.*)_.*$", "\\1", set_level)
#     stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/", file_name))
#     sampleID <- unique(stmobj$sampleID)
#     theta <- stmobj$theta
#     rownames(theta) <- stmobj$DocName
#     time <- stmobj$settings$covariates$X[,2]
#     names(time) <- stmobj$DocName
#     
#     allsamp.Y <- vector(mode = "list")
#     for (sample in sampleID){
#       t1 <- theta[time==1 & stmobj$sampleID == sample,]
#       t2 <- theta[time==2 & stmobj$sampleID == sample,]
#       Y <- rbind(colMeans(t1),colMeans(t2))
#       #temp <- CompDTUReg(genename = sample, Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL)
#       #res <- rbind(res, temp)
#       allsamp.Y[[sample]] <- Y
#     }
#     if(length(all.Y)==0){all.Y <- allsamp.Y}
#     all.Y <- mapply(rbind, all.Y, allsamp.Y, SIMPLIFY = FALSE) # combine by each batch
#   }
#   res <- lapply(all.Y, FUN = compDTU, level = level, opt = opt)
#   
#   res <- lapply(names(res), function(batch) {
#     df <- res[[batch]]
#     df$Batch <- batch  # Add a new column with the batch name
#     return(df)
#   })
#   
#   # Combine all the data frames into a single data frame
#   res <- do.call(rbind, res)
#   return(res)
# }
# 
# compDTU <- function(Y, level, opt){
#   Y <- compositions::ilr(Y)
#   Group <- rep(c(1,2), times = nrow(Y)/2)
#   Group <- factor(Group)
#   if(level %in% c("L2","L3")){sample = paste0(level,"_",opt)} else(sample = level)
#   res <- CompDTUReg(genename = sample, Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL)
#   return(res)
# }
# 
# power_error_plot <- function(dat, threshold = 0.05,
#                              title1, title2){
#   dat <- dat %>% 
#     mutate(padj = p.adjust(pval_CompDTU, method = "fdr")) %>% # not using adjusted pvalue
#     mutate(sig_change = ifelse(pval_CompDTU < 0.05, 1, 0))
#   
#   dat <- dat %>%
#     mutate(truth = case_when(
#       gene_id %in% paste0("L", 4:8) & Batch %in% paste0("Batch", 1:3) ~ 0,
#       sub("_.*", "", gene_id) %in% paste0("L", 2:3) & sub("^[^_]*_", "", gene_id) == "neg" ~ 0,
#       gene_id == "L9" & Batch %in% paste0("Batch", 1:5) ~ 0,
#       TRUE ~ 1  # Default case if none of the above conditions are met
#     ))
#   
#   power <- dat %>%
#     filter(truth == 1) %>%  # Subset where the null hypothesis is false
#     group_by(gene_id) %>%
#     summarise(
#       total_cases = n(),  # Total cases where truth = 1
#       successful_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
#       power = successful_rejections / total_cases  # Proportion of successful rejections
#     ) %>%
#     mutate(level = sub("_.*", "", gene_id))
#   
#   typeI <- dat %>%
#     filter(truth == 0) %>%  # Subset where the null hypothesis is false
#     group_by(gene_id) %>%
#     summarise(
#       total_cases = n(),  # Total cases where truth = 1
#       false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
#       alpha = false_rejections / total_cases  # Proportion of successful rejections
#     ) %>%
#     mutate(level = sub("_.*", "", gene_id))
#   
#   par(mfrow = c(1, 2))
#   p <- ggplot(power, aes(x = level, y = power, group = 1)) +
#     geom_line(color = "blue") +  # Connect points with lines
#     geom_point(size = 3, color = "red") +  # Highlight each point
#     labs(title = title1,
#          x = "Level",
#          y = "Power") +
#     theme_minimal()
#   
#   t1 <- ggplot(typeI, aes(x = level, y = alpha, group = 1)) +
#     geom_line(color = "blue") +  # Connect points with lines
#     geom_point(size = 3, color = "red") +  # Highlight each point
#     labs(title = title2,
#          x = "Level",
#          y = "Type I Error") +
#     theme_minimal()
#   
#   return(grid.arrange(p, t1))
# }
# dat <- data.frame()
# level <- paste0("L",2:9)
# for(l in level){
#   if (l %in% c("L2","L3")){
#     temp1 <- compDTU_rep(l,neg=TRUE)
#     temp2 <- compDTU_rep(l,neg=FALSE)
#     temp <- rbind(temp1,temp2)
#   } else{
#     temp <-  compDTU_rep(l,neg=NULL)
#   }
#   dat <- rbind(dat, temp)
# }
# write.csv(dat, file = "res/res_ilr_compDTU_BatchasSample.csv")
# 
# dat <- read.csv("res/res_ilr_compDTU_BatchasSample.csv")
# power_error_plot(dat, title1 = "Power using CompDTU, using posterior mean as replicates",
#                  title2 = "Type I error using CompDTU, using posterior mean as replicates")
# 
# 
# library(stats)
# threshold <- 0.05
# dat <- dat %>% 
#   mutate(padj = p.adjust(pval_CompDTU, method = "fdr")) %>% # not using adjusted pvalue
#   mutate(sig_change = ifelse(pval_CompDTU < 0.05, 1, 0))
# 
# dat <- dat %>%
#   mutate(truth = case_when(
#     gene_id %in% paste0("L", 4:8) & Batch %in% paste0("Batch", 1:3) ~ 0,
#     sub("_.*", "", gene_id) %in% paste0("L", 2:3) & sub("^[^_]*_", "", gene_id) == "neg" ~ 0,
#     gene_id == "L9" & Batch %in% paste0("Batch", 1:5) ~ 0,
#     TRUE ~ 1  # Default case if none of the above conditions are met
#   ))
# 
# power <- dat %>%
#   filter(truth == 1) %>%  # Subset where the null hypothesis is false
#   group_by(gene_id) %>%
#   summarise(
#     total_cases = n(),  # Total cases where truth = 1
#     successful_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
#     power = successful_rejections / total_cases  # Proportion of successful rejections
#   )
# 
# ggplot(power, aes(x = level, y = power, group = 1)) +
#   geom_line(color = "blue") +  # Connect points with lines
#   geom_point(size = 3, color = "red") +  # Highlight each point
#   labs(title = "Statistical Power by Level",
#        x = "Level",
#        y = "Power") +
#   theme_minimal()
# 
# typeI <- dat %>%
#   filter(truth == 0) %>%  # Subset where the null hypothesis is false
#   group_by(gene_id) %>%
#   summarise(
#     total_cases = n(),  # Total cases where truth = 1
#     false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
#     alpha = false_rejections / total_cases  # Proportion of successful rejections
#   )
# 
# ggplot(typeI, aes(x = level, y = alpha, group = 1)) +
#   geom_line(color = "blue") +  # Connect points with lines
#   geom_point(size = 3, color = "red") +  # Highlight each point
#   labs(title = "Type I Error by Level",
#        x = "Level",
#        y = "Type I Error") +
#   theme_minimal()
# 
# ### first try to apply comDTU on all cells, and treat each cell as a sample ###
# dat <- data.frame()
# for (file_name in files){
#   set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
#   seed <- sub("(.*)_.*$", "\\1", set_level)
#   
#   stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/scSTM_combat_filterGenes_noContent_",set_level,".rds"))
#   sampleID <- unique(stmobj$sampleID)
#   theta <- stmobj$theta
#   rownames(theta) <- stmobj$DocName
#   time <- stmobj$settings$covariates$X[,2]
#   names(time) <- stmobj$DocName
#   res <- data.frame()
#   
#   for (sample in sampleID){
#     t1 <- theta[time==1 & stmobj$sampleID == sample,]
#     t2 <- theta[time==2 & stmobj$sampleID == sample,]
#     Y <- rbind(t1,t2)
#     Y <- compositions::ilr(Y) # using ilr
#     Group <- time[names(time) %in% rownames(Y)]
#     Group <- factor(Group)
#     temp <- CompDTUReg(genename = sample, Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL)
#     res <- rbind(res, temp)
#   }
#   res$Batch <- res$gene_id
#   L <- sub(".*_", "", set_level)
#   opt <- sub("_.*", "", set_level)
#   if(opt %in% c("pos", "neg")){res$gene_id <- paste0(L,"_",opt)} else{res$gene_id=L}
#   res$seed = seed = seed
#   dat <- rbind(dat,res)
# }
# write.csv(dat, file = "res/res_ilr_compDTU_CellasSample.csv")
# 
# dat <- read.csv("res/res_ilr_compDTU_CellasSample.csv")
# 
# 
# ##### Using compDTU, but use mean cell composition as samples ####
# #### extract samples from posterior distribution ######
# nsims = 10
# 
# r.file <- paste0("R/",list.files("R/"))
# sapply(r.file, source)
# 
# dat <- data.frame()
# for (file_name in files){
#   set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
#   seed <- sub("(.*)_.*$", "\\1", set_level)
#   
#   stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/scSTM_combat_filterGenes_noContent_",set_level,".rds"))
#   thetasims <- theta_posterior(stmobj, nsims)
#   sampleID <- unique(stmobj$sampleID)
#   res <- data.frame()
#   for (sample in sampleID){
#     Y <- lapply(thetasims, calc_irl, stmobj = stmobj, sample = sample)
#     Y <- do.call(rbind, Y)
#     Group <- rep(c(1,2), nrow(Y)/2)
#     Group <- factor(Group)
#     temp <- CompDTUReg(genename = sample, Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL)
#     res <- rbind(res, temp)
#   }
#   res$Batch <- res$gene_id
#   L <- sub(".*_", "", set_level)
#   opt <- sub("_.*", "", set_level)
#   if(opt %in% c("pos", "neg")){res$gene_id <- paste0(L,"_",opt)} else{res$gene_id=L}
#   res$seed = seed = seed
#   dat <- rbind(dat,res)
# }
# 
# write.csv(dat, file = "res/res_ilr_compDTU_PosteriorasSample.csv")
# power_error_plot(dat, title1 = "Power using CompDTU, using posterior mean as replicates",
#                  title2 = "Type I error using CompDTU, using posterior mean as replicates")
# 
