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
library(stats)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", pattern = "sims")

# file_name <- files[1]

# calculate proportion changes using aitchison test
source("R/CompositionChange.R")
power_error_plot <- function(dat, threshold = 0.05,
                             title1, title2){
  dat <- dat %>% 
    mutate(padj = p.adjust(`p-value`, method = "fdr")) %>% # not using adjusted pvalue
    mutate(sig_change = ifelse(padj < 0.05, 1, 0))
  
  dat <- dat %>%
    mutate(truth = case_when(
      id %in% paste0("L", 4:8) & Batch %in% paste0("Batch", 1:3) ~ 0,
      sub("_.*", "", id) %in% paste0("L", 2:3) & sub("^[^_]*_", "", id) == "neg" ~ 0,
      id == "L9" & Batch %in% paste0("Batch", 1:5) ~ 0,
      TRUE ~ 1  # Default case if none of the above conditions are met
    ))
  
  power <- dat %>%
    filter(truth == 1) %>%  # Subset where the null hypothesis is false
    group_by(id) %>%
    summarise(
      total_cases = n(),  # Total cases where truth = 1
      successful_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
      power = successful_rejections / total_cases  # Proportion of successful rejections
    ) %>%
    mutate(level = sub("_.*", "", id))
  
  typeI <- dat %>%
    filter(truth == 0) %>%  # Subset where the null hypothesis is false
    group_by(id) %>%
    summarise(
      total_cases = n(),  # Total cases where truth = 1
      false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
      alpha = false_rejections / total_cases  # Proportion of successful rejections
    ) %>%
    mutate(level = sub("_.*", "", id))
  
  par(mfrow = c(1, 2))
  p <- ggplot(power, aes(x = level, y = power, group = 1)) +
    geom_line(color = "blue") +  # Connect points with lines
    geom_point(size = 3, color = "red") +  # Highlight each point
    labs(title = title1,
         x = "Level",
         y = "Power") +
    theme_minimal()
  
  t1 <- ggplot(typeI, aes(x = level, y = alpha, group = 1)) +
    geom_line(color = "blue") +  # Connect points with lines
    geom_point(size = 3, color = "red") +  # Highlight each point
    labs(title = title2,
         x = "Level",
         y = "Type I Error") +
    theme_minimal()
  
  return(grid.arrange(p, t1))
}

##### test for each sample ######

dat <- data.frame()
for(file_name in files){
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  seed <- sub("(.*)_.*$", "\\1", set_level)
  stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/scSTM_combat_filterGenes_noContent_",set_level,".rds"))
  
  res <- composition_change(stmobj)
  res <- res %>%
    mutate(
           Batch = rownames(res),
           seed = seed)
  opt <- sub("_.*", "", set_level)
  L <- sub(".*_", "", set_level)
  if(opt %in% c("pos", "neg")){res$id <- paste0(L,"_",opt)} else{res$id=L}
  dat <- rbind(dat,res)
}

write.csv(dat, file = "res/res_ilr_aitchison_combat_f_nc.csv")

dat <- read.csv("res/res_ilr_CompChange_combat_f_nc.csv")
power_error_plot(dat, title1 = "Power for Aitchison", title2 = "TypeI Error for Aitchison")

#### Using colmeans and compare group ######
files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", pattern = "sims")

composition_change <- function(Y, level, opt){
  res <- data.frame()
  x1 <- Y[(1:nrow(Y)) %% 2 == 1, ]
  x2 <- Y[(1:nrow(Y)) %% 2 == 0, ]
  
  y1 <- compositions::ilr(x1)
  y2 <- compositions::ilr(x2)
  m1 <- colMeans(y1)
  m2 <- colMeans(y2)
  d <- dim(y1)[2]
  n1 <- dim(y1)[1] 
  n2 <- dim(y2)[1]
  s1 <- (crossprod(y1) - n1 * tcrossprod(m1) ) / n1
  s2 <- (crossprod(y2) - n2 * tcrossprod(m2) ) / n2
  
  i <- 1
  s1h <- s1
  s2h <- s2
  s1inv <- solve(s1h)  ;  s2inv <- solve(s2h)
  mha <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2)
  s1h <- s1h + tcrossprod(m1 - mha)
  s2h <- s2h + tcrossprod(m2 - mha)
  s1inv <- solve(s1h) 
  s2inv <- solve(s2h)
  mhb <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2 )
  while ( sum( abs(mha - mhb) ) > 1e-6 ) {
    i <- i + 1
    mha <- mhb
    mhb <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2 )
    s2h <- s1h + tcrossprod(m1 - mhb)
    s2h <- s1h + tcrossprod(m2 - mhb)
  }
  dof <- d
  stat <- n1 * log( det(s1h) / det(s1) ) + n2 * log( det(s2h) / det(s2) )
  pvalue <- pchisq(stat, dof, lower.tail = FALSE)
  
  # summarize all data
  if(level %in% c("L2","L3")){sample = paste0(level,"_",opt)} else(sample = level)
  res <- rbind(res, c(sample,stat,pvalue,dof))
  colnames(res) <- c("id", "Statistics", "p-value","DegreeOfFreedom")
  return(res)
}

comp_change <- function(level, neg=NULL){
  
  if(level %in% c("L2", "L3")){
    if (neg){opt = "neg"} else{opt = "pos"}
    stmFiles <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/",
                           pattern = paste0(".*", opt, ".*", level, ".*"))
  } else{
    stmFiles <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/",
                           pattern = level)
  }
  
  all.Y <- vector(mode = "list")
  for (file_name in stmFiles){
    set_level <- sub("scSTM_combat_filterGenes_noContent_([^.]*)\\.rds", "\\1",  file_name)
    seed <- sub("(.*)_.*$", "\\1", set_level)
    stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/", file_name))
    sampleID <- unique(stmobj$sampleID)
    theta <- stmobj$theta
    rownames(theta) <- stmobj$DocName
    time <- stmobj$settings$covariates$X[,2]
    names(time) <- stmobj$DocName
    
    allsamp.Y <- vector(mode = "list")
    for (sample in sampleID){
      t1 <- theta[time==1 & stmobj$sampleID == sample,]
      t2 <- theta[time==2 & stmobj$sampleID == sample,]
      Y <- rbind(colMeans(t1),colMeans(t2))
      #temp <- CompDTUReg(genename = sample, Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL)
      #res <- rbind(res, temp)
      allsamp.Y[[sample]] <- Y
    }
    if(length(all.Y)==0){all.Y <- allsamp.Y}
    all.Y <- mapply(rbind, all.Y, allsamp.Y, SIMPLIFY = FALSE) # combine by each batch
  }
  
  res <- lapply(all.Y, FUN = composition_change, level = level, opt = opt)
  
  res <- lapply(names(res), function(batch) {
    df <- res[[batch]]
    df$Batch <- batch  # Add a new column with the batch name
    return(df)
  })
  
  # Combine all the data frames into a single data frame
  res <- do.call(rbind, res)
  return(res)
}

dat <- data.frame()
level <- paste0("L",2:9)
for(l in level){
  if (l %in% c("L2","L3")){
    temp1 <- comp_change(l,neg=TRUE)
    temp2 <- comp_change(l,neg=FALSE)
    temp <- rbind(temp1,temp2)
  } else{
    temp <-  comp_change(l,neg=NULL)
  }
  dat <- rbind(dat, temp)
}

write.csv(dat, file = "res/res_ilr_aitchison_combinedBatch_combat_f_nc.csv")

dat <- read.csv("res/res_ilr_aitchison_combinedBatch_combat_f_nc.csv")
power_error_plot(dat, title1 = "Power for Aitchison", title2 = "TypeI Error for Aitchison")

#### try using propeller ######
library(speckle)



source("R/thetaPosterior.R")

for(file_name in files){
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  seed <- sub("(.*)_.*$", "\\1", set_level)
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", file_name))
  stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/scSTM_combat_filterGenes_noContent_",set_level,".rds"))
  
  
  truth <- colData(sims) %>%
    as.data.frame() %>%
    count(Batch, Group, time) %>%
    pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
    mutate(ratio = `1` / `2`) %>% #calcualte ratio between time for each group
    select(Batch, Group, ratio)
  
  tolerance <- 0.05
  
  # Mutate the new 'truth' column
  truth <- truth %>%
    group_by(Batch) %>%
    mutate(truth = if (all(abs(ratio - 1) <= tolerance)) 0 else 1) %>%
    ungroup() %>%
    group_by(Batch) %>%
    slice(1)
  
  res <- composition_change(stmobj)
  L <- sub(".*_", "", set_level)
  temp <- data.frame(Batch = rownames(res),
                     level = L,
                     seed = seed,
                     pvalue = res$`p-value`,
                     truth = truth$truth)
  dat <- rbind(dat,temp)
  
  
  
  }

###### Try CompDTU #######
library(CompDTUReg)

file_name <- files[1]
set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
seed <- sub("(.*)_.*$", "\\1", set_level)
sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", file_name))
stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/scSTM_combat_filterGenes_noContent_",set_level,".rds"))

theta <- stmobj$theta
rownames(theta) <- stmobj$DocName
time <- stmobj$settings$covariates$X[,2]
s1_t1 <- theta[time==1 & stmobj$sampleID == "Batch1",]
s2_t1 <- theta[time==1 & stmobj$sampleID == "Batch2",]
s1_t2 <- theta[time==2 & stmobj$sampleID == "Batch1",]
s2_t2 <- theta[time==2 & stmobj$sampleID == "Batch2",]
Y <- rbind(s1_t1,s1_t2)
# Y <- rbind(colMeans(s1_t1),colMeans(s1_t2), colMeans(s2_t1),colMeans(s2_t1))
# Y <- rbind(colMeans(s1_t1),colMeans(s1_t2))
Y <- compositions::ilr(Y)
# Group <- c(1,2)
Group <- time[names(time) %in% rownames(Y)] 
Group <- factor(Group)
res <- CompDTUReg(genename = "test", Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL)
