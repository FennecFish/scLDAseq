# this script is directly calculating composition change using simulated data
# the purpose of this script is to evaluate different compositional data DE testing methods

setwd("/proj/milovelab/wu/scLDAseq")
library(compositions)
library(Matrix)
library(dplyr)
library(tidyverse)
library(SingleCellExperiment)
library(CompDTUReg)
library(stats)
library(ggplot2)
library(gridExtra)
power_error_compDTU_plot <- function(dat, threshold = 0.05,
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

power_error_ait_plot <- function(dat, threshold = 0.05,
                                     title1, title2){
  dat <- dat %>%
    mutate(padj = p.adjust(dat$`p-value`, method = "fdr")) %>% # not using adjusted pvalue
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

all.Y <- readRDS("allY_V3_sims.rds")
all.Y <- lapply(all.Y, function(x) x[1:100, ])

##### CompDTU ########
compDTU_v3 <- function(Y){
  Y <- Y[,2:ncol(Y)]
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

compDTU_v1 <- function(Y){
  Batch <- unique(Y$Batch)
  res <- data.frame()
  for (b in Batch){
    subY <- Y %>%
      filter(Batch == b) 
    subY <- subY[,3:ncol(subY)]
    subY <- compositions::ilr(subY)
    Group <- rep(c(1,2), times = nrow(subY)/2)
    Group <- factor(Group)
    temp <- CompDTUReg(genename = b, Y = subY, Group = Group, runWithME = FALSE, YInfRep = NULL)
    res <- rbind(res, temp)
  }
  return(res)
}

res <- lapply(all.Y, FUN = compDTU_v3)
res <- Map(function(df, name) {
  df$truth <- str_split(name, "_", simplify = TRUE)[1]
  df$nCellType <- str_split(name, "_", simplify = TRUE)[2]
  # df$nCellType <- name
  df
}, res, names(res))
res <- do.call(rbind, res) %>%
  as.data.frame() %>%
  # mutate(truth = case_when(
  #   nCellType %in% c("L4", "L5", "L6", "L7", "L8") & gene_id %in% c("Batch1", "Batch2", "Batch3") ~ 0,
  #   nCellType == "L9" & gene_id %in% c("Batch1", "Batch2", "Batch3", "Batch4", "Batch5") ~ 0,
  #   TRUE ~ 1))
  mutate(truth = ifelse(truth == "neg", 0, 1))

power_error_plot(res, title1 = "Power Plot for CompDTU", title2 = "Type I error plot for CompDTU")

######## using aitchison test ##########
# this doesn't work right now....
# this function is the account for singularity matrix
inv_cov <- function(sigma, epsilon = 1e-10) {
  sigobj <- try(solve(sigma), silent = TRUE)
  if (inherits(sigobj, "try-error")) {
    regularized_sigma <- sigma + diag(epsilon, nrow(sigma))
    siginv <- solve(regularized_sigma)
  } else {
    siginv <- solve(sigma)
  }
  return(siginv)
}

aittest <- function(y1, y2){
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
  
  s1inv <- inv_cov(s1h)
  s2inv <- inv_cov(s2h)
  mha <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2)
  s1h <- s1h + tcrossprod(m1 - mha)
  s2h <- s2h + tcrossprod(m2 - mha)
  s1inv <- inv_cov(s1h) 
  s2inv <- inv_cov(s2h)
  mhb <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2 )
  while ( sum( abs(mha - mhb) ) > 1e-6 ) {
    i <- i + 1
    mha <- mhb
    mhb <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2 )
    s2h <- s1h + tcrossprod(m1 - mhb)
    s2h <- s1h + tcrossprod(m2 - mhb)
  }
  dof <- d
  stat <- n1 * log(abs(det(s1h) / det(s1))) + n2 * log( abs(det(s2h) / det(s2)))
  pvalue <- pchisq(stat, dof, lower.tail = FALSE)
  # summarize all data
  res <- c(stat,pvalue,dof)
  names(res) <- c("Statistics", "p-value","DegreeOfFreedom")
  return(res)
}
comp_change <- function(Y){
  x1 <- Y[Y$time==1,2:ncol(Y)]
  x2 <- Y[Y$time==2,2:ncol(Y)]
  y1 <- compositions::ilr(x1)
  y2 <- compositions::ilr(x2)
  rep <- nrow(y1)/10
  
  res <- data.frame()
  for (n in 1:rep){
    start <- n*10-9
    end <- n*10
    sub_y1 <- y1[start:end,]
    sub_y2 <- y2[start:end,]
    temp <- aittest(sub_y1, sub_y2)
    res <- rbind(res, temp)
    colnames(res) <- c("Statistics", "p-value","DegreeOfFreedom")
  }
  return(res)
}
res <- lapply(all.Y, FUN = comp_change)
res <- Map(function(df, name) {
  df$truth <- str_split(name, "_", simplify = TRUE)[1]
  df$nCellType <- str_split(name, "_", simplify = TRUE)[2]
  # df$nCellType <- name
  df
}, res, names(res))
res <- do.call(rbind, res) %>%
  as.data.frame() %>%
  # mutate(truth = case_when(
  #   nCellType %in% c("L4", "L5", "L6", "L7", "L8") & gene_id %in% c("Batch1", "Batch2", "Batch3") ~ 0,
  #   nCellType == "L9" & gene_id %in% c("Batch1", "Batch2", "Batch3", "Batch4", "Batch5") ~ 0,
  #   TRUE ~ 1))
  mutate(truth = ifelse(truth == "neg", 0, 1))
power_error_ait_plot(res, title1 = "Power Plot for Aitchison Test", title2 = "Type I error plot for Aitchison Test")




######### construct dataset for sims V1 ###########
# control <- c("neg", "pos")
# nCellType <- paste0("c",c(3,5,8))
# level <- expand.grid(control, nCellType) %>% mutate(level = paste0(Var1,"_",Var2)) %>% select(level)
# level <- level$level
level <- paste0("L", 4:9)

all.Y <- vector(mode = "list")
for (l in level){
  files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", pattern = l)
  res <- data.frame()
  for (file_name in files){
    sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", file_name))
    dat <- colData(sims) %>%
      as.data.frame() %>%
      dplyr::count(Batch, Group, time) %>%
      dplyr::group_by(time) %>%  
      dplyr::mutate(Proportion = n / sum(n)) %>%
      dplyr::select(Batch, time, Proportion) %>%
      dplyr::group_by(time, Batch) %>%
      dplyr::mutate(id = row_number()) %>%
      tidyr::pivot_wider(
        names_from = id,
        values_from = Proportion,
        names_prefix = "K_")
    
    cat(file_name, "\n")
    res <- rbind(res, dat)
  }
  all.Y[[l]] <- res
}
saveRDS(all.Y, file = "allY_V1_sims.rds")


######## construct the dataset for sims V3 ########
control <- c("neg", "pos")
nCellType <- paste0("c",c(3,5,8))
level <- expand.grid(control, nCellType) %>% mutate(level = paste0(Var1,"_",Var2)) %>% select(level)
level <- level$level

all.Y <- vector(mode = "list")
for (l in level){
  files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", pattern = l)
  res <- data.frame()
  for (file_name in files){
    sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/V1/sims/", file_name))
    dat <- colData(sims) %>%
      as.data.frame() %>%
      dplyr::count(Batch, Group, time) %>%
      dplyr::group_by(time) %>%  
      dplyr::mutate(Proportion = n / sum(n)) %>%
      dplyr::select(time, Proportion) %>%
      dplyr::group_by(time) %>%
      dplyr::mutate(id = row_number()) %>%
      tidyr::pivot_wider(
        names_from = id,
        values_from = Proportion,
        names_prefix = "K_")
    
    cat(file_name, "\n")
    res <- rbind(res, dat)
  }
  all.Y[[l]] <- res
}
saveRDS(all.Y, file = "all_Y.rds")

