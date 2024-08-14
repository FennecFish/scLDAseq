# compare composition change 
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(stm)
library(scater)
library(scran)

max_scSTM <- function(scSTMobj){
  all_values <- unlist(scSTMobj$bound)
  max_value <- max(all_values)
  max_position_in_vector <- which(all_values == max_value)
  scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
  return(scSTMobj)
}
process_scSTM <- function(scSTMobj) {
  if(class(scSTMobj) == "selectModel") {
    all_values <- unlist(scSTMobj$bound)
    max_value <- max(all_values)
    max_position_in_vector <- which(all_values == max_value)
    scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
  }
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- colnames(scSTMobj$mu$mu)
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/single_sample_benchmark/scSTM_Content_Prevalance/"
scSTMobj.neg <- readRDS(paste0(path, "scSTM_1721108559_neg_L3_c5_nsample1.rds"))

scSTMobj.neg <- max_scSTM(scSTMobj.neg)

res <- list()
# sum of normal distribution of time 1
for (k in seq(scSTMobj.neg$settings$dim$K-1)){
  eta.t1 <- scSTMobj.neg$eta[scSTMobj.neg$settings$covariates$betaindex==1,k]
  nu.t1 <- scSTMobj.neg$nu[scSTMobj.neg$settings$covariates$betaindex==1]
  mean.t1 <- sum(eta.t1)
  var.t1 <- sum(sapply(nu.t1, function(mat) mat[k, k]))
  
  eta.t2 <- scSTMobj.neg$eta[scSTMobj.neg$settings$covariates$betaindex==2,k]
  nu.t2 <- scSTMobj.neg$nu[scSTMobj.neg$settings$covariates$betaindex==2]
  mean.t2 <- sum(eta.t2)
  var.t2 <- sum(sapply(nu.t2, function(mat) mat[k, k]))
  
  diff.mean <- mean.t1- mean.t2
  diff.var <- var.t1 + var.t2
  z <- diff.mean/sqrt(diff.var)
  pval <- 2 * (1 - pnorm(abs(z)))
  # use a t-test
  sample1 <- rnorm(1000, mean = mean.t1, sd = var.t1)
  sample2 <- rnorm(1000, mean = mean.t2, sd = var.t2)
  
  # Perform a two-sample t-test
  res[[k]] <- t.test(sample1, sample2)
  
}

z.test(
  x,
  y = NULL,
  alternative = "two.sided",
  mu = 0,
  sigma.x = NULL,
  sigma.y = NULL,
  conf.level = 0.95
)
