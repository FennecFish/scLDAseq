setwd("/proj/milovelab/wu/scLDAseq")
set.seed(1)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(mclust)
library(lme4)
library(lmerTest)
library(dplyr)
library(tibble)
library(stats)
library(limma)

thetaPosteriorSample <- function(model, nsims=100) {
  lambda <- model$eta
  nu <- model$nu
  out <- vector(mode="list",length=nrow(lambda)) 
  for (i in 1:length(out)) {
    sigma <- nu[[i]]
    choleskydecomp <- chol(sigma)
    mat <- rmvnorm(nsims, lambda[i,],nu[[i]],choleskydecomp)
    mat <- cbind(mat, 0)
    out[[i]] <- exp(mat - row.lse(mat))
  }
  return(out)
}

rmvnorm<-function(n,mu,Sigma,chol.Sigma=chol(Sigma)) {
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol.Sigma) +c(mu))
}

row.lse <- function(mat) {
  matrixStats::rowLogSumExps(mat)
}

limma_test <- function(scSTMobj, cov) {
  
  # design matrix construction
  metadata <- data.frame(scSTMobj$settings$covariates$X) %>% select(-X.Intercept.)
  colnames(metadata) <- sub(".*\\.", "", colnames(metadata))
  xmat <- factor(metadata[[cov]])
  levels(xmat) <- c("cov1", "cov2")
  design <- model.matrix(~0+xmat)
  colnames(design) <- levels(xmat)
  
  thetasims <- thetaPosteriorSample(scSTMobj, nsims=1) # draw from posterior distribution
  thetasims <- do.call(rbind, thetasims)
  thetasims <- t(thetasims)
  # estimate within sample correlation 
  corfit <- duplicateCorrelation(thetasims,design,block=scSTMobj$sampleID)
  # fit linear model
  fit <- lmFit(thetasims,design = design, block=scSTMobj$sampleID,correlation=corfit$consensus)
  cm <- makeContrasts(cov_contrast = cov2-cov1, levels=design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, coef="cov_contrast", adjust.method = "fdr")
  res <- res[ order(row.names(res), decreasing = F),]
  return(res)
}

