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

lm_est <- function(scSTMobj, fixed_effect, nsims = 20) {
  K <- scSTMobj$settings$dim$K
  metadata <- data.frame(scSTMobj$settings$covariates$X) %>% select(-X.Intercept.)
  colnames(metadata) <- sub(".*\\.", "", colnames(metadata))
  formula <- paste("1:K ~", paste(fixed_effect, collapse = "+"))
  formula = as.formula(formula)
  
  response <- as.character(formula)[2]
  K <- eval(parse(text=response))
  formula <- formula(paste(as.character(formula)[c(1,3)], collapse = " "))
  termobj <- terms(formula, data=metadata)
  
  mf <- model.frame(termobj, data=metadata)
  xmat <- model.matrix(termobj,data=metadata)
  varlist <- all.vars(termobj)
  metadata <- metadata[, varlist, drop=FALSE]
  xmat <- as.data.frame(xmat)
  xmat$sample <- scSTMobj$sampleID
  
  res <-vector(mode = "list", length = nsims)
  for (n in 1:nsims) {
    thetasims <- thetaPosteriorSample(scSTMobj, nsims=1) # draw from posterior distribution
    thetasims <- do.call(rbind, thetasims)
    
    thetaLogit <- log(thetasims + 0.1/(1-thetasims + 0.1))
    
    output <- data.frame()
    for(k in K) {
      y <- thetaLogit[,k]
      full_formula <- paste("y ~", paste(fixed_effect, collapse = "+"), "+ (1|sample)")
      full_formula <- as.formula(full_formula)
      
      lm.mod = lmer(full_formula, REML = F, data = xmat)
      lm_0 <-lmer(y ~ (1|sample), REML = F, data = xmat)
      lrt <- anova(lm.mod,lm_0) %>% data.frame() %>% select(Chisq, Df, Pr..Chisq.) %>% filter(!is.na(Pr..Chisq.))
      
      est <- summary(lm.mod)$coefficients %>% 
        data.frame() %>% 
        rownames_to_column("covariate") %>% 
        filter(!covariate=="(Intercept)") %>% 
        column_to_rownames("covariate") %>%
        select(Estimate, Std..Error)
      est <- cbind(est, lrt)
      colnames(est) <- c("Estimate", "Std", "Chisq", "df", "p-value")
      est$topic <- paste0("topic",k)
      # PBmodcomp(lm.mod,lm_0,nsim=200)
      output <- rbind(output, est)
    }
    res[[n]] <- output
    # output$fdr <- p.adjust(output$`p-value`, method = "fdr")
    aggregate_results <- do.call(rbind, res) %>%
      group_by(topic) %>%
      summarise(estimate = mean(Estimate),
                Std = mean(Std),
                Chisq = mean(Chisq),
                df = mean(df)) %>%
      mutate(LRT_pValue = 1 - pchisq(Chisq, df))
    aggregate_results$fdr <- p.adjust(aggregate_results$LRT_pValue, method = "fdr")
  }
  return(aggregate_results)
}


