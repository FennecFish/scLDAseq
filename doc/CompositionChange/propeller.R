setwd("/proj/milovelab/wu/scLDAseq")
library(compositions)
library(Matrix)
library(dplyr)
library(tidyverse)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(speckle)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", pattern = "sims")

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
nsims = 2

dat <- data.frame()
for (file_name in files){
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  seed <- sub("(.*)_.*$", "\\1", set_level)
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", file_name))
  
  tolerance <- 0.05
  truth <- colData(sims) %>%
    as.data.frame() %>%
    count(Batch, Group, time) %>%
    pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0)) %>%
    mutate(ratio = `1` / `2`) %>% #calcualte ratio between time for each group
    select(Batch, Group, ratio) 
  
  truth <- truth %>%
    group_by(Batch) %>%
    mutate(truth = if (all(abs(ratio - 1) <= tolerance)) 0 else 1) %>%
    ungroup() %>%
    group_by(Batch) %>%
    slice(1)

  stmobj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/scSTM_combat_f_nc/scSTM_combat_filterGenes_noContent_",set_level,".rds"))
  
  cluster <- vector(mode = "list", length = nsims)
  thetasims <- thetaPosterior(stmobj, nsims=nsims, type="Global")
  new_lists <- vector("list", nsims)
  for (i in 1:nsims) {
    new_lists[[i]] <- lapply(thetasims, function(x) x[i, ])
  }
  new_matrices <- lapply(new_lists, function(lst) {
    do.call(rbind, lst)
  })
  
  thetasims <- new_matrices
  cluster <- lapply(thetasims, function(mat) {
    process_scSTM(stmobj, theta = mat)
  })
  
  vec.cluster <- unlist(cluster)
  grp <- rep(stmobj$settings$covariates$X[,2], nsims)
  
  transform <- "logit"
  prop.list <- getTransformedProps(clusters, sample, transform)
  
  sampleID <- unique(stmobj$sampleID)
  res <- data.frame()
  for (sample in sampleID){
    CellName <- stmobj$DocName[stmobj$sampleID==sample]
    sub.cluster <- vec.cluster[names(vec.cluster) %in% CellName]
    sub.grp <- grp[names(grp) %in% CellName]
    samp <- rep(1:nsims, each = length(CellName))
    samp <- paste0("s",samp,"-",sub.grp)
    
    transform <- "arcsin"
    prop.list <- getTransformedProps(sub.cluster, samp, transform)
    
    group.coll <- table(samp, sub.grp)
    design <- matrix(as.integer(group.coll != 0), ncol=ncol(group.coll))
    colnames(design) <- colnames(group.coll)
    
    contrasts <- c(1,-1)
    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions
    
    fit <- lmFit(prop.trans, design)
    fit.cont <- contrasts.fit(fit, contrasts=contrasts)
    fit.cont <- eBayes(fit.cont, robust=TRUE, trend=FALSE)
    
    temp <- propeller(clusters = sub.cluster, group = sub.grp, sample = samp, trend = FALSE, robust = TRUE) %>%
      mutate(CellType = BaselineProp.clusters,
             Batch = sample)
    res <- rbind(res, temp)
  }
  L <- sub(".*_", "", set_level)
  temp <- data.frame(Batch = res$Batch,
                     level = L,
                     seed = seed,
                     pValue = res$P.Value,
                     FDR = res$FDR)
  temp$truth <- truth$truth[match(temp$Batch,truth$Batch)]
  dat <- rbind(dat,temp)
}

write.csv(dat, file = "res/res_propeller_PosteriorasSample.csv")


library(stats)
threshold <- 0.05
dat <- dat %>% 
  mutate(padj = p.adjust(pvalue, method = "fdr")) %>%
  mutate(sig_change = ifelse(padj < 0.05, 1, 0))

power <- dat %>%
  filter(truth == 1) %>%  # Subset where the null hypothesis is false
  group_by(level) %>%
  summarise(
    total_cases = n(),  # Total cases where truth = 1
    successful_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
    power = successful_rejections / total_cases  # Proportion of successful rejections
  )

ggplot(power, aes(x = level, y = power, group = 1)) +
  geom_line(color = "blue") +  # Connect points with lines
  geom_point(size = 3, color = "red") +  # Highlight each point
  labs(title = "Statistical Power by Level",
       x = "Level",
       y = "Power") +
  theme_minimal()

typeI <- dat %>%
  filter(truth == 0) %>%  # Subset where the null hypothesis is false
  group_by(level) %>%
  summarise(
    total_cases = n(),  # Total cases where truth = 1
    false_rejections = sum(sig_change == 1),  # Cases where the null was correctly rejected
    alpha = false_rejections / total_cases  # Proportion of successful rejections
  )

ggplot(typeI, aes(x = level, y = alpha, group = 1)) +
  geom_line(color = "blue") +  # Connect points with lines
  geom_point(size = 3, color = "red") +  # Highlight each point
  labs(title = "Type I Error by Level",
       x = "Level",
       y = "Type I Error") +
  theme_minimal()

