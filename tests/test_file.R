library(slam)
library(SeuratObject)
library(Seurat)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
# sce <- readRDS("data/p10_p15_sub.rds")
# #subset for 300 genes, and 200 cells
# sub <- sce[500:799,10696:10895] 
# saveRDS(sub,file = "data/toydata.rds")
# source("R/STMfunctions.R")
# source("R/asSTMCorpus.R")
# source("R/readCountMatrix.R")

##### STM ########
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")


sce <- readRDS("../scLDAseq/data/toydata.rds")
dat <- prepsce(sce)
K <- 10
prevalence <- as.formula(~dat$meta$timepoint)
content <- NULL
documents  <- dat$documents
vocab <- dat$vocab
data <- dat$meta
sce <- dat$sce
sample <- "patient_id"

res <- multi_stm(documents = documents, vocab = vocab,
    K = K, prevalence = prevalence, content = NULL,
    data = data, sce = sce,
    sample = sample,
    init.type= "Spectral",
    gamma.prior= "Pooled",
    kappa.prior= "L1")

res$convergence$bound
# convergence
plot(res$convergence$bound, type = "l", ylab = "Approximate Objective", main = "Convergence")
# ###### STMinit.R ############
# source("R/spectral.R")
# 
# ###### stm.control.R ############
# source("R/STMinit.R")
# 
# 
# ### in STM.R
# source("R/stm.R")
# source("R/stm.contol.R")
# source("R/STMlncpp.R")
# 
# 
# ## E-M
# sourceCpp("src/STMCfuns.cpp")
# source("R/STMestep.R")
# source("R/STMmu.R")
# source("R/STMoptbeta.R")
# source("R/STMsigma.R")
# source("R/STMsigs.R")
# source
# # after E-M
# source("R/STMreport.R")
# 
# 
# ###### for stm function stm.R #########
# # init.type=c("Spectral", "LDA", "Random", "Custom")
# init.type= "Spectral"
# seed=NULL
# max.em.its=500
# emtol=1e-5
# verbose=TRUE
# reportevery=5
# LDAbeta=TRUE
# interactions=TRUE
# ngroups=1
# model=NULL
# # gamma.prior=c("Pooled", "L1")
# gamma.prior= "Pooled"
# sigma.prior=0
# #kappa.prior=c("L1", "Jeffreys")
# kappa.prior="L1"
# control=list()
# 
# 
# 
# ###### STMinit.R ############
# source("R/spectral.R")
# 
# ###### stm.control.R ############
# source("R/STMinit.R")
# 
# 
# ### in STM.R
# model = NULL
# source("R/STMlncpp.R")
# library(Rcpp)
# sourceCpp("src/STMCfuns.cpp")
# source("R/STMestep.R")
# 
# ## STM.control
# source("R/STMmu.R")
# source("R/STMoptbeta.R")
# source("R/STMsigma.R")
# source("R/STMsigs.R")
# 
# # after E-M
# source("R/STMreport.R")
