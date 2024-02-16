library(slam)
library(SeuratObject)
library(Seurat)
library(SingleCellExperiment)

# sce <- readRDS("data/p10_p15_sub.rds")
# #subset for 300 genes, and 200 cells
# sub <- sce[500:799,10696:10895] 
# saveRDS(sub,file = "data/toydata.rds")
source("R/STMfunctions.R")
source("R/asSTMCorpus.R")
source("R/readCountMatrix.R")
sce <- readRDS("data/toydata.rds")
dat <- prepsce(sce)
K <- 3
prevalence <- as.formula(~dat$meta$timepoint)
content <- NULL
documents  <- dat$documents
vocab <- dat$vocab
data <- dat$meta
sample <- "patient_id"


###### for stm function stm.R #########
# init.type=c("Spectral", "LDA", "Random", "Custom")
init.type= "Spectral"
seed=NULL
max.em.its=500
emtol=1e-5
verbose=TRUE
reportevery=5
LDAbeta=TRUE
interactions=TRUE
ngroups=1
model=NULL
# gamma.prior=c("Pooled", "L1")
gamma.prior= "Pooled"
sigma.prior=0
#kappa.prior=c("L1", "Jeffreys")
kappa.prior="L1"
control=list()



###### STMinit.R ############
source("R/spectral.R")

###### stm.control.R ############
source("R/STMinit.R")
