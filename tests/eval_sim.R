setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
set.seed(1)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

sce <- readRDS("data/sim_single_sample.rds")
prop.table(table(sce$Group, sce$time))
dat <- prepsce(sce)
K <- 5
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
data <- dat$meta
# sample <- "patient_id"


res <- multi_stm(documents = documents, vocab = vocab,
           K = K, prevalence = prevalence, content = NULL,
           data = data, 
           sample = NULL,
           init.type= "Spectral",
           gamma.prior= "Pooled",
           kappa.prior= "L1")


# plot(PrevFit, type="summary", xlim=c(0,1))
# plot(PrevFit, type="labels", topics=c(1,2,3))
# plot(PrevFit, type="hist")
# plot(PrevFit, type="perspectives", topics=c(1,2))


prep <- estimateEffect(1:2 ~ time, PrevFit, meta=stm_out$meta, 
                       uncertainty="Global")
plot(prep,"time")
summary(prep)

theta <- PrevFit$theta
rownames(theta) <- colnames(sim_new)
theta1 <- theta[sim_new$time==1,]
theta2 <- theta[sim_new$time==2,]

cp1 <- apply(theta1, 1, which.max) %>% 
    data.frame(topic = .) %>%
    rownames_to_column("cell")
prop.table(table(cp1$topic))

cp2 <- apply(theta2, 1, which.max) %>% 
    data.frame(topic = .) %>%
    rownames_to_column("cell")
prop.table(table(cp2$topic))

pheatmap::pheatmap(theta1)


# dist1 <- distHellinger(theta1)
# dist1 <- as.dist(dist1)
# 
# dist2 <- distHellinger(theta2)
# dist2 <- as.dist(dist2)
# distHellinger(t(theta1),t(theta2))
# dt<-make.dt(PrevFit,meta=stm_out$meta)
# dist <- distHellinger(PrevFit$theta)
# dist <- as.dist(dist)




kc <- pam(dist, 2)
kc_res <- kc$clustering
names(kc_res) <- dt$cellname


adjustedRandIndex(kc_res,sim$Group)

sim$topic <- factor(kc_res)
plotPCA(sim, colour_by = "Step")
plotPCA(sim, colour_by = "topic")
plotPCA(sim, colour_by = "seurat_cluster")


pca <- sim@int_colData$reducedDims@listData$PCA

sim_new <- logNormCounts(sim_new)
sim_new <- runPCA(sim_new)
plotPCA(sim_new, colour_by = "time")
An error occurred
Please try again later
Contact SupportClose
Full-text Access 