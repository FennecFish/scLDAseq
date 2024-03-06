setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(mclust)
library(Seurat)
library(RaceID)
library(cidr)
library(cluster)

set.seed(1)

###########################################################
######## adjusted rand index & silhouette  ################
###########################################################

sc_eval <- function(sims, dat) {
    
    # compute silhouette score
    scSTM.sil <- silhouette(as.numeric(as.factor(dat$scSTM_cluster)), dist(t(counts(sims))))
    seurat.sil <- silhouette(as.numeric(as.factor(dat$seurat_cluster)), dist(t(counts(sims))))
    raceid.sil <- silhouette(as.numeric(as.factor(dat$raceID_cluster)), dist(t(counts(sims))))
    cidr.sil <- silhouette(as.numeric(as.factor(dat$cidr_cluster)), dist(t(counts(sims))))
    
    res <- data.frame(
        scSTM_adjR = adjustedRandIndex(dat$scSTM_cluster,sims$Group),
        Seurat_adjR = adjustedRandIndex(dat$seurat_cluster, sims$Group), 
        raceID_adjR = adjustedRandIndex(dat$raceID_cluster,sims$Group),
        CIDR_adjR = adjustedRandIndex(dat$cidr_cluster,sims$Group),
        scSTM_sil = mean(scSTM.sil[,3]),
        seurat_sil = mean(seurat.sil[,3]),
        raceID_sil = mean(raceid.sil[,3]),
        CIDR_sil = mean(cidr.sil[,3])
    )
    return(res)
}




res_dat <- colData(sims) %>% 
    data.frame() %>% 
    mutate(assigned_cluster = res_cluster)

res_prop <- res_dat %>% 
    group_by(time, Batch) %>%
    count(assigned_cluster) %>%
    mutate(Proportion = n / sum(n))




summary(cidr.sil)
plot(cidr.sil)
# png("../res/2sample_5type_assigned_prop.png", height = 2000, width = 2000, res = 300)
ggplot(res_prop, aes(x = assigned_cluster, y = Proportion, fill = as.factor(time))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "assigned_cluster", y = "Proportion", fill = "time", 
         title = "Cell Type Assigned Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    facet_grid(~Batch)
# dev.off()
# access the true proportion distribution
true_prop <- dat %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))


