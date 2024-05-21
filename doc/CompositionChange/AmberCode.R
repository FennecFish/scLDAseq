library(stageR)
library(dplyr)
library(tidyr)
library(stringr)

val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

load('changed_transcripts.RData')

alpha <- seq(0,1,.001)

# load(paste0('TranscriptResChange',2,'.RData'))
# AllPvals <- results
# load(paste0('TranscriptResChange',2,'T.RData'))
AllPvalsT <- results

AllPvals <- results

AllPvalsT <- AllPvalsT %>% filter(change==2)


for(m in c('CompDTU','CompDTUme')){
  res <- list()
  for(combo in 1:100){
    pvals <- AllPvals[AllPvals$combo==combo,]
    pScreen <- pvals$pval_CompDTU
    names(pScreen) <- pvals$gene_id
    
    pvalsT <- AllPvalsT[AllPvalsT$combo==combo,]
    pConfirmation <- matrix(pvalsT$pval_CompDTU,ncol=1,
                            dimnames = list(pvalsT$tx_id,'transcript'))
    tx2gene <- data.frame(row.names = pvalsT$tx_id,
                          transcript = pvalsT$tx_id,
                          gene=pvalsT$gene_id)
    
    
    pConfirmation <- pConfirmation[!is.na(pConfirmation[,'transcript']),,drop=FALSE]
    tx2gene <- tx2gene[tx2gene$transcript %in% rownames(pConfirmation),]
    temp <- tx2gene %>% group_by(gene) %>% summarise(n=n()) %>% filter(n==1)
    tx2gene <- tx2gene[!(tx2gene$gene %in% temp$gene),]
    pScreen <- pScreen[names(pScreen) %in% tx2gene$gene]
    pConfirmation <- pConfirmation[rownames(pConfirmation) %in% tx2gene$transcript,,drop=FALSE
    ]
    stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
    stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=alpha[val])
    padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)
    padj$combo <- combo
    
    results <- padj
    
    results$reject <- ifelse(!is.na(results$transcript),results$transcript < alpha[val],0)
    results$reject_gene <- ifelse(!is.na(results$transcript),1,0)
    results$changed <- results$txID %in% trans
    
    res[[combo]] <- results
  }
  results <- do.call('rbind',res)
  
  
  TPR <- sum(results$changed & results$reject) / sum(results$changed)
  FPR <- sum(!results$changed & results$reject) / sum(!results$changed)
  
  results_in <- results[results$reject_gene==1,]
  
  TPR_in <- sum(results_in$changed & results_in$reject) / sum(results_in$changed)
  FPR_in <- sum(!results_in$changed & results_in$reject) / sum(!results_in$changed)
  
  results <- data.frame(TPR=TPR,FPR=FPR,alpha=alpha[val],method=m,TPR_in=TPR_in,FPR_in=FPR_in)
  
  save(results,file=paste0('results/AllPvals',m,alpha[val],'.RData'))
  
}

