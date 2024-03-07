
params <- newSplatParams()
params <- setParams(params, group.prob = c(0.5,0.3,0.2),
                    de.prob = c(0.3, 0.3, 0.3), 
                    nGenes = 600, batchCells=c(500,500,500),  
                    lib.loc = 15)
sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)# create change of proportions

# we assume that the cells pre and post treatment are equal
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/length(unique(sims$Group)))
pre_prp <- matrix(c(0.5,0.2,0.3))
pre_count <- cell_count %*% t(pre_prp)


sampled_data <- colData(sims) %>%
    data.frame() %>%
    group_by(Group, Batch) %>%
    mutate(time = case_when(
        Group == "Group1" & Batch == "Batch1" & 
            row_number() %in% sample(row_number(), min(pre_count[1,1], n())) ~ 1,
        Group == "Group2" & Batch == "Batch1" & 
            row_number() %in% sample(row_number(), min(pre_count[1,2], n())) ~ 1,
        Group == "Group1" & Batch == "Batch2" & 
            row_number() %in% sample(row_number(), min(pre_count[2,1], n())) ~ 1,
        Group == "Group2" & Batch == "Batch2" & 
            row_number() %in% sample(row_number(), min(pre_count[2,2], n())) ~ 1,
        Group == "Group1" & Batch == "Batch3" & 
            row_number() %in% sample(row_number(), min(pre_count[3,1], n())) ~ 1,
        Group == "Group2" & Batch == "Batch3" & 
            row_number() %in% sample(row_number(), min(pre_count[3,2], n())) ~ 1,
        TRUE ~ 2 
    )) %>%
    ungroup() 


time_prop <- sampled_data %>% 
    data.frame() %>%
    group_by(time, Batch) %>%
    count(Group, Batch, time) %>%
    mutate(Proportion = n / sum(n))


sims$time <- sampled_data$time

##### eval 
r.file <- paste0("R/",list.files("R/"))
# r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
# sourceCpp("../stm/src/STMCfuns.cpp")
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims)
K <- 3
prevalence <- as.formula(~dat$meta$time)
content <- NULL
sce <- dat$sce
documents  <- dat$documents
vocab <- dat$vocab
meta <- dat$meta
sample <- "Batch"

res <- multi_stm(documents = documents, vocab = vocab,
                 K = K, prevalence = prevalence, content = NULL,
                 data = meta, 
                 sce = sce,
                 sample = sample,
                 init.type= "Spectral",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 control = list(gamma.maxits=3000))


# b_prep <-estimateEffect(1:5 ~ time + batch, stmobj = b, meta= meta,uncertainty = "Global")

form <- as.formula(1:5~time)
form

