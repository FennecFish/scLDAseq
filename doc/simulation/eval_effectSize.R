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

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/",
                    pattern = "*_3cellTypes_")
file_name <- files[1]
sim_name <- sub("\\_sims.rds$", "", file_name)
truth <- readRDS(file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", sim_name, "_sims.rds"))
stmobj <-  readRDS(file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/scLDAseq_", sim_name, ".rds"))
K <- stmobj$settings$dim$K
meta <- colData(truth)
# metadata <- stmobj$settings$dim
effM <- estimateMixedEffect(alt.formula = 1:K ~ time, 
                            null.formula = NULL, metadata = meta,
                            stmobj = stmobj, uncertainty = "Global", nsims = 10)

plot.fixedEffect(effM, "time")

multi_eff <- estimateEffect(1:K ~ time, stmobj = stmobj, metadata = meta, uncertainty = "Global")
summary(multi_eff)

logiteff <- estimateEffect(1:K ~ time, stmobj = stmobj, metadata = meta, uncertainty = "Global")
summary(logiteff)
plot(multi_eff, "time", model=stmobj,
     method="difference",cov.value1=1,cov.value2=2)

# truth proportion change
time_prop <- colData(truth) %>% 
  data.frame() %>%
  group_by(time) %>%
  count(Group, time) %>%
  mutate(Proportion = n / sum(n))

time_prop_sample <- colData(truth) %>% 
  data.frame() %>%
  group_by(time, Batch) %>%
  count(Group, Batch, time) %>%
  mutate(Proportion = n / sum(n))

ggplot(time_prop_sample, aes(x = Group, y = Proportion, fill = as.factor(time))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Group", y = "Proportion", fill = "time", 
       title = "Cell Type True Proportion") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") + 
  facet_grid(~Batch)
