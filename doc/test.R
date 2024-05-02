r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

library(stm)
prep <- prepDocuments(poliblog5k.docs, poliblog5k.voc,
                      poliblog5k.meta,subsample=500,
                      lower.thresh=20,upper.thresh=200)
heldout <- make.heldout(prep$documents, prep$vocab)
documents <- heldout$documents
vocab <- heldout$vocab
meta <- prep$meta

stm1<- stm(documents, vocab, 100,
           prevalence =~ rating+ s(day),
           init.type="Random",
           data=meta, max.em.its=5)
eval.heldout(stm1, heldout$missing)

res <- searchK(prep$documents, prep$vocab, K = c(5,10,15,20,30,50,100), prevalence =~ rating + s(day), data = prep$meta)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
res$results %>% as.data.frame() %>%
    transmute(K,
              `Lower bound` = lbound,
              Residuals = map_dbl(residual, "dispersion"),
              `Semantic coherence` = map_dbl(semantic_coherence, mean),
              `Held-out likelihood` = map_dbl(eval_heldout, "expected.heldout")) %>%
    gather(Metric, Value, -K) %>%
    ggplot(aes(K, Value, color = Metric)) +
    geom_line(size = 1.5, alpha = 0.7, show.legend = FALSE) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(x = "K (number of topics)",
         y = NULL,
         title = "Model diagnostics by number of topics",
         subtitle = "These diagnostics indicate that a good number of topics would be around 60")
