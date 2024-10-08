#' Assists the user in selecting the best STM model.
#' 
#' Discards models with the low likelihood values based on a small number of EM
#' iterations (cast net stage), then calculates semantic coherence,
#' exclusivity, and sparsity (based on default STM run using selected
#' convergence criteria) to allow the user to choose between models with high
#' likelihood values. 
#' 
#' @param documents The documents to be modeled.  Object must be a list of with
#' each element corresponding to a document.  Each document is represented as
#' an integer matrix with two rows, and columns equal to the number of unique
#' vocabulary words in the document.  The first row contains the 1-indexed
#' vocabulary entry and the second row contains the number of times that term
#' appears.
#' 
#' This is similar to the format in the \pkg{lda} package except that
#' (following R convention) the vocabulary is indexed from one. Corpora can be
#' imported using the reader function and manipulated using the
#' \code{\link{prepDocuments}}.
#' @param vocab Character vector specifying the words in the corpus in the
#' order of the vocab indices in documents. Each term in the vocabulary index
#' must appear at least once in the documents.  See
#' \code{\link{prepDocuments}} for dropping unused items in the vocabulary.
#' @param K A positive integer (of size 2 or greater) representing the desired
#' number of topics. Additional detail on choosing the number of topics in
#' details.
#' @param prevalence A formula object with no response variable or a matrix
#' containing topic prevalence covariates.  Use \code{s()}, \code{ns()} or
#' \code{bs()} to specify smooth terms. See details for more information.
#' @param content A formula containing a single variable, a factor variable or
#' something which can be coerced to a factor indicating the category of the
#' content variable for each document.
#' @param runs Total number of STM runs used in the cast net stage.
#' Approximately 15 percent of these runs will be used for running a STM until
#' convergence.
#' @param data Dataset which contains prevalence and content covariates.
#' @param init.type The method of initialization.  Must be either Latent
#' Dirichlet Allocation (LDA), Dirichlet Multinomial Regression Topic Model
#' (DMR), a random initialization or a previous STM object.
#' @param seed Seed for the random number generator. \code{stm} saves the seed
#' it uses on every run so that any result can be exactly reproduced.  Setting
#' the seed here simply ensures that the sequence of models will be exactly the
#' same when respecified.  Individual seeds can be retrieved from the component
#' model objects.
#' @param max.em.its The maximum number of EM iterations.  If convergence has
#' not been met at this point, a message will be printed.
#' @param emtol Convergence tolerance.  EM stops when the relative change in
#' the approximate bound drops below this level.  Defaults to .001\%.
#' @param verbose A logical flag indicating whether information should be
#' printed to the screen.
#' @param frexw Weight used to calculate exclusivity
#' @param net.max.em.its Maximum EM iterations used when casting the net
#' @param netverbose Whether verbose should be used when calculating net
#' models.
#' @param M Number of words used to calculate semantic coherence and
#' exclusivity.  Defaults to 10.
#' @param N Total number of models to retain in the end. Defaults to .2 of
#' runs.
#' @param to.disk Boolean. If TRUE, each model is saved to disk at the current
#' directory in a separate RData file.  This is most useful if one needs to run
#' \code{multiSTM()} on a large number of output models.
#' @param \dots Additional options described in details of stm.
#' @return \item{runout}{List of model outputs the user has to choose from.
#' Take the same form as the output from a stm model.} \item{semcoh}{Semantic
#' coherence values for each topic within each model in runout}
#' \item{exclusivity}{Exclusivity values for each topic within each model in
#' runout.  Only calculated for models without a content covariate}
#' \item{sparsity}{Percent sparsity for the covariate and interaction kappas
#' for models with a content covariate.} 
#' @examples
#' 
#' \dontrun{
#' 
#' temp<-textProcessor(documents=gadarian$open.ended.response, metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' docs<-out$documents
#' vocab<-out$vocab
#' meta <-out$meta
#' set.seed(02138)
#' mod.out <- selectModel(docs, vocab, K=3, prevalence=~treatment + s(pid_rep), 
#'                        data=meta, runs=5)
#' plotModels(mod.out)
#' selected<-mod.out$runout[[1]]
#' }
#' @export
selectModel <- function(sce , K, sample = NULL,
                        prevalence=NULL, content=NULL,
                        max.em.its=100, verbose=TRUE, # init.type = "TopicScore",
                        emtol= 1e-06, seed=NULL, 
                        ts_runs = 10, random_run = 20, frexw=.7, 
                        net.max.em.its=5, netverbose=TRUE, M=10, N=NULL,
                        to.disk=F, control=list(), ...){
  if(!is.null(seed)) set.seed(seed)

  runs <- ts_runs + random_run + 1

  if(is.null(N)){
    N <-  round(.2*runs)
  }
  
  if(runs < 2){
    stop("Number of runs must be two or greater.")
  }
  
  if(runs < N){
    stop("Number in the net must be greater or equal to the number of final models.")
  }
    
    args <- prepsce(sce)
    documents <- args$documents
    vocab <- args$vocab
    data <- args$meta
    sce <- args$sce
    # divide runs between using Poisson NMF (Random), TopicScore and Spectral
    # if (runs >=3 ) {
    #     random_run <- ceiling(runs/2)-1
    #     ts_run <- runs-random_run-1
    # }else {
    #     ts_run <- runs - 1
    #     random_run <- 0
    # }
    # 
    
  seedout <- NULL
  likelihood <- NULL
  cat("Casting net \n")
  if(ts_runs > 0){
      for(i in 1:ts_runs){
          cat(paste(i, "models in net \n"))
          tryCatch({
              mod.out <- scSTMseq(sce = sce, documents = documents, vocab = vocab, data = data, 
                                  sample = sample, K = K, 
                                  prevalence=prevalence, content=content, init.type="TopicScore",
                                  max.em.its=net.max.em.its, emtol=emtol, verbose=netverbose,...)
              seedout[i] <- mod.out$settings$seed
              likelihood[i] <- mod.out$convergence$bound[length(mod.out$convergence$bound)]
          }, error = function(e) {
              cat("An error occurred in TopicScore model", i + ts_runs + 1, ": ", e$message, "\n")
          })
      }
  }

  
  # running spectral
  cat(paste(ts_runs + 1, "models in net \n"))
  tryCatch({
      mod.out <- scSTMseq(sce = sce, documents = documents, vocab = vocab, data = data,
                          sample = sample, K,
                          prevalence=prevalence, content=content, init.type="Spectral",
                          max.em.its=net.max.em.its, emtol=emtol, verbose=netverbose,...)
      likelihood[ts_runs + 1] <- mod.out$convergence$bound[length(mod.out$convergence$bound)]
      seedout[ts_runs + 1] <- mod.out$settings$seed
  }, error = function(e) {
      cat("An error occurred in spectral model", ts_runs + 1, ": ", e$message, "\n")
  })
  
  # Random run
  if(random_run > 0){
      for(i in 1:random_run){
          cat(paste(i + ts_runs + 1, "models in net \n"))
          tryCatch({
              mod.out <- scSTMseq(sce = sce, documents = documents, vocab = vocab, data = data,
                                  sample = sample, K,
                                  prevalence=prevalence, content=content, init.type="Random",
                                  max.em.its=net.max.em.its, emtol=emtol, verbose=netverbose,...)
              likelihood[ts_runs + i + 1] <- mod.out$convergence$bound[length(mod.out$convergence$bound)]
              seedout[ts_runs + i + 1] <- mod.out$settings$seed
              #TRUE # return TRUE
          }, error = function(e) {
              cat("An error occurred in random model", i + ts_runs + 1, ": ", e$message, "\n")
              #FALSE # Return FALSE if an error occurred
          })
          # if (!result) {
          #     next  # Skip to the next iteration if there was an error
          # }
      }
  }

  keep <- order(likelihood, decreasing=T)[1:N]
  keepseed <- seedout[keep]
  cat("Keep the following seed", keepseed, "\n")
  cat("Running select models \n")
  runout <- list()
  semcoh <- list()
  exclusivity <- list()
  sparsity <- list()
  bound <- list()
  
  for(i in 1:length(keepseed)){
    cat(paste(i, "select model run \n"))
    initseed <- keepseed[i]
    if(is.na(initseed)){next} # if seed is NA, then skip this run
    initseed_index <- which(seedout == initseed)
    if (initseed_index <= ts_runs) {
        # If the index is the last of the run, do Spectral
      init_type <- "TopicScore"
    }  else if (initseed_index == ts_runs + 1) {
        # otherwise do default
        init_type <- "Spectral"
    }
    else{
      init_type <- "Random"
    }
    

    tryCatch({
        mod.out <- scSTMseq(sce = sce, documents = documents, vocab = vocab, data = data, 
                            sample = sample, K = K, 
                            prevalence = prevalence, content = content, init.type = init_type, 
                            seed = initseed, max.em.its = max.em.its, emtol = emtol, 
                            verbose = verbose, ...)
        runout[[i]] <- mod.out
        bound[[i]] <- max(mod.out$convergence$bound)
        
        if(to.disk==T){
            mod <- mod.out
            save(mod, file=paste("runout", i, ".RData", sep=""))
        }
        semcoh[[i]] <- semanticCoherence(mod.out, documents, M)
        if(length(mod.out$beta$logbeta)<2){
            exclusivity[[i]] <- exclusivity(mod.out, M=M, frexw=.7)
            sparsity[[i]] = "Sparsity not calculated for models without content covariates"
        }
        if(length(mod.out$beta$logbeta)>1){
            exclusivity[[i]] = "Exclusivity not calculated for models with content covariates"
            kappas <- t(matrix(unlist(mod.out$beta$kappa$params), ncol=length(mod.out$beta$kappa$params)))
            topics <-mod.out$settings$dim$K
            numsparse = apply(kappas[(K+1):nrow(kappas),], 1,function (x) sum(x<emtol))
            sparsity[[i]] = numsparse/ncol(kappas)
        }
    }, error = function(e) {
        cat("An error occurred in final model with seed", initseed, "and initial type", 
            init_type, ": ", e$message, "\n")
        runout[[i]] <- keepseed[i]
        bound[[i]] <- NA
        semcoh[[i]] <- NA
        exclusivity[[i]] <- NA
        sparsity[[i]] <- NA
    })
  }
  
  out <- list(runout=runout, bound = bound, semcoh=semcoh, exclusivity=exclusivity, sparsity=sparsity)
  class(out) <- "selectModel"
  return(out)
}