#' Heldout Likelihood by Document Completion
#' 
#' Tools for making and evaluating heldout datasets.
#' 
#' These functions are used to create and evaluate heldout likelihood using the
#' document completion method.  The basic idea is to hold out some fraction of
#' the words in a set of documents, train the model and use the document-level
#' latent variables to evaluate the probability of the heldout portion. See the
#' example for the basic workflow.
#' 
#' @aliases make.heldout eval.heldout
#' @param documents the documents to be modeled (see \code{\link{stm}} for format).
#' @param vocab the vocabulary item
#' @param N number of docs to be partially held out
#' @param proportion proportion of docs to be held out.
#' @param seed the seed, set for replicability
#' @sample 

#' @export
make.heldout.sce <- function(sce , N=200, 
                         proportion=.1, seed=NULL, sample) {
  if(!is.null(seed)) set.seed(seed)
  
    if(is.null(sce)) stop("Please provide a SingleCellExperiment Object")
    
    # use whichever give more samples N or proportion
    ndoc <- ncol(sce)
    pie <- ifelse(proportion * ndoc >= N, proportion, N/ndoc)
    
    # Split the meta data by sample ID
    groups <- split(seq_along(sce[[sample]]), sce[[sample]])
    
    # Sample from each group certain proportion
    sampled_indices <- lapply(groups, function(x) 
        {sample(x, size = floor(length(x) * pie))})
    sampled_indices <- sort(unlist(sampled_indices, use.names = FALSE))
    
    sce_sub <- sce[,sampled_indices]
  return(list(sce_sub = sce_sub, keep = sampled_indices))
}

#' #' @export
#' eval.heldout <- function(model, missing) {
#'   heldout <- vector(length=length(missing$index))
#'   ntokens <- vector(length=length(missing$index))
#'   beta <- lapply(model$beta$logbeta, exp)
#'   bindex <- model$settings$covariates$betaindex[missing$index]
#'   for(i in 1:length(missing$index)) {
#'     docid <- missing$index[i]
#'     words <- missing$docs[[i]][1,]
#'     probs <- model$theta[docid,]%*%beta[[bindex[i]]][,words]
#'     probs <- rep(probs, missing$docs[[i]][2,])
#'     heldout[i] <- mean(log(probs)) 
#'     ntokens[i] <- sum(missing$docs[[i]][2,])
#'   }
#'   out <- list(expected.heldout=mean(heldout, na.rm=TRUE), doc.heldout=heldout,
#'               index=missing$index, ntokens=ntokens) #the mean here to deal with 0 length docs
#'   return(out)
#' }