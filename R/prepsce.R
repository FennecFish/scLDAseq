# from stm
read.slam.doc <- function(corpus) {
    #convert a simple triplet matrix to list format.
    if(!inherits(corpus, "simple_triplet_matrix")) stop("corpus is not a simple triplet matrix")
    if (inherits(corpus,"TermDocumentMatrix")) {
        non_empty_docs <- which(slam::col_sums(corpus) != 0)
        documents <- ijv.to.doc(corpus[,non_empty_docs]$j, corpus[,non_empty_docs]$i, corpus[,non_empty_docs]$v) 
        names(documents) <- corpus[,non_empty_docs]$dimnames[[1]]
    } else {
        non_empty_docs <- which(slam::row_sums(corpus) != 0)
        documents <- ijv.to.doc(corpus[non_empty_docs,]$i, corpus[non_empty_docs,]$j, corpus[non_empty_docs,]$v) 
        names(documents) <- corpus[non_empty_docs,]$dimnames[[1]]
    }
    return(documents=documents)
}

# require SingleCellExperiment, with meta data in ColData
prepsce <- function(sce, sample = NULL, lower.thresh=1, upper.thresh=Inf, 
                    verbose=TRUE){
    #Functions:
    # 1) Detect and renumber zero-indexed data.
    # 2) Detect and renumber missing terms
    # 3) Remove words appearing only in [lower.thresh] documents.
    # 4) Remove words appear in upper.thresh or more of the documents.
    
    # Can also optionally subsample the data.
    
    # check if used SingleCellExperiment
    if((!class(sce)[1]=="SingleCellExperiment")){
        stop("The input data needs to be a SingleCellExperiment Object")
    }
    
    # extract documents (cell) and vocab (genes) 
    documents <- read.slam.doc(slam::as.simple_triplet_matrix(t(assays(sce)$counts)))
    # browser()
    vocab <- rownames(sce)
    
    #error check for inputs
    if((is.null(documents))){
        stop("No cells found in your data")
    }
    if(is.null(vocab)){
        stop("No genes found in your data")
    }
    if(!inherits(documents,"list")) {
        stop("documents must be a list in stm() format.  See ?stm() for format.  
          See ?readCorpus() for tools for converting from popular formats")
    }
    # 
    # if(!is.null(subsample)) {
    #     index <- sample(1:length(documents), subsample)
    #     documents <- documents[index]
    #     if(!is.null(meta)) meta <- meta[index, , drop = FALSE] 
    # }
    # browser()
    #check that there are no 0 length documents
    len <- unlist(lapply(documents, length))
    if(any(len==0)) {
        stop("Some documents have 0 length.  Please check input. 
          See ?prepDocuments() for more info.")
    } 
    
    triplet <- doc.to.ijv(documents) #this also fixes the zero indexing.
    nms <- names(documents) 
    documents <- ijv.to.doc(triplet$i, triplet$j, triplet$v)
    names(documents) <- nms
    meta <- colData(sce)[nms,]
    docs.removed <- c()
    
    #Detect Missing Terms
    miss.vocab <- NULL
    vocablist <- sort(unique(triplet$j))
    wordcounts <- tabulate(triplet$j)
    if(length(vocablist)>length(vocab)) {
        stop("Your documents object has more unique features than your vocabulary file has entries.")
    } 
    if(length(vocablist)<length(vocab)) {
        if(verbose) cat("Detected Missing Terms, renumbering \n")
        miss.vocab <- vocab[-vocablist]
        vocab <- vocab[vocablist]
        new.map <- cbind(vocablist, 1:length(vocablist))
        documents <- lapply(documents, function(d) {
            nm <- names(d)
            d[1,] <- new.map[match(d[1,], new.map[,1]),2]
            names(d) <- nm
            return(d)
        })
        wordcounts <- wordcounts[vocablist]
    }
    
    #Remove Words Appearing Only n Times
    toremove <- which(wordcounts <= lower.thresh | wordcounts >= upper.thresh)
    keepers <- which(wordcounts > lower.thresh & wordcounts < upper.thresh)
    droppedwords <- c(miss.vocab,vocab[toremove])
    if(length(toremove)) {
        if(verbose) cat(sprintf("Removing %i of %i genes (%i of %i tokens) due to frequency \n", 
                                length(toremove), length(wordcounts), sum(wordcounts[toremove]), sum(wordcounts)))
        vocab <- vocab[-toremove]
        remap <- 1:length(keepers)
        for(i in 1:length(documents)) {
            doc <- documents[[i]]
            dockeep <- doc[1,]%in%keepers
            doc <- doc[,dockeep,drop=FALSE]
            doc[1,] <- remap[match(doc[1,], keepers)]
            documents[[i]] <- doc
            if(ncol(doc)==0) docs.removed <- c(docs.removed,i)
        }
        if(length(docs.removed)) {
            if(verbose) cat(sprintf("Removing %i Cells with No Gene Counts \n", length(docs.removed)))
            documents <- documents[-docs.removed]
        }
        toprint <- sprintf("Your Dataset now has %i cells, %i genes and %i tokens.", 
                           length(documents), length(vocab), sum(wordcounts[keepers]))
        if(verbose) cat(toprint)
    }
    
    # re-organize meta data to remove cells that has been removed
    if(!is.null(docs.removed) & !is.null(meta)){
        meta <- meta[-docs.removed, , drop = FALSE]
    }
    # browser()
    # re-organize sce to remove vocab and docs
    sce <- sce[vocab, rownames(meta)]
    #recast everything as an integer
    documents <- lapply(documents, function(x) matrix(as.integer(x), nrow=2))
    
    if(!is.null(sample)) {
        sample_list <- meta[sample]
        return(list(documents=documents, vocab=vocab, 
                    sample = sample_list, meta=meta, sce=sce,
                    words.removed=droppedwords, docs.removed=docs.removed, 
                    tokens.removed=sum(wordcounts[toremove]), wordcounts=wordcounts))
    } else {
        return(list(documents=documents, vocab=vocab, 
                    sample = NULL, meta=meta, sce=sce, 
                    words.removed=droppedwords, docs.removed=docs.removed, 
                    tokens.removed=sum(wordcounts[toremove]), wordcounts=wordcounts))
    }
    
    
   
}
