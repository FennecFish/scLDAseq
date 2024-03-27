thetaPosteriorSample <- function(model, nsims=100) {
    lambda <- model$eta
    nu <- model$nu
    out <- vector(mode="list",length=nrow(lambda)) 
    for (i in 1:length(out)) {
        sigma <- nu[[i]]
        choleskydecomp <- chol(sigma)
        mat <- rmvnorm(nsims, lambda[i,],nu[[i]],choleskydecomp)
        mat <- cbind(mat, 0)
        out[[i]] <- exp(mat - row.lse(mat))
    }
    return(out)
}

estimateMixedEffect <- function(alt.formula, null.formula = NULL,
                           stmobj, metadata=NULL, 
                           # slope = FALSE,
                           # sampleNames = NULL, sampleIDs = NULL, # if id = null, then combine all samples
                           uncertainty=c("Global", "None"), 
                           nsims=25, prior=NULL) {
    origcall <- match.call()
    thetatype <- match.arg(uncertainty)
    if(thetatype=="None") nsims <- 1 #override nsims for no uncertainty
    
    #Step 1: Extract the fixed effect formula and do some error checking
    ##
    if(!inherits(alt.formula,"formula")) stop("formula must be a formula object.")
    if(!is.null(null.formula) & !inherits(null.formula,"formula")) stop("formula must be a formula object.")
    if(!is.null(metadata) & !is.data.frame(metadata)) metadata <- as.data.frame(metadata)
    termobj <- terms(alt.formula, data=metadata)
    if(attr(termobj, "response")==1){
        #if a response is specified we have to parse it and remove it.
        # as.character of a formula turns
        # dv ~ iv
        # into:  c("~", "dv", "iv")
        response <- as.character(alt.formula)[2] #second object is the response in this cases
        K <- eval(parse(text=response))
        # browser()
        if(!(posint(K) && max(K)<=stmobj$settings$dim$K)) stop("Topics specified as response in formula must be a set of positive integers equal to or less than the number of topics in the model.")   
        #now we reconstruct the formula removing the response
        alt.formula <- formula(paste(as.character(alt.formula)[c(1,3)], collapse = " "))
 
        #the above used to be the below code but the use got deprecated.
        #as.formula(as.character(formula)[c(1,3)])
        termobj <- terms(alt.formula, data=metadata)
    } else {
        K <- 1:stmobj$settings$dim$K
    }
    
    # now we subset the metadata to include only sample specified
    # if(!is.null(sampleIDs)) {
    #     if(is.null(sampleNames)) stop("Please specify the colname name for sampleIDs in the input metadata")
    #     if(!sampleIDs %in% stmobj$sampleID) stop("Sample ids specified must be exactly the same as the ones used in the model")
    #     metadata = metadata[metadata[,sampleNames] == sampleIDs,]
    #     subDocName <- stmobj$DocName[stmobj$sampleID %in% sampleIDs]
    #     stmobj <- STMsubset(stmobj, subDocName) 
    # } 
    
    mf <- model.frame(termobj, data=metadata)
    xmat <- model.matrix(termobj,data=metadata)
    varlist <- all.vars(termobj)
    if(!is.null(metadata)) {
        data <- metadata[, varlist, drop=FALSE]
    } else {
        templist <- list()
        for(i in 1:length(varlist)) {
            templist[[i]] <- get(varlist[i])
        }
        data <- data.frame(templist)
        names(data) <- varlist
        rm(templist)
    }
    metadata <- data
    rm(data)
    
    
    ##
    #Step 2: Compute the QR decomposition
    ##
    # all the models here are essentially just OLS regressions
    # becuase we will run them many times we want to cache the 
    # expensive components in advance.
    # if(!is.null(prior)) {
    #     if(!is.matrix(prior)) {
    #         prior <- diag(prior, nrow=ncol(xmat))
    #     } 
    #     if(ncol(prior)!=ncol(xmat)) stop("number of columns in prior does not match columns in design matrix")
    #     prior.pseudo <- chol(prior)
    #     xmat <- rbind(xmat,prior.pseudo)
    # }
    # repeat the design matrix nsims times
     xmat <- do.call(rbind, replicate(nsims, xmat, simplify=FALSE))
    # qx <- qr(xmat)
    # if(qx$rank < ncol(xmat)) {
    #     prior <- diag(1e-5, nrow=ncol(xmat))
    #     prior.pseudo <- chol(prior)
    #     xmat <- rbind(xmat,prior.pseudo)
    #     qx <- qr(xmat)
    #     warning("Covariate matrix is singular.  See the details of ?estimateEffect() for some common causes.
    #          Adding a small prior 1e-5 for numerical stability.")
    # }
    xmat <- as.data.frame(xmat)
    xmat$sample <- rep(stmobj$sampleID, nsims)
    xmat <- xmat[,-1]
    ##  
    #Step 3: Calculate Coefficients
    ##
    
    # first simulate theta 
    storage <- vector(mode="list", length=nsims)
    for(i in 1:nsims) {
        # 3a) simulate theta
        if(thetatype=="None") thetasims <- stmobj$theta
        else {
            thetasims <- thetaPosteriorSample(stmobj, nsims=1)
            thetasims <- do.call(rbind, thetasims)
            storage[[i]] <- thetasims
        }
    }
    storage <- do.call(rbind, storage)
    # transform theta into logit of theta
    thetaLogit <- log(storage/(1-storage))
    
    # 3b) perform linear mixed effect model
    output <- vector(mode="list", length=length(K))
    for(k in K) {
        y <- thetaLogit[,k]
         lm.mod <- mix.lm(y, xmat, alt.formula, null.formula)
         output[[k]] <- lm.mod 
         names(output) <- paste0("topic",K)
    }
    ##
    #Step 4: Return Values
    ##

    # Process 'output' for 'plrt'
    plrt <- lapply(K, function(k) as.data.frame(output[[k]]$LRT))
    names(plrt) <- paste0("topic", K)
    
    # Process 'output' for 'param'
    param <- lapply(K, function(k) list(est = output[[k]]$coef[, 1], # Assuming you want all rows, first column for estimates
                                        std = output[[k]]$coef[, 2], # Assuming you want all rows, second column for std
                                        vcov = output[[k]]$vcov))
    names(param) <- paste0("topic", K)
    
    # Process 'output' for 'fixEffect'
    fixEffect <- lapply(K, function(k) output[[k]]$coef)
    names(fixEffect) <- paste0("topic", K)
    
    toreturn <- list(parameters=param, LRT = plrt,
                     fixedEffect = fixEffect, topics=K,
                     call=origcall, uncertainty=thetatype, 
                     fixed.formula=alt.formula, data=metadata,
                     modelframe=mf, varlist=varlist)
    class(toreturn) <- "estimateMixedEffect"
    return(toreturn)
}

# a function for performing mixed linear effect model
mix.lm <- function(y, xmat, alt.formula, null.formula){
    xdat <- xmat
    xdat$y <- y
    # varlist <- c(varlist, "(1|sample)")
    alt.form <- as.formula(paste("y", 
                                 paste(c(alt.formula,"(1|sample)"), 
                                       collapse = "+")))
    # random intercept 
    int.lmer = lmer(alt.form, REML = F, data = xdat)
    
    # null model
    if (!is.null(null.formula)) {
        null.formula <- formula(paste(as.character(null.formula)[c(1,3)], collapse = " "))
        null.form <- as.formula(paste("y", 
                                     paste(c(null.formula,"(1|sample)"), 
                                           collapse = "+")))
    } else{
        null.form <- as.formula("y ~ (1|sample)")
    }
    # null.x <- varlist[!(varlist %in% "time")]
    # null.formula <- as.formula(paste("y~", 
    #                                  paste(null.x, collapse = " + ")))
    null.fit <- lmer(null.form, REML = F, data = xdat)
    # perform LRT 
    res.lrt <- anova(null.fit, int.lmer)
    res.lm <- summary(int.lmer, ddf = "Satterthwaite")
    out <- list(coef = res.lm$coefficients, 
         resduals = res.lm$residuals,
         vcov = res.lm$vcov,
         LRT = res.lrt)
    out
}
    # slp.lmer = lmer(y ~ xmat$time + (1 + xmat$time|xmat$sample), REML = F)
    # testing if random slope is necessary
    # lm.res <- anova(int.lmer, slp.lmer)
    # if (lm.res$`Pr(>Chisq)`[2] < 0.05) {slope = TRUE}
    # if (slope){
    #     lmer.fit <- slp.lmer
    #     null.fit <- lmer(y ~ (1|xmat$sample), REML = F)
    #     res <- anova(null.fit, slp.lmer)
    # } else{
    #     lmer.fit <- int.lmer
    #     null.fit <- lmer(y ~ (1|xmat$sample), REML = F)
    #     res <- anova(null.fit, int.lmer)
    # }

# A function for performing simple linear regression with a cached QR decomposition
# this should be lighter weight than lm().  Works with summary.qr.lm() to give
# vcov calculations etc.
# qr.lm <- function(y, qx) {
#     if(length(y)!=nrow(qx$qr)) {
#         #probably don't match because of a prior
#         if(length(y)!=(nrow(qx$qr)-ncol(qx$qr))) stop("number of covariate observations does not match number of docs")
#         #if it passes this check its the prior. thus
#         y <- c(y,rep(0, ncol(qx$qr)))
#     }
#     beta <- solve.qr(qx, y)
#     residuals <- qr.resid(qx,y)
#     fitted.values <- qr.fitted(qx,y)
#     df.residual <- length(fitted.values) - qx$rank
#     out <- list(coefficients=beta, residuals=residuals, 
#                 fitted.values=fitted.values, 
#                 df.residual=df.residual, rank=qx$rank, qr=qx)
#     out 
# }
#this function rewrites the summary.lm() function
# to calculate from our reduced regression
# summary.qr.lm <- function (object) {
#     z <- object
#     p <- z$rank
#     rdf <- z$df.residual
#     
#     Qr <- object$qr
#     n <- nrow(Qr$qr)
#     p1 <- 1L:p
#     r <- z$residuals
#     f <- z$fitted.values
#     
#     mss <- ifelse(attr(z$terms, "intercept"), sum((f - mean(f))^2), sum(f^2)) 
#     rss <- sum(r^2)
#     
#     resvar <- rss/rdf
#     R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
#     se <- sqrt(diag(R) * resvar)
#     est <- z$coefficients[Qr$pivot[p1]]
#     sigma <- sqrt(resvar)
#     list(est=est, vcov=(sigma^2 * R))
# }
