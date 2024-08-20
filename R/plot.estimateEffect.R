plot.estimateEffect <- function(x, covariate, model=NULL,
                                topics=x$topics,
                                method=c("pointestimate", "difference","continuous"),
                                cov.value1=NULL, cov.value2=NULL,
                                moderator=NULL, moderator.value=NULL,
                                npoints=100, nsims=100, ci.level=.95,
                                xlim=NULL, ylim=NULL, xlab="",ylab=NULL,
                                main="", printlegend=T,
                                labeltype="numbers", n=7, frexw=.5,
                                add=F, linecol=NULL, width=25,
                                verbose.labels=T, family=NULL,
                                custom.labels=NULL, omit.plot=FALSE,...){
  
  method <- match.arg(method)
  if(method=="difference" && (is.null(cov.value1) | is.null(cov.value2))) {
    stop("For method='difference' both cov.value1 and cov.value2 must be specified.")
  }
  #Produce cdata (data in original form) and
  #cmatrix (data in design matrix form)
  cthis <- produce_cmatrix(prep=x, covariate=covariate, method=method,
                           cov.value1=cov.value1,
                           cov.value2=cov.value2, npoints=npoints,
                           moderator=moderator, moderator.value=moderator.value)
  cdata <- cthis$cdata
  cmat <- cthis$cmatrix

  #Simulate betas
  simbetas <- simBetas(x$parameters, nsims=nsims)

  #Find offset for confidence level
  offset <- (1-ci.level)/2

  
  #Plot for each method
  if(method=="continuous"){
    toreturn <- plotContinuous(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                   offset=offset,xlab=xlab, ylab=ylab, main=main,
                   xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                   labeltype=labeltype,n=n,custom.labels=custom.labels,model=model,frexw=frexw,printlegend=printlegend,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
  if(method=="pointestimate"){
    toreturn <- plotPointEstimate(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                      offset=offset,xlab=xlab, ylab=ylab, main=main,
                      xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                      labeltype=labeltype,n=n,
                                  custom.labels=custom.labels,model=model,frexw=frexw,width=width,
                                  verbose.labels=verbose.labels,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
  if(method=="difference"){
    if(missing(cov.value1)) stop("Missing a value for cov.value1. See documentation.")
    if(missing(cov.value2)) stop("Missing a value for cov.value2. See documentation.")
    toreturn <- plotDifference(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                   offset=offset,xlab=xlab, ylab=ylab, main=main,
                   xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                   labeltype=labeltype,n=n,
                   custom.labels=custom.labels,model=model,frexw=frexw,width=width,
                   cov.value1=cov.value1,
                               cov.value2=cov.value2,verbose.labels=verbose.labels,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
}