#plotContinuous: function to plot stm output as a function of a
#continuous variable
#prep: output of plot.estimateEffect
#covariate: covariate of interest
#topics: topics of interest
#cdata: output of produce_cmatrix, data for simulation
#cmat: cmatrix output of procude_cmatrix, model matrix for simulation
#simbetas: simulated betas from stimBetas
#offset: offset for quantiles from ci.level
#linecol: colors for lines
#add: whether to plot new graph or add to old one
#labeltype: type of label
#n: number of words to print in label
#custom.labels: custom labels for labeltype="custom"
#model: stm model (used for labeling)
#rest of arguments are graphical, see par

plotContinuous <- function(prep,covariate,topics, cdata, cmat, simbetas, offset,xlab=NULL,
                           ylab=NULL, main=NULL, xlim=NULL, ylim=NULL,
                           linecol=NULL, add=F, labeltype,n=7, custom.labels=NULL,model=NULL,frexw=.5,printlegend=T,
                           omit.plot=FALSE,
                           ...){
  #What are the unique values of the covariate we are going to plot over?
  uvals <- cdata[,covariate]
  
  #For each topic, 1. Simulate values, 2. Find means and cis
  means = list()
  cis = list()
  for(i in 1:length(topics)){
    #Simulate values
    sims <- cmat%*%t(simbetas[[which(prep$topics==topics[i])]])
    #Find means and cis
    means[[i]] <- rowMeans(sims)
    cis[[i]] = apply(sims,1, function(x) quantile(x, c(offset,1-offset)))
  }

  #Get the plot set up
  if(!omit.plot) {
    if (is.null(xlim)) xlim <- c(min(uvals), max(uvals))
    if (is.null(ylim)) ylim <- c(min(unlist(cis), na.rm=T),
                                 max(unlist(cis), na.rm=T))
    if (is.null(ylab)) ylab <- "Expected Topic Proportion"
    if (add==F) plot(0, 0,col="white",xlim=xlim, ylim=ylim, main=main,
                             xlab=xlab, ylab=ylab,  ...)
    if (is.null(linecol)) cols = grDevices::rainbow(length(topics))
    if (!is.null(linecol)) cols=linecol
  
    #Plot everything
    for(i in 1:length(topics)){
      lines(uvals, means[[i]], col=cols[i])
      lines(uvals, cis[[i]][1,], col=cols[i], lty=2)
      lines(uvals, cis[[i]][2,], col=cols[i], lty=2)
    }
  }
  #Create legend
  if(printlegend==T){
    labels = createLabels(labeltype=labeltype, covariate=covariate, method="continuous",
      cdata=cdata, ref=NULL, alt=NULL,model=model,n=n,
      topics=topics,custom.labels=custom.labels, frexw=frexw)
    if(!omit.plot) legend(xlim[1], ylim[2], labels, cols)
    return(invisible(list(x=uvals, topics=topics,means=means, ci=cis, labels=labels)))
  }else{
    return(invisible(list(x=uvals, topics=topics,means=means, ci=cis)))
  }
  
}


#plotPointEstimate: function to plot stm output for unique values of a
#variable
#prep: output of plot.estimateEffect
#covariate: covariate of interest
#topics: topics of interest
#cdata: output of produce_cmatrix, data for simulation
#cmat: cmatrix output of procude_cmatrix, model matrix for simulation
#simbetas: simulated betas from stimBetas
#offset: offset for quantiles from ci.level
#labeltype = type of label
#n: number of words to print in labels
#custom.labels: labels if labeltype="custom"
#model: model output
#frexw: frex weight
#width: width adjustment for labels
#rest of arguments are graphical, see par

plotPointEstimate <- function(prep,covariate,topics, cdata, cmat, simbetas, offset,xlab=NULL,
                           ylab=NULL, main=NULL, xlim=NULL, ylim=NULL,
                           linecol=NULL, add=F, labeltype,n=7, custom.labels=NULL,model=NULL,frexw=.5,width=25,verbose.labels=T,
                           omit.plot=FALSE, ...){
  #What are the unique values of the covariate we are going to plot over?
  uvals <- cdata[,covariate]
  
  #For each topic, 1. Simulate values, 2. Find means and cis
  means = list()
  cis = list()
  for(i in 1:length(topics)){
    #Simulate values
    sims <- cmat%*%t(simbetas[[which(prep$topics==topics[i])]])
    #Find means and cis
    means[[i]] <- rowMeans(sims)
    cis[[i]] = apply(sims,1, function(x) quantile(x, c(offset,1-offset)))
  }

  
  if(!omit.plot) {
    #Get the plot set up
    if (is.null(ylim)) ylim <- c(0, length(topics)*length(uvals)+1)
    if (is.null(xlim)) xlim <- c(min(unlist(cis), na.rm=T),
                                 max(unlist(cis), na.rm=T))
    if(is.null(xlab)) xlab <- "Estimated Topic Proportions"
    if (is.null(ylab)) ylab <- ""
    if (add==F) plot(0, 0,col="white",xlim=xlim, ylim=ylim, main=main,
                             xlab=xlab, ylab=ylab,yaxt="n",...)
  
    
    #Add a line for the 0 x-axis
    lines(c(0,0), c(0, length(topics)*length(uvals)+2), lty=2)
  }
  #Create the labels:
  labels = createLabels(labeltype=labeltype, covariate=covariate, method="pointestimate",
      cdata=cdata, ref=NULL, alt=NULL,model=model,n=n,
      topics=topics,custom.labels=custom.labels, frexw=frexw, verbose.labels=verbose.labels)

  if(!omit.plot) {
    #Plot everything
    #Start at the top, move down
    it <- length(topics)*length(uvals)
    lab <- 1
    for(i in 1:length(uvals)){
      for(j in 1:length(topics)){
        points(means[[j]][i], it, pch=16)
        lines(c(cis[[j]][1,i],cis[[j]][2,i]),c(it,it))
        axis(2,at=it, labels=stringr::str_wrap(labels[lab],width=width),las=1, cex=.25, tick=F, pos=cis[[j]][1,i])
        it = it-1
        lab = lab+1
      }
    }
  }
  return(invisible(list(uvals=uvals, topics=topics, means=means, cis=cis, labels=labels)))
}

#plotDifference: function to plot stm output for difference between
#two unique values of a variable
#prep: output of plot.estimateEffect
#covariate: covariate of interest
#topics: topics of interest
#cdata: output of produce_cmatrix, data for simulation
#cmat: cmatrix output of procude_cmatrix, model matrix for simulation
#simbetas: simulated betas from stimBetas
#offset: offset for quantiles from ci.level
#rest of arguments are graphical, see par

plotDifference <- function(prep,covariate,topics, cdata, cmat, simbetas, offset,xlab=NULL,
                           ylab=NULL, main=NULL, xlim=NULL, ylim=NULL,
                           ref=NULL, alt=NULL,
                           linecol=NULL, add=F,labeltype,n=7, custom.labels=NULL, printlegend = T,
                           model=NULL,frexw=.5,width=25,verbose.labels=T,
                           omit.plot=FALSE,...){
  #What are the unique values of the covariate we are going to plot over?
  uvals <- cdata[,covariate]
  #For each topic, 1. Simulate values, 2. Find means and cis
  means = list()
  cis = list()
  for(i in 1:length(topics)){
    #Simulate values
    sims <- cmat%*%t(simbetas[[which(prep$topics==topics[i])]])
    #Take difference
    diff <- sims[2,]-sims[1,]
    #Find means and cis
    means[[i]] <- mean(diff)
    cis[[i]] = quantile(diff, c(offset,1-offset))
  }
 
  if(!omit.plot) {
    #Get the plot set up
    if (is.null(ylim)) ylim <- c(0, length(topics)+1)
    if (is.null(xlim)) xlim <- c(min(unlist(cis), na.rm=T),
                                 max(unlist(cis), na.rm=T))
    if(is.null(xlab)) xlab <- "Estimated Topic Proportions"
    if (is.null(ylab)) ylab <- ""
    if (add==F) plot(0, 0,col="white",xlim=xlim, ylim=ylim, main=main,
                             xlab=xlab, ylab=ylab,yaxt="n",...)
    
    
    #Add a line for the 0 x-axis
    lines(c(0,0), c(0, length(topics)+2), lty=2)
  }
  #Create labels
    labels = createLabels(labeltype=labeltype, covariate=covariate, method="difference",
      cdata=cdata, ref=ref, alt=alt,model=model,n=n,
      topics=topics,custom.labels=custom.labels, frexw=frexw, verbose.labels=verbose.labels)

  #Plot everything
  #Start at the top, move down
  if(!omit.plot) {
    it <- length(topics)
    for(i in 1:length(topics)){
        if(add){
            # Shift the line for better visualization
            shift_amount <- 0.1  
            shifted_y_values <- it - shift_amount
            points(means[[i]], shifted_y_values, pch = 16, col = linecol)
            lines(c(cis[[i]][1],cis[[i]][2]),c(shifted_y_values,shifted_y_values), col = linecol)
        }else{
            points(means[[i]], it, pch=16, col = linecol)
            lines(c(cis[[i]][1],cis[[i]][2]),c(it,it), col = linecol)
        }
      if(printlegend){
          axis(2,at=it, labels=stringr::str_wrap(labels[i],width=width),las=1, cex=.25, tick=F, pos=cis[[i]][1])
      }
      it = it-1
    }
  }
  return(invisible(list(topics=topics, means=means, cis=cis, labels=labels)))
}

#CreateLabels: function that creates labels for each plotting function.
#labeltype: type of label (numbers, prob, lift, frex, score, custom)
#covariate: covariate of interest
#method: plotting method
#cdata: cdata output from produce_cmatrix
#ref: for method difference, first level of the covariate
#alt: for method difference, second level of the covariate
#model: model output for labeling topics
#n: number of words to print for prob, lift, frex, score
#topics: topics of interest
#custom.labels: user-inputted customlabels
createLabels <- function(labeltype,covariate, method, cdata, ref,
                         alt, model, n, topics, custom.labels, 
                         frexw=.5, verbose.labels=T){
  #Final Output
  labelsout <- NULL

  #For each type of label, create a preliminary label (without
  #covariate specific labeling, just with topic labeling

  #Label type numbers, just topic numbers:
  if(labeltype=="numbers"){
      labelsout=paste("Topic", topics)
    }
  if(labeltype!="numbers" & labeltype!="custom" & is.null(model)){
    stop("Model output is needed in order to use words as labels.  Please enter the model output for the argument 'model'.")
  }

  #Labeltype prob, probability labeling
  if(labeltype=="prob"){
    for(i in 1:length(topics)){
      labelsout[i] <- commas(labelTopics(model, n=n)$prob[topics[i],])
    }
  }
  #Labeltype lift, lift labeling
  if(labeltype=="lift"){
    for(i in 1:length(topics)){
      labelsout[i] <- commas(labelTopics(model, n=n)$lift[topics[i],])
    }
  }
  #Labeltype frex, frex labeling
  if(labeltype=="frex"){
    for(i in 1:length(topics)){
      labelsout[i] <- commas(labelTopics(model, n=n, frexweight=frexw)$frex[topics[i],])
    }
  }
  #Labeltype score, score labeling
  if(labeltype=="score"){
    for(i in 1:length(topics)){
      labelsout[i] <- commas(labelTopics(model, n=n)$score[topics[i],])
    }
  }
  if(labeltype=="custom"){
    labelsout = custom.labels
  }else{
    if(verbose.labels==T){
                                        #For methods point estimate and difference, add covariate
                                        #information
                                        #Method: Point estimate
      if(method=="pointestimate"){
                                        #For each topic, add the
                                       #covariate level
        labels <- labelsout
        labelsout <- NULL
        uvals <- cdata[,covariate]
        it <- 1
        for(i in 1:length(uvals)){
          for(j in 1:length(labels)){
            labelsout[it] <- paste(labels[j], "(Covariate Level: ", uvals[i],
                                   ")", sep="")
            it = it +1
          }
        }        
      }
                                        #Method: Difference
      if(method=="difference"){
        labels <- labelsout
        labelsout <- NULL
        labelsout <-  paste(labels, " (Covariate Level ", ref, " Compared to ",
                            alt, ")", sep="")
      }
    }
  }
  return(labelsout)
}
