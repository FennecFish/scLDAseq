#E-Step for a Document Block
#[a relatively straightforward rewrite of previous
# code with a focus on avoiding unnecessary computation.]

#Input: Documents and Key Global Parameters
#Output: Sufficient Statistics

# Approach:
# First we pre-allocate memory, and precalculate where possible.
# Then for each document we:
#  (1) get document-specific priors, 
#  (2) infer doc parameters, 
#  (3) update global sufficient statistics
# Then the sufficient statistics are returned.

#Let's start by assuming its one beta and we may have arbitrarily subset the number of docs.
estep <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                       beta, lambda.old, mu, sigma,
                       samples, pi.old, sigs,
                       verbose) {
  
  #quickly define useful constants
  V <- ncol(beta$beta[[1]])
  K <- nrow(beta$beta[[1]])
  N <- length(documents)
  A <- length(beta$beta)
  I <- length(unique(samples))
  ctevery <- ifelse(N>100, floor(N/100), 1)
  # browser()
  if(!update.mu) mu.l <- as.numeric(mu)
  
  # 1) Initialize Sufficient Statistics 
  sigma.ss <- diag(0, nrow=(K-1))
  beta.ss <- vector(mode="list", length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  
  # 2) Precalculate common components
    # calculate inverse and entropy of sigmat
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if(inherits(sigobj,"try-error")) {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
      sigmaentropy <- sum(log(diag(sigobj)))
      siginv <- chol2inv(sigobj)
  }
  
  # 3) Document Scheduling
  # For right now we are just doing everything in serial.
  # the challenge with multicore is efficient scheduling while
  # maintaining a small dimension for the sufficient statistics.
  
  
  omega <- matrix(0, nrow = I, ncol = I)
  for (i in 1:I) {
  psi.i <- rep(pi.old[i], ncol(lambda.old)) # repeat pi into a K-1 dimensional vector
  sigs.i <- diag(sigs)[i]
  omega.i <- 1/(sum(diag(siginv)) + 1/sigs.i)
  Ni <- which(samples == unique(samples)[i])
      for(l in Ni) {
          #update components
          doc <- documents[[l]]
          words <- doc[1,]
          aspect <- beta.index[l]
          init <- lambda.old[l,]
          if(update.mu) mu.l <- mu[,l]
          beta.l <- beta$beta[[aspect]][,words,drop=FALSE]
          
          #infer the document
          # browser()
          doc.results <- logisticnormalcpp(eta=init, mu=mu.l, psi = psi.i,
                                           siginv=siginv, 
                                           sigs = sigs.i,
                                           beta=beta.l, 
                                           doc=doc, 
                                           sigmaentropy=sigmaentropy)
          # update sufficient statistics 
          sigma.ss <- sigma.ss + doc.results$eta$nu
          beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
          bound[l] <- doc.results$bound
          lambda[[l]] <- c(doc.results$eta$lambda)
          if(verbose && l%%ctevery==0) cat(".")
      }
  pi[i] <- doc.results$pi
  omega[i,i] <- omega.i
  }
  
  if(verbose) cat("\n") #add a line break for the next message.

  #4) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, 
              lambda=lambda, pi = pi, omega = omega))
}

