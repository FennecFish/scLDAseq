#Optimization for the Global Covariance Matrix
opt.sigma <- function(nu, omega, lambda, pi, mu, sigprior, samples) {  
    
    sigma <- list()
    for (i in 1:length(unique(samples))) {

        Ni <- which(samples == unique(samples)[i])
        # browser()
        psi <- matrix(rep(pi[i], times = length(Ni)*ncol(lambda)), nrow = length(Ni))
        sigs <- diag(diag(omega)[i], nrow = nrow(mu), ncol = nrow(mu))
        #find the covariance
        if(ncol(mu)==1) {
          browser()
            covariance <- crossprod(sweep(lambda[Ni,], 2, STATS=as.numeric(mu[,Ni] + psi), FUN="-"))
        } else {
            covariance <- crossprod(matrix(lambda[Ni,], nrow = length(Ni))-matrix(t(mu)[Ni,], nrow = length(Ni)) - psi) 
        }
        sigma[[i]] <- (covariance + nu + sigs)/length(Ni) #add to estimation variance
    }
    sigma <- Reduce("+", sigma)
    sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
    return(sigma)
}


