#Optimization for the Global Covariance Matrix
opt.sigma <- function(nu, omega, lambda, pi, mu, sigprior, samples) {  
    
    sigma <- list()
    for (i in 1:length(unique(samples))) {
        Ni <- which(samples == unique(samples)[i])
        psi <- rep(pi[i], times = length(Ni))
        sigs <- diag(diag(omega)[i], nrow = nrow(mu), ncol = nrow(mu))
        #find the covariance
        if(ncol(mu)==1) {
            covariance <- crossprod(sweep(lambda[Ni,], 2, STATS=as.numeric(mu[,Ni] - psi), FUN="-"))
        } else {
            covariance <- crossprod(lambda[Ni,]-t(mu[,Ni]) - psi) 
        }
        sigma[[i]] <- (covariance + nu + sigs)/length(Ni) #add to estimation variance
    }
    sigma <- Reduce("+", sigma)
    sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
    return(sigma)
}


