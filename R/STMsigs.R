#Optimization for the patient hereogeneity Covariance Matrix
opt.sigs <- function(pi, omega, samples) {  
    
    sigs <- matrix(0, nrow = nrow(omega), ncol = ncol(omega))
    for (i in 1:length(unique(samples))) {
        Ni <- which(samples == unique(samples)[i])
        sigs[i,i] <- pi[i]^2 + diag(omega)[i]
    }
    return(sigs)
}


