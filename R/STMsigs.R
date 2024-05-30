#Optimization for the patient hereogeneity Covariance Matrix
opt.sigs <- function(pi, omega, samples) {  
    sigs <- matrix(0, nrow = nrow(omega), ncol = ncol(omega))
    psi <- numeric(length(unique(samples)))
    for (i in 1:length(unique(samples))) {
        Ni <- which(samples == unique(samples)[i])
        psi[i] <- mean(colMeans(pi[Ni,]))
        sigs[i,i] <- psi[i]^2 + diag(omega)[i]
    }
    return(list(sigs = sigs, psi = psi))
}


