library(compositions)

# use multivariate normal
mean <- sample(200:300, 5)

# create sigma
n <- 5
A <- matrix(runif(n^2)*2-1, ncol=n) 
Sigma <- t(A) %*% A

x1 <- mvrnorm(100, mean, Sigma)
x2 <- mvrnorm(100, mean, Sigma)

# use ilr
y1 <- clr(x1)
y2 <- clr(x2)

n1 <- nrow(y1)
n2 <- nrow(y2)

m1 <- colMeans(y1)
m2 <- colMeans(y2)

pcov <- (n1* cov(y1) + n2 *cov(y2))/(n1+n2)
cmu <- (n1 * m1 + n2 * m2)/(n1 + n2)
ccov <- pcov + (n1*n2*tcrossprod(m1-m2))/((n1 + n2)^2)

# test 1
Q1 <- n1*log(det(ccov)/det(cov(y1))) + n2*log(det(ccov)/det(cov(y2))) 
pvalue <- pchisq(Q, df=4*5/2, lower.tail=FALSE)
