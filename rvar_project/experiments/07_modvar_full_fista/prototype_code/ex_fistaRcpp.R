

install.packages("RcppArmadillo")
library("RcppArmadillo")
Rcpp::sourceCpp("C:/Users/joses/Dropbox/UC-Riverside Documents/Research/OSGB-Working-Group/rvar_project/experiments/07_modvar_full_fista/fista.cpp")

set.seed(1)

n <- 200
p <- 100

A <- matrix(rnorm(n*p), n, p)
beta_true <- c(rep(2, 5), rep(0, p-5))
b <- A %*% beta_true + rnorm(n)

# Lipschitz constant = largest eigenvalue of A'A
L <- max(eigen(t(A) %*% A, symmetric = TRUE, only.values = TRUE)$values)

beta_hat <- fista_lasso(A, b, lambda = 30, L = L)

plot(beta_true, beta_hat)




