

#################################################
#################################################
#################################################
#################################################
## simulate_rvar_data_1: 
##    For an R-VAR model $\{Phi_k\}_{k=0}^p$ and the individual
##    covariate vectors Y, we generate R-VAR data. 
##
##  INPUTS
##    Phi_list  :  
##    X_covar   : total number of variables in time series.
##    Y         : numeric N x p matrix of covariates.
##    N         : number of individuals sampled.
##    T         : number of time measurements per individual.
##
##  OUTPUT
##    Y         : matrix of covariates for each individual.
##    X_list    : list of length N containing d x T matrices.
##    B_list    : list of length N with dxd matrices with 
##                  individual time series parameters.
##    B_decmp   : list with the decomposed parameters
##                  of a time series.
##    x_cov     : covariate matrix for the error.
##    d         : time series dimension.
##    p         : covariate dimension.
##
simulate_rvar_data_1 <- function(Phi_list, X_covar, Y, N, T) {
  
  ## Simulate y values.
  d <- ncol(Phi_list$Phi_c)
  p <- ncol(Y)
  
  B_list <- list()
  for (n_ind in 1:nrow(Y)) {

    y_ind <- Y[n_ind, ]
    phi_ind_list <- Phi_list$Phi_i
    
    for (l in 1:p) { phi_ind_list[[l]] <- y_ind[l] * phi_ind_list[[l]] }
    
    B_list[[n_ind]] <- Phi_list$Phi_c + Reduce("+", phi_ind_list)
    
  }
  
  X_list <- list()
  
  for (n_ind in 1:N) {
    
    X <- matrix(0, ncol = T + 400, nrow = d)
    E <- t(rmvnorm(T + 400, sigma = X_covar))
    
    for (i in 2:(T + 400)) {
      X[, i] <- B_list[[n_ind]] %*% X[, i-1] + E[, i]
    }
    
    X_list[[n_ind]] <- X[, -(1:400)]
    rm(X)
  }
  
  output <- list(
    Y       = Y,
    X_list  = X_list,
    B_list  = B_list,
    B_decmp = Phi_list,
    x_covar = X_covar,
    d       = d,
    p       = p)
    
  return(output)

}
