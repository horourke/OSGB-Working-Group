#################################################
#################################################
#################################################
##
## In the following document, we will
## define a method for generating positive
## definite matrices with hubs, which is 
## based on first generating the eigen-
## structure.
##
#################################################
#################################################



#################################################
#################################################
## generate_rvar_model: 
##    For an R-VAR model of the form:
##    $$ X(t) = [\Phi_0  + \sum_{i=1}^p y_i\Phi_k]  X(t-1) + E(t),$$
##    we generate the matrices $\{Phi_k\}_{k=0}^p$.
##
##  INPUTS
##    d         : total number of variables in time series.
##    p         : number of individual-level covariates.
##    prob_tot  : total proportion of non-zero entries.
##    prob_c    : total proportion of common entries.
##    entry_min : minimum value of non-zero entries.
##    entry_max : maximum value of non-zero entries.
##
##  OUTPUT
##    Phi_c   : d x d matrix corresponding to common structure.
##    Phi_i   : list containing p matrices of dimension d x d.
##
generate_rvar_matrices_1 <- function(
  d, p, prob_tot, prob_c, 
  entry_min, entry_max, signed = FALSE) {
    
  entryct_com  <- ceiling(d * d * prob_tot * prob_c)
  entryct_ind  <- ceiling(d * d * prob_tot * (1 - prob_c) / p)
  entryct_null <- d * d - entryct_com - entryct_ind * p
    
  entryvec_reps <- rep(
    c(-1,          0,            1:p),
    c(entryct_com, entryct_null, rep(entryct_ind, p)))

  ind_mat <- matrix(0, d, d)
  ind_mat[] <- sample(entryvec_reps, d*d, replace = FALSE)

  signs <- NULL
  if (signed) {
    signs <- c(-1, 1)
  } else {
    signs <- c(1)
  } 

  mat_list <- lapply(
    c(-1, 1:p),
    function(ind, ind_mat, signs) {
      nmat <- matrix(
        runif(d*d, min = entry_min, max = entry_max) * 
          sample(signs, d*d, replace = TRUE),
        d, d)
      nmat 
      nmat <- nmat * (ind_mat == ind)

      return(nmat)
      
    },
    ind_mat, signs)

  output <- list(
    Phi_c = mat_list[[1]],
    Phi_i = mat_list[-1]
  )
  return(output)
}

example = FALSE
if (example) {

  library(plot.matrix)

  d = 5
  p = 2
  prob_tot = 0.5
  prob_c   = 0.5

  entry_min = 0.1
  entry_max = 0.9

  mats <- generate_rvar_matrices_1(d,p,prob_tot,prob_c,entry_min,entry_max)

  par(mfrow = c(3,1))
  plot(mats$Phi_c != 0)
  plot(mats$Phi_i[[1]] != 0)
  plot(mats$Phi_i[[2]] != 0)
  dev.off()
}
rm(example)


#################################################
#################################################
#################################################
#################################################
## generate_rvar_models_1: 
##    For an R-VAR model of the form:
##    $$ X(t) = [\Phi_0  + \sum_{i=1}^p y_i\Phi_k]  X(t-1) + E(t),$$
##    we generate the matrices $\{Phi_k\}_{k=0}^p$ and the individual
##    covariate vectors Yk. 
##    In the generating process, we verify that the generated sample
##    satisfies the condition for causal stationary solution:
##      det(PHI(z)) != 0 for  |z| <= 1 for z in C,
##      where PHI(z) = I - Phi_1 * z - Phi_2 * z^2 - ... - Phi_q * z^q. 
##
##
##  INPUTS
##    d         : total number of variables in time series.
##    p         : number of individual-level covariates.
##    prob_tot  : total proportion of non-zero entries.
##    prob_c    : total proportion of common entries.
##    entry_min : minimum value of non-zero entries.
##    entry_max : maximum value of non-zero entries.
##    N         : number of individuals sampled.
##
##  OUTPUT
##    Phi_list  : list containing output of generate_rvar_matrices_1()
##    Y         : matrix of covariates for each individual.
##
generate_rvar_models_1 <- function(
  d, p, prob_tot, prob_c,
  entry_min, entry_max,  signed = FALSE,
  N, max_iter) {

  count <- 0 
  is.unstable <- TRUE
  while(count < max_iter & is.unstable) {
    
    Phi_list  <- generate_rvar_matrices_1(d, p, prob_tot, prob_c, entry_min, entry_max, signed) 
    Y         <- rmvnorm(N, mean = rep(0,p)) / 2

    unstable_count <- sum(verify_stab_1(Phi_list, Y))

    is.unstable <- (unstable_count != 0)
    count       <- count + 1

  }

  if(is.unstable) {
    stop("Error: no stable sampling obtained before the max.iteration")
  }

  output <- list(
    Phi_list = Phi_list,
    Y = Y)
  return(output)

}


#################################################
## .verify_stab_1: 
##    Auxiliary function. 
##    Verifies if, for a given set of covariates "y"
##    and a list Phi_list, it verifies if the lag-1
##    time series model:
##
##    $$ \Phi = \Phi_0 + \sum_{k=1}^p y_k \cdot \Phi_k$$
##
##    satisfies the eigenvalue stability condition.
##
##
##  INPUTS
##    Phi_list  : list resulting as output of function
##                  generate_rvar_matrices_1().
##    y         : numeric vector of length p of covariates.
##
##  OUTPUT
##    unstab_eig_count : number of eigenvalues of $\Phi$ with 
##                        complex module greater or equal than 1.
##
.verify_stab_1 <- function(Phi_list, y) {
  
  p         <- length(y) 
  mat_aggr  <- 1:p %>% 
  
    lapply(
      function(k, mat_list, y) {
        return(y[k] * mat_list[[k]])
      }, 
      mat_list = Phi_list$Phi_i, y = y) %>% 
    
    {Reduce('+', .)} %>%

    {(.) + Phi_list$Phi_c}

  #print(mat_aggr)

  eig_aggr  <- eigen(mat_aggr)$values
  unstab_eig <- ifelse(Mod(eig_aggr) >= 1, TRUE, FALSE)

  return(sum(unstab_eig))

}


#################################################
## verify_stab_1: 
##    Auxiliary function. 
##    Verifies if, for a given set of covariates "y"
##    and a list Phi_list, it verifies if the lag-1
##    time series model:
##
##    $$ \Phi = \Phi_0 + \sum_{k=1}^p y_k \cdot \Phi_k$$
##
##    satisfies the eigenvalue stability condition.
##
##
##  INPUTS
##    Phi_list  : list resulting as output of function
##                  generate_rvar_matrices_1().
##    Y         : numeric N x p matrix of covariates.
##
##  OUTPUT
##    is_unstable : logical vector N with TRUE entries for
##                    observations that are unstable.
##
verify_stab_1 <- function(Phi_list, Y) {

  N <- nrow(Y)
  is_unstable <- rep(FALSE, N)
  for(obs_ind in 1:N) {
    y <- Y[obs_ind,]
    count_unstable <- sum(.verify_stab_1(Phi_list, y))
    is_unstable[obs_ind] <- (count_unstable > 0)
  }

  return(is_unstable)
}






example = FALSE
if (example) {

  library(plot.matrix)
  library(mvtnorm)
  library(magrittr)

  N <- 10
  d <- 10
  p <- 2

  prob_tot <- 0.2
  prob_c   <- 0.5

  entry_min <- 0.1
  entry_max <- 0.5

  mats <- generate_rvar_matrices_1(d,p,prob_tot,prob_c,entry_min,entry_max, FALSE)
  mats <- generate_rvar_matrices_1(d,p,prob_tot,prob_c,entry_min,entry_max, TRUE)
  Y    <- rmvnorm(N, mean = rep(0,p)) / 2

  par(mfrow = c(3,1))
  plot(mats$Phi_c != 0)
  plot(mats$Phi_i[[1]] != 0)
  plot(mats$Phi_i[[2]] != 0)
  dev.off()

  par(mfrow = c(3,1))
  plot(mats$Phi_c)
  plot(mats$Phi_i[[1]])
  plot(mats$Phi_i[[2]])
  dev.off()

  .verify_stab_1(mats, Y[1,])
  verify_stab_1(mats, Y)



  #################################################
  #################################################

  max_iter <- 1000

  model <- generate_rvar_models_1(
    d, p, prob_tot, prob_c,
    entry_min, entry_max, signed = FALSE,
    N, max_iter)

  par(mfrow = c(3,1))
  plot(model$Phi_list$Phi_c != 0)
  plot(model$Phi_list$Phi_i[[1]] != 0)
  plot(model$Phi_list$Phi_i[[2]] != 0)
  dev.off()

  par(mfrow = c(3,1))
  plot(model$Phi_list$Phi_c)
  plot(model$Phi_list$Phi_i[[1]])
  plot(model$Phi_list$Phi_i[[2]])
  dev.off()


  #################################################
  #################################################
  max_iter <- 1000

  model <- generate_rvar_models_1(
    d, p, prob_tot, prob_c,
    entry_min, entry_max, signed = TRUE,
    N, max_iter)

  par(mfrow = c(3,1))
  plot(model$Phi_list$Phi_c != 0)
  plot(model$Phi_list$Phi_i[[1]] != 0)
  plot(model$Phi_list$Phi_i[[2]] != 0)
  dev.off()

  par(mfrow = c(3,1))
  plot(model$Phi_list$Phi_c)
  plot(model$Phi_list$Phi_i[[1]])
  plot(model$Phi_list$Phi_i[[2]])
  dev.off()



}
rm(example)

