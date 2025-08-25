if (!require(mvtnorm)) {
  install.packages("mvtnorm")
  library(mvtnorm)
}
if (!require(plot.matrix)) {
  install.packages("plot.matrix")
  library(plot.matrix)
}
if (!require(tidyverse)) {
  install.packages("tidyverse")
  library(tidyverse)
}
if (!require(magrittr)) {
  install.packages("magrittr")
  library(magrittr)
}
if (!require(glmnet)) {
  install.packages("glmnet")
  library(glmnet)
}



######################################################################
######################################################################
## Generates the parameters B0,B1 of a RVAR model where 
##  there is lag-1 dependency: q = 1.
##
##  INPUTS:
## d          : dimension of the time series Xt in R^d.
## p          : dimension of subject covariates Y in R^p.
## prob_phi0  : connection prob for "joint model" phi0 in R^(d x d)
## prob_phip  : connection prob for "covariate effect" phip in R^(d x d)
## phi0_min, phi0_max : magnitude of non-zero entries in phi0.
## phip_min, phip_max : magnitude of non-zero entries in phip.
## vmin, vmax : magnitude of variances in cov-mat of X-error.
##
##  OUTPUTS:
## phi0       : d x d matrix of joint model.
## phip_list  : list of d x d matrices corresponding to covariate effects.
## p          : number of covariates Y.
## covD       : d x d matrix for error covariance.
##
generate_rvar_pars <- function(d, p, prob_phi0, prob_phip,
                               phi0_min, phi0_max,
                               phip_min, phip_max,
                               vmin, vmax) {
  
  sort <- sample(1:d, d, replace = FALSE)
  
  n_entries <- d * (d-1) / 2
  phi0 <- matrix(0, d, d)
  phi0[upper.tri(phi0)] <- rbinom(n_entries, 1, prob_phi0) * 
    runif(n_entries, phi0_min, phi0_max) *
    sample(c(-1,1), n_entries, TRUE)
  phi0 <- phi0[sort, sort]
  
  phip_list <- list()
  for (ind in 1:p) {
    phip_temp <- matrix(0, d, d)
    phip_temp[upper.tri(phip_temp)] <- rbinom(n_entries, 1, prob_phip) * 
      runif(n_entries, phip_min, phip_max) *
      sample(c(-1,1), n_entries, TRUE)
    
    phip_list[[ind]] <- phip_temp[sort, sort]
  }
  
  covD <- diag(runif(d, vmin, vmax))
  
  output <- list(phi0 = phi0, phip_list = phip_list, p = p, covD = covD)
  return(output)
}


#########################
## EXAMPLES:
#########################
example <- FALSE
if (example) {
  
  
  d         <- 10
  p         <- 2
  prob_phi0 <- 0.3
  prob_phip <- 0.1
  min0      <- 0.3
  max0      <- 0.5
  minp      <- 0.2
  maxp      <- 0.3
  vmin      <- 0.3
  vmax      <- 0.5
  
  output <- generate_rvar_pars(d, p, 
                               prob_phi0, prob_phip, 
                               min0, max0, minp, maxp,
                               vmin, vmax)
  
  par(mfrow = c(3,1), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0)
  plot(output$phip_list[[1]])
  plot(output$phip_list[[2]])
  
  
}


######################################################################
######################################################################
## simulate_exogenous_vars:
##  Function to generate the matrix of time-invariant exogenous 
##  covariates X. 
##
##  INPUTS: 
##    p       : number of exogenous time-invariant covariates.
##    N       : number of subjects.
##    type    : if "unif", it generates independent uniform values in
##                $[-u_max.-u_min]\cup[u_min,u_max]$, or $[u_min,u_max]$,
##                depending on the parameter "signs." If "gaussian", it
##                takes a normal (or absolute value of normal) with sd
##                u_sd. 
##    u_min   : used if type = "unif". 
##    u_max   : used if type = "unif".  
##    u_sd    : used if type = "gaussian". 
##    signed  : whether we want only positive effects or signed effects.
##    nz_x_prob : probability of non-zero entries.
##
##  OUTPUT:
##    X       : N x p matrix of numeric exogenous time invariant covariates.
##
simulate_exogenous_vars <- function(p, N, type,
                                    u_min = 0.5, u_max = 1, g_sd = 0.25,
                                    signed = FALSE, nz_x_prob) {
  
  ## Generate the signed 
  signs <- NULL
  nz_x_probs <- NULL
  if (signed) {
    signs <- c(0, -1, 1)
    nz_x_probs  <- c(1 - nz_x_prob, nz_x_prob / 2, nz_x_prob / 2)
  } else {
    signs <- c(0, 1)
    nz_x_probs  <- c(1 - nz_x_prob, nz_x_prob)
  }
  mat_signs <- sample(signs, p * N, prob = nz_x_probs, TRUE)
  

  ## Generate entry magnitudes:
  mat_entries <- NULL
  if (type == "unif") {
    mat_entries <- runif(p * N, u_min, u_max)
  } else mat_entries <- abs(rnorm(p * N, mean = 0, sd = g_sd))

  X <- matrix(mat_signs * mat_entries, ncol = p, nrow = N)

  return(X)
}




######################################################################
######################################################################
## simulate_rvar1:
##  Simulates data according to the R-VAR model:
##    Yk_t = [Phi0 + Y_k1 Phi1 + ... + Y_kp Phip] Yk_{t-1} + Ek_{t-1}.
##
##  INPUTS:
##    rvar_pars1 : parameters for the R-VAR model generated from the
##                  function "generate_rvar_pars."
##    X          : Matrix of covariates 
##    N          : Number of subjects.
##    Ti         : Vector of time points for each subject.
##
##  OUTPUTS:
##    X           : Nxp matrix of time-invariant exogenous covariates.
##    Y_list      : list of length N containing Ti[k] x d matrices.
##    B_list      : list of length N with dxd matrices with 
##                  individual time series parameters.
##    B_decmp     : list with the decomposed parameters
##                    of a time series.
##    err_covmat  : covariate matrix for the time series error.
##    d           : time series dimension.
##    p           : covariate dimension.
##
##
simulate_rvar1 <- function(rvar_pars1, X, N, Ti) {
  
  ## Simulate y values.
  d <- ncol(rvar_pars1$phi0)
  p <- ncol(X)
  
  ## Calculate subject specific parameters:
  B_list <- list()
  for (k_ind in 1:nrow(X)) {
    xk_ind <- X[k_ind,]
    phip_list_mod <- rvar_pars1$phip_list
    
    for (l in 1:p) phip_list_mod[[l]] <- xk_ind[l] * phip_list_mod[[l]] 
    
    B_list[[k_ind]] <- rvar_pars1$phi0 + Reduce("+", phip_list_mod)
    
  }
  
  ## Generate subject-specific time series data:
  Y_list <- list()
  for (k_ind in 1:N) {
    
    Y <- matrix(0, ncol = Ti[k_ind] + 400, nrow = d)
    
    E <- t(rmvnorm(Ti[k_ind] + 400, sigma = rvar_pars1$covD))
    
    for (i in 2:(Ti[k_ind] + 400)) {
      Y[, i] <- B_list[[k_ind]] %*% Y[, i - 1] + E[, i]
    }
    
    Y_list[[k_ind]] <- t(Y[, -(1:400)])
    rm(Y)
  }
  
  ## Gather results and return:
  output <- list(
    X           = X,
    Y_list      = Y_list,
    B_list      = B_list,
    B_dcmp      = rvar_pars1[1:2],
    err_covmat  = rvar_pars1$covD,
    d           = d,
    p           = p)
  
  return(output)
  
}

#########################
## EXAMPLES:
#########################

example <- FALSE
if (example) {
  
  #########################
  ## Generate parameters:
  #########################
  library(plot.matrix)
  library(tidyverse)
  library(magrittr)
  
  set.seed(20)
  ## R-VAR parameters:
  d         <- 5
  p         <- 2
  prob_phi0 <- 0.35
  prob_phip <- 0.15
  min0      <- 0.3
  max0      <- 0.5
  minp      <- 0.2
  maxp      <- 0.3
  vmin      <- 0.3
  vmax      <- 0.5
  
  ## Exogenous parameters:
  type      <- "unif"
  u_min     <- 0.5
  u_max     <- 1
  signed    <- TRUE
  nz_x_prob   <- 0.7

  output <- generate_rvar_pars(d, p, 
                               prob_phi0, prob_phip, 
                               min0, max0, minp, maxp,
                               vmin, vmax)
  
  par(mfrow = c(3,1), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0)
  plot(output$phip_list[[1]])
  plot(output$phip_list[[2]])
  
  #########################
  ## Generate Data:
  #########################
  N <- 50
  X         <- simulate_exogenous_vars(p = p, N = N, type = "unif",
                                       u_min = u_min, u_max = u_max,
                                       signed = signed, nz_x_prob = nz_x_prob)
  sims_data <- simulate_rvar1(output, X = X, N = N, Ti = rep(100, N))
  
  
  lapply(sims_data$Y_list, dim)
  
  lapply(sims_data$B_list, dim)
  lapply(sims_data$B_list, function(x) {sum(x != 0)})
  
  
  #########################
  ## Visualizing parameters:
  #########################
  
  par(mfrow = c(3,3))
  col_lims <- seq(-0.6, 0.6, length.out = 10)
  plot(sims_data$B_dcmp$phi0, main = "Joint Effect", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[1]], main = "Individual Effect Y1", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[2]], main = "Individual Effect Y2", breaks = col_lims)
  
  xrange <- c(-max(abs(sims_data$X)) - 0.1, 
              max(abs(sims_data$X))+ 0.1)
  
  plot(sims_data$X[,1], sims_data$X[,2], 
       xlab = "X1", ylab = "X2", main = "Exogenous Covariates X",
       xlim = xrange, ylim = xrange, col = rep(c("red","black"), c(5, N-5)))
  text(x = sims_data$Y[,1], y = 0.3 + sims_data$Y[,2],  # Fine-tune the position
       label = c(1:5, rep("", N-5)), col = rep(c("red","black"), c(5, N-5))) 
  
  plot(sims_data$B_list[[1]], main = "Sample 1", breaks = col_lims)
  plot(sims_data$B_list[[2]], main = "Sample 2", breaks = col_lims)
  plot(sims_data$B_list[[3]], main = "Sample 3", breaks = col_lims)
  plot(sims_data$B_list[[4]], main = "Sample 4", breaks = col_lims)
  plot(sims_data$B_list[[5]], main = "Sample 5", breaks = col_lims)
  
  
  
  
  #########################
  ## Visualizing data:
  #########################
  
  sims_data$Y_mat <- lapply(
    seq_along(sims_data$Y_list), 
    function(k, data) {
      x <- t(data[[k]])
      colnames(x) <- paste0("t", 1:ncol(x))
      x <- as_tibble(x) %>%
        mutate(subject = k, 
               var = 1:nrow(x), 
               .before = 1)
      return(x)},
    data = sims_data$Y_list) %>%
    
    {Reduce(rbind, .)}
  
  sims_data$Y_mat %>%
    
    as_tibble() %>%
    
    filter(subject < 11) %>%
    
    pivot_longer(cols = t1:t100,
                 names_to = "time",
                 names_prefix = "t", 
                 values_to = "value") %>%
    
    mutate(time = as.numeric(time),
           var = factor(var),
           subject = factor(subject)) %>%
    
    ggplot(aes(x= time, y = value)) +
    geom_line(aes(col = var)) +
    facet_grid(subject ~ var)
  
  
  
  
}

