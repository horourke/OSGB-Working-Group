
######################################################################
######################################################################
## Generates the parameters B0,B1, Deltas of a Multi+MOD VAR model where 
##  there is lag-1 dependency: q = 1.
## This new version ensures that nomothetic, idiographic and moderator entries
## are disjoint. This way, it is easier to evaluate the performance
## of all components
##
##  INPUTS:
## d          : dimension of the time series Yt in R^d.
## p          : dimension of subject covariates X in R^p.
## prob_phi0  : connection prob for "joint model" phi0 in R^(d x d)
## prob_phip  : connection prob for "moderate effect" phip in R^(d x d)
## prob_delta : connection prob for "idiographic effects" Delta in R^(d x d)
## phi0_min, phi0_max   : magnitude of non-zero entries in phi0.
## phip_min, phip_max   : magnitude of non-zero entries in phip.
## delta_min, delta_max : magnitude of non-zero entries in Delta.
## vmin, vmax : magnitude of variances in cov-mat of X-error.
##
##  OUTPUTS:
## phi0       : d x d matrix of joint model.
## phip_list  : list of d x d matrices corresponding to moderator effects.
## delta_list : list of d x d matrices corresponding to fully idiographic effects.
## p          : number of covariates X.
## covD       : d x d matrix for error covariance.
##
generate_modvar_pars <- function(
  d, p, n, prob_phi0, prob_phip, prob_delta,
  phi0_min, phi0_max,
  phip_min, phip_max,
  delta_min, delta_max,
  vmin, vmax,
  signed = TRUE) {
  
  sort_var <- sample(1:d, d, replace = FALSE)
  
  n_entries <- d * (d - 1) / 2
  n_nonzero0 <- as.integer(d * d * prob_phi0)
  n_nonzerop <- as.integer(d * d * prob_phip)
  n_nonzerod <- as.integer(d * d * prob_delta)
  n_zero <- n_entries - n_nonzero0 - p * n_nonzerop

  ## Determine which entries are non-zero to which effects.
  nzero_inds <- rep(c(-1:p), c(n_nonzero0, n_zero, rep(n_nonzerop, p))) %>%
    sample(n_entries, replace = FALSE)

  ## Create nomothetic effects:
  phi0 <- matrix(0, d, d)
  
  phi0[upper.tri(phi0)] <- (nzero_inds == -1) * 
    runif(n_entries, phi0_min, phi0_max) 
    
  if (signed) {
    phi0[upper.tri(phi0)] <- phi0[upper.tri(phi0)] * 
      sample(c(-1,1), n_entries, TRUE)
  }
  
  phi0 <- phi0[sort_var, sort_var]
  
  ## Create moderator-dependent effects:
  phip_list <- list()
  for (ind in 1:p) {
    phip_temp <- matrix(0, d, d)
  
    phip_temp[upper.tri(phip_temp)] <- (nzero_inds == ind)  * 
      runif(n_entries, phip_min, phip_max) 
      
    if (signed) {
      phip_temp[upper.tri(phip_temp)] <- phip_temp[upper.tri(phip_temp)] * 
        sample(c(-1,1), n_entries, TRUE)
    }
      
    phip_list[[ind]] <- phip_temp[sort_var, sort_var]
  }

  ## Create fully idiographic effects:
  delta_list <- list()
  for (ind in 1:n) {
    zero_loc  <- which(nzero_inds == 0) 
    delta_loc <- sample(zero_loc, n_nonzerod, replace = FALSE) 
    dzero_inds <- (1:n_entries) %in% delta_loc

    delta <- matrix(0, d, d)
    delta[upper.tri(delta)] <- dzero_inds * 
      runif(n_entries, delta_min, delta_max) 

    if (signed) {
      delta[upper.tri(delta)] <- delta[upper.tri(delta)] * 
        sample(c(-1,1), n_entries, TRUE)
    }
    delta_list[[ind]] <- delta[sort_var, sort_var]
  }
      
  covD <- diag(runif(d, vmin, vmax))
  
  output <- list(
    phi0 = phi0, 
    phip_list = phip_list, 
    delta_list = delta_list, 
    p = p, n = n, covD = covD)
  return(output)
}


#########################
## EXAMPLES:
#########################
example <- FALSE
if (example) {
  
  wd <- getwd()
  req_lib_dir <- paste0(wd,"/req_lib")
  .libPaths(req_lib_dir)

  library(plot.matrix)
  library(magrittr)
  library(tidyverse)

  set.seed(10)
  d           <- 10
  p           <- 2
  n           <- 3
  prob_phi0   <- 0.3
  prob_phip   <- 0.05
  prob_delta  <- 0.02
  min0        <- 0.1
  max0        <- 0.2
  minp        <- 0.3
  maxp        <- 0.4
  mind        <- 0.5
  maxd        <- 0.6
  vmin        <- 0.3
  vmax        <- 0.5
  
  output <- generate_modvar_pars(
    d, p, n, 
    prob_phi0, prob_phip, prob_delta,
    min0, max0, 
    minp, maxp, 
    mind, maxd,
    vmin, vmax, signed = FALSE)


  
  pdf("plot.pdf", width = 7, height = 3)
  par(mfrow = c(2,3), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0)
  plot(output$phip_list[[1]])
  plot(output$phip_list[[2]])
  plot(output$delta_list[[1]])
  plot(output$delta_list[[2]])
  plot(output$delta_list[[3]])
  dev.off()
  

  sum(output$phi0 != 0)
  sum(output$phip_list[[1]] != 0)
  sum(output$phip_list[[2]] != 0)
  sum(output$delta_list[[1]] != 0)
  sum(output$delta_list[[2]] != 0)
  sum(output$delta_list[[3]] != 0)
  
  sum(output$phi0 != 0 & output$phip_list[[1]] != 0)
  sum(output$phi0 != 0 & output$phip_list[[2]] != 0)

  sum(output$phi0 != 0 & output$delta_list[[1]] != 0)
  sum(output$phi0 != 0 & output$delta_list[[2]] != 0)
  sum(output$phi0 != 0 & output$delta_list[[3]] != 0)

  sum(output$phip_list[[1]] != 0 & output$phip_list[[2]] != 0)

  sum(output$phip_list[[1]] != 0 & output$delta_list[[1]] != 0)
  sum(output$phip_list[[1]] != 0 & output$delta_list[[2]] != 0)
  sum(output$phip_list[[1]] != 0 & output$delta_list[[3]] != 0)
  sum(output$phip_list[[2]] != 0 & output$delta_list[[1]] != 0)
  sum(output$phip_list[[2]] != 0 & output$delta_list[[2]] != 0)
  sum(output$phip_list[[2]] != 0 & output$delta_list[[3]] != 0)

  ## Verifying the signed parameter:
  set.seed(20)
  d           <- 10
  p           <- 2
  n           <- 3
  prob_phi0   <- 0.3
  prob_phip   <- 0.05
  prob_delta  <- 0.02
  min0        <- 0.3
  max0        <- 0.5
  minp        <- 0.2
  maxp        <- 0.3
  mind        <- 0.5
  maxd        <- 0.6
  vmin        <- 0.3
  vmax        <- 0.5
  
  output <- generate_modvar_pars(
    d, p, n, 
    prob_phi0, prob_phip, prob_delta,
    min0, max0, 
    minp, maxp, 
    mind, maxd,
    vmin, vmax, signed = TRUE)
  
  pdf("plot.pdf", width = 7, height = 3)
  par(mfrow = c(1,3), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0)
  plot(output$phip_list[[1]])
  plot(output$phip_list[[2]])
  plot(output$delta_list[[1]])
  plot(output$delta_list[[2]])
  plot(output$delta_list[[3]])
  dev.off()
  

  sum(output$phi0 != 0)
  sum(output$phip_list[[1]] != 0)
  sum(output$phip_list[[2]] != 0)
  sum(output$delta_list[[1]] != 0)
  sum(output$delta_list[[2]] != 0)
  sum(output$delta_list[[3]] != 0)
  
  sum(output$phi0 != 0 & output$phip_list[[1]] != 0)
  sum(output$phi0 != 0 & output$phip_list[[2]] != 0)

  sum(output$phi0 != 0 & output$delta_list[[1]] != 0)
  sum(output$phi0 != 0 & output$delta_list[[2]] != 0)
  sum(output$phi0 != 0 & output$delta_list[[3]] != 0)

  sum(output$phip_list[[1]] != 0 & output$phip_list[[2]] != 0)

  sum(output$phip_list[[1]] != 0 & output$delta_list[[1]] != 0)
  sum(output$phip_list[[1]] != 0 & output$delta_list[[2]] != 0)
  sum(output$phip_list[[1]] != 0 & output$delta_list[[3]] != 0)
  sum(output$phip_list[[2]] != 0 & output$delta_list[[1]] != 0)
  sum(output$phip_list[[2]] != 0 & output$delta_list[[2]] != 0)
  sum(output$phip_list[[2]] != 0 & output$delta_list[[3]] != 0)

}


######################################################################
######################################################################
## simulate_exogenous_vars:
##  Function to generate the matrix of time-invariant exogenous 
##  covariates X. 
##
##  INPUTS: 
##    p       : number of exogenous time-invariant covariates.
##    n       : number of subjects.
##    type    : if "unif", it generates independent uniform values in
##                $[-u_max.-u_min]\cup[u_min,u_max]$, or $[u_min,u_max]$,
##                depending on the parameter "signs." If "gaussian", it
##                takes a normal (or absolute value of normal) with sd
##                u_sd. 
##    u_min   : used if type = "unif". 
##    u_max   : used if type = "unif".  
##    g_sd    : used if type = "gaussian". 
##    signed  : whether we want only positive effects or signed effects.
##    nz_x_prob : probability of non-zero entries. If nz_x_prob = 1, we obtain
##                a pure uniform or gaussian distribution. If nz_x_prob, it is
##                a mixture of the distribution and an atom distribution at 
##                zero.
##
##  OUTPUT:
##    X       : n x p matrix of numeric exogenous time invariant covariates.
##
simulate_exogenous_vars <- function(p, n, type,
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
  mat_signs <- sample(signs, p * n, prob = nz_x_probs, TRUE)
  

  ## Generate entry magnitudes:
  mat_entries <- NULL
  if (type == "unif") {
    mat_entries <- runif(p * n, u_min, u_max)
  } else mat_entries <- abs(rnorm(p * n, mean = 0, sd = g_sd))

  X <- matrix(mat_signs * mat_entries, ncol = p, nrow = n)

  return(X)
}




######################################################################
######################################################################
## simulate_modvar1:
##  Simulates data according to the R-VAR model:
##    Yk_t = [Phi0 + Y_k1 Phi1 + ... + Y_kp Phip] Yk_{t-1} + Ek_{t-1}.
##
##  INPUTS:
##    modvar_pars1  : parameters for the MOD-VAR model generated from the
##                    function "generate_modvar_pars."
##    X             : Matrix of covariates 
##    n             : Number of subjects.
##    Ti            : Vector of time points for each subject.
##
##  OUTPUTS:
##    X             : n x p matrix of time-invariant exogenous covariates.
##    Y_list        : list of length n containing Ti[k] x d matrices.
##    B_list        : list of length n with d x d matrices with 
##                    individual time series parameters.
##    B_decmp       : list with the decomposed parameters
##                      of a time series.
##    err_covmat    : covariate matrix for the time series error.
##    d             : time series dimension.
##    p             : covariate dimension.
##
##
simulate_modvar1 <- function(modvar_pars1, X, n, Ti) {
  
  ## Simulate y values.
  d <- ncol(modvar_pars1$phi0)
  p <- ncol(X)
  
  ## Calculate subject specific parameters:
  B_list <- list()
  for (k_ind in 1:n) {
    xk_ind <- X[k_ind,]
    phip_list_mod <- modvar_pars1$phip_list
    
    for (l in 1:p) phip_list_mod[[l]] <- xk_ind[l] * phip_list_mod[[l]] 
    
    B_list[[k_ind]] <- modvar_pars1$phi0 + Reduce("+", phip_list_mod)
    B_list[[k_ind]] <- B_list[[k_ind]] + modvar_pars1$delta_list[[k_ind]]
  }
  
  ## Generate subject-specific time series data:
  Y_list <- list()
  for (k_ind in 1:n) {
    
    Y <- matrix(0, ncol = Ti[k_ind] + 400, nrow = d)
    
    E <- t(rmvnorm(Ti[k_ind] + 400, sigma = modvar_pars1$covD))
    
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
    B_dcmp      = modvar_pars1[1:3],
    err_covmat  = modvar_pars1$covD,
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
  d           <- 5
  p           <- 2
  n           <- 3
  prob_phi0   <- 0.20
  prob_phip   <- 0.05
  prob_delta  <- 0.05
  min0        <- 0.1
  max0        <- 0.2
  minp        <- 0.3
  maxp        <- 0.4
  mind        <- 0.5
  maxd        <- 0.6
  vmin        <- 0.3
  vmax        <- 0.5
  
  ## Exogenous parameters:
  type      <- "unif"
  u_min     <- 0.5
  u_max     <- 1
  signed    <- TRUE
  nz_x_prob   <- 0.7

  output <- generate_modvar_pars(
    d, p, n,
    prob_phi0, prob_phip, prob_delta,
    min0, max0, 
    minp, maxp, 
    mind, maxd,
    vmin, vmax)
  
  par(mfrow = c(3,1), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0)
  plot(output$phip_list[[1]])
  plot(output$phip_list[[2]])
  plot(output$delta_list[[1]])
  plot(output$delta_list[[2]])
  plot(output$delta_list[[3]])
  
  #########################
  ## Generate Data:
  #########################
  X         <- simulate_exogenous_vars(p = p, n = n, type = "unif",
                                       u_min = u_min, u_max = u_max,
                                       signed = signed, nz_x_prob = nz_x_prob)
  sims_data <- simulate_modvar1(output, X = X, n = n, Ti = rep(100, n))
  
  
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
       xlim = xrange, ylim = xrange, col = rep(c("red","black"), c(5, n-5)))
  text(x = sims_data$Y[,1], y = 0.3 + sims_data$Y[,2],  # Fine-tune the position
       label = c(1:5, rep("", n-5)), col = rep(c("red","black"), c(5, n-5))) 
  
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





######################################################################
######################################################################
## rvar_to_blist:
##  Given the p + n + 1 matrices d x d that determine the MOD-VAR model, this 
##  function constructs the n matrices d x d of VAR parameters.
##
##  INPUTS:
##    modvar_pars1  : parameters for the R-VAR model generated from the
##                    function "generate_modvar_pars."
##    X             : Matrix of moderators
##
##  OUTPUTS:
##    B_list        : list of length n with d x d matrices with 
##                    individual time series parameters.
##
modvar_to_blist <- function(modvar_pars1, X) {

  B_list <- list()
  for (k_ind in 1:nrow(X)) {
    xk_ind <- X[k_ind,]
    phip_list_mod <- modvar_pars1$phip_list
    
    for (l in 1:ncol(X)) {
      phip_list_mod[[l]] <- xk_ind[l] * phip_list_mod[[l]] 
    }

    B_list[[k_ind]] <- modvar_pars1$phi0 + Reduce("+", phip_list_mod)
    B_list[[k_ind]] <- B_list[[k_ind]] + modvar_pars1$delta_list[[k_ind]]
  }
  return(B_list)
}
  