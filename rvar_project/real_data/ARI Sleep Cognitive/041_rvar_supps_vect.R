######################################################################
######################################################################
## rearrange rVAR: 
##  Function to rearrange data to the regression format:
## 
##  INPUTS:
##    Y_list      : list of length N containing Tk x d matrices.
##    X           : N x p matrix of exogenous time-invariant covariates.
##    p           : number of time-invariant covariates.
##
##  OUTPUT:
##    response    : sum(T) x (d + 2) matrix with the multicolumn response.
##    covariates  : sum(T) x (d * (p + 1) + 2) matrix of multicolumn
##                    covariates.
##
rearrange_rvar_data <- function(Y_list, X, p) {
  
  d <- ncol(Y_list[[1]])
  p <- ncol(X)
  N <- length(Y_list)

  ## To generate full covariance matrix... 
  ## we first generate for each K = 1,2,...,N:
  covariate_mat <- lapply(
    1:N, 
    function(k, Y_list, Y) {
      p   <- ncol(Y)
      xk  <- X[k,]
      Yk  <- Y_list[[k]]
      
      ## Setting names of variables for Y
      y_var_names <- NULL
      if (is.null(colnames(Yk))) {
        y_var_names <- paste0("var", 1:ncol(Yk))
      } else {
        y_var_names <- colnames(Yk)
      }
      ## Setting names of variables for X
      x_var_names <- NULL 
      if (is.null(colnames(X))) {
        x_var_names <- paste0("x", 0:ncol(X))
      } else {
        x_var_names <- c("x0", colnames(X))
      }

      
      ## Construct covariate matrix for
      ## K-th subject:
      Yk_x <- lapply(
        1:(p+1),
        function(i, Y, x, names) {
          Yx <- x[i] * Y
          colnames(Yx) <- paste0(y_var_names, "_",names[i])
          return(Yx)
        }, 
        Y = Yk, x = c(1,xk), names = x_var_names) %>%
        {Reduce(cbind, .)}
      
      ## Add subject and time information:
      Yk_x <- Yk_x %>%
        as_tibble() %>%
        {mutate(.,
                subject = k,
                time = 1:nrow(.) + 1, ## Index for the time t that corresponds
                .before = 1)} %>%     ##  to its corresponding response time.
        filter(time != max(time))
      
      return(Yk_x)
    },
    Y_list, X) %>%

    ## Finally, merge all matrices:
    {Reduce(rbind, .)}
  
  
  response_mat <- lapply(
    1:N, 
    function(k, Y_list) {
      Yk <- Y_list[[k]]
      
      ## Setting names of variables
      y_var_names <- NULL
      if (is.null(colnames(Yk))) {
        y_var_names <- paste0("var", 1:ncol(Yk))
      } else {
        y_var_names <- colnames(Yk)
      }

      Y_clean <- 1 * Yk 
      colnames(Y_clean) <- y_var_names
      Y_clean <- Y_clean %>%
        as_tibble() %>%
        {mutate(., 
                subject = k,
                time = 1:nrow(.),
                .before = 1)} %>%
        filter(time != min(time))
      
      return(Y_clean)
    }, Y_list) %>%
    {Reduce(rbind, .)}
  
  output <- list(response = response_mat, covariates = covariate_mat)
  
  return(output)
}


######################################################################
######################################################################
## vectorize_rvar_data:
##  Function to rearrange data to the single linear regression format.
## 
##  INPUTS:
##    X_list      : list of length N containing d x T matrices.
##    Y           : N x p matrix of individual-level covariates.
##    p           : number of time-invariant covariates.
##
##  OUTPUT:
##    covariates_vectorized : 
##    response_vectorized   : 
##                    covariates.
##
vectorize_rvar_data <- function(Y_list, X, y_var_names) {
  
  N    <- length(Y_list)
  p    <- ncol(X)
  d    <- ncol(Y_list[[1]])
  data <- rearrange_rvar_data(Y_list, X, p)

  ## Vectorizing the response!
  response_vectorized <- data$response %>%
    pivot_longer(
        cols = -c(subject, time),
        names_to = "var_name",
        values_to = "response") %>% 

    mutate(var_name = factor(var_name, levels = y_var_names, ordered = TRUE)) %>% 
    arrange(var_name, subject, time) %>%
    mutate(var_name = as.character(var_name))

    ## Building meta-data for covariates.
    data$covariates <- data$covariates %>% 
      arrange(subject, time)
    meta_data <- tibble()
    for (i in 1:d) {
        meta_data <- bind_rows(meta_data, data$covariates %>% select(subject, time))
    }
    var_names <- data$covariates %>% select(-subject, -time) %>% colnames()
    
    new_names <- paste0(rep(var_names, d), "_", "task", rep(1:d, rep((p+1)*d, d)))

    ## Vectorizing with Kronecker product.
    covariates_vectorized <- data$covariates %>% 
      select(-subject, -time) %>%
      as.matrix() %>%
      { diag(d) %x% (.) } ## expand!

  colnames(covariates_vectorized) <- new_names
  covariates_vectorized <- covariates_vectorized %>% 
    as.data.frame() %>%
    tibble() %>%
    {bind_cols(meta_data, .)}

  OUTPUT <- list(
    covariates_vectorized = covariates_vectorized,
    response_vectorized = response_vectorized)
  return(OUTPUT)

}



#################################################
## process_coeffs:
##  Function that inputs a glmnet model, and rearranges
##  them to be according to the d x (p + 1) matrix format.
##
##  INPUTS:
##    d             : time series dimension.
##    p             : covariate dimension.
##    B_tibble      : M x (d * (p + 1))
##    glmnet_sparse : output GLMNET model.
##    pf_val        : penalty factor value corresponding to
##                     current run.
##
##  OUTPUTS:
##    B_update : (M + d) x (d * (p + 1) + 2) updated version of
##                the input B_tibble object.
##
process_coeffs <- function(d, p, B_tibble, glmnet_sparse, pf_val, 
                           y_var_names, x_var_names) {
  
  nlambda <- length(glmnet_sparse$lambda)
  
  B_update <- B_tibble
  c_names <- paste0(
    rep(y_var_names, p + 1), "_",
    rep(x_var_names, rep(d, p + 1)))
  r_names <- y_var_names

  for(lambda_ind in 1:nlambda) {
    
    ## Calculating penalty values:
    alpha <- (p + 1) / (p * pf_val + 1)
    beta <- (p * pf_val + pf_val) / (p * pf_val + 1)
    lambda1 <- glmnet_sparse$lambda[lambda_ind]
    lambda2 <- lambda1 * pf_val
    
    ## Merging B
    B_update <- glmnet_sparse$beta[, lambda_ind] %>%
      matrix(nrow = d, byrow = TRUE, 
             dimnames = list(NULL, c_names)) %>%
      as_tibble() %>%
      
      mutate(lambda1 = lambda1,
             lambda2 = lambda2,
             var     = r_names,
             .before = 1) %>%
      
      rbind(B_update)
  }
  
  return(B_update)
}




#########################
#########################
## Examples:
examples <- FALSE

if (examples) {
  
  #########################
  ## Import requirements:
  #########################
  
  source("003_Generating_RvarData.R")
  source("041_rvar_supps_vect.R")

  wd <- getwd()
  req_lib_dir <- paste0(wd,"/req_lib")
  .libPaths(req_lib_dir)
  
  library(plot.matrix)
  library(tidyverse)
  library(magrittr)
  library(mvtnorm)
  
  #########################
  ## Generate parameters:
  #########################
  
  set.seed(20)
  ## R-VAR parameters:
  d         <- 5
  p         <- 2
  prob_phi0 <- 0.1
  prob_phip <- 0.1
  min0      <- 0.3
  max0      <- 0.5
  minp      <- 0.2
  maxp      <- 0.3
  vmin      <- 0.3
  vmax      <- 0.5
  
  ## Exogenous parameters:
  type      <- "unif"
  u_min     <- 0.2
  u_max     <- 1
  signed    <- TRUE
  nz_x_prob   <- 0.9

  output <- generate_rvar_pars(d, p, 
                               prob_phi0, prob_phip, 
                               min0, max0, minp, maxp,
                               vmin, vmax)
  
  pdf("plot.pdf", width = 7, height = 3)
  par(mfrow = c(1,3), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0)
  plot(output$phip_list[[1]])
  plot(output$phip_list[[2]])
  dev.off()
  
  #########################
  ## Generate Data:
  #########################

  N <- 50
  X         <- simulate_exogenous_vars(p = p, N = N, type = "unif",
                                       u_min = u_min, u_max = u_max,
                                       signed = signed, nz_x_prob = nz_x_prob)
  sims_data <- simulate_rvar1(output, X = X, N = N, Ti = rep(100, N))
  
  
  ## Preparing Inputs of vectorize_rvar_data
  Y_list <- sims_data$Y_list
  X      <- sims_data$X
  y_var_names <- NULL
  if (is.null(colnames(Y_list[[1]]))) {
    y_var_names <- paste0("var", 1:ncol(Y_list[[1]]))
  } else {
    y_var_names <- colnames(Y_list[[1]])
  }
  ## You can use this environment for debugging.

}

rm(examples)