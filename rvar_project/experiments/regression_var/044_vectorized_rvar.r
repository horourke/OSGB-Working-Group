
source("041_rvar_supps.R")

source("rvar_project/experiments/regression_var/041_rvar_supps.R")

######################################################################
######################################################################
## rearrange rVAR: 
##  Function to rearrange data to the regression format:
## 
##  INPUTS:
##    X_list      : list of length N containing d x T matrices.
##    Y           : N x p matrix of individual-level covariates.
##    p           : number of time-invariant covariates.
##
##  OUTPUT:
##    response    : sum(T) x (d + 2) matrix with the multicolumn response.
##    covariates  : sum(T) x (d * (p + 1) + 2) matrix of multicolumn
##                    covariates.
##
vectorize_rvar_data <- function(X_list, Y) {
  
  N    <- length(X_list)
  p    <- ncol(Y)
  d    <- dim(X_list[[1]])[1]
  data <- rearrange_rvar_data(X_list, Y, p)

  dim(data$response)
  dim(data$covariates)
  

  ## Vectorizing the response!
  response_vectorized <- data$response %>% 
    pivot_longer(
        cols = -c(subject, time),
        names_to = "var_name",
        values_to = "response"
        ) %>% 
    arrange(var_name, subject, time)
    
    ## Building meta-data for covariates.
    meta_data <- tibble()
    for (i in 1:d) {
        meta_data <- bind_rows(meta_data, data$covariates %>% select(subject, time))
    }
    var_names <- data$covariates %>% select(-subject, -time) %>%
        colnames()
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
        covariates = covariates_vectorized,
        response = response_vectorized)
    return(OUTPUT)

}


biccv.solve_rvar_glmnet <- function(
    d, p, X_list, Y, 
    lambda.seq, penalty.factor , 
    type.measure = "bic",
    nfolds = 5, verbose = FALSE, ...) { ## nfolds < N-individuals.
  
  ###################
  ## Initializing:
  N     <- nrow(Y)
  rdata <- vectorize_rvar_data(X_list, Y)
  n_pf  <- length(penalty.factor)
  
  ###################
  ## CV-fold selection/setup:
  bic_fit_error <- matrix(NA, 
                          ncol = length(penalty.factor),
                          nrow = length(lambda.seq))
  
  for (pf_val_ind in seq_along(penalty.factor)) {
    
    pf_val <- penalty.factor[pf_val_ind]
    alpha <- (p + 1) / (p * pf_val + 1)
    beta <- (p * pf_val + pf_val) / (p * pf_val + 1)
    
    names <- colnames(rdata$covariates[-c(1,2)])
    pf.vectorized <- ifelse(str_detect(names, "y0"), alpha, beta)

    x_train <- as.matrix(rdata$covariates[-c(1,2)])
    y_train <- as.matrix(rdata$response[-c(1:3)])

    glmnet_sparse <- glmnet(
      x = x_train,
      y = y_train, family = "gaussian",
      lambda  = lambda.seq,
      penalty.factor = pf.vectorized, 
      intercept = FALSE)
    
    n_obs <- glmnet_sparse$nobs
    dfs   <- glmnet_sparse$df
    dev   <- glmnet_sparse$nulldev * glmnet_sparse$dev.ratio

    bic_fit_error[pf_val_ind, ] <- dfs * log(n_obs) - dev
           
    if(verbose) print(pf_val_ind)

  }
  
  if(verbose) {
    plot(log(bic_fit_error  - min(bic_fit_error) + 1e-4) , 
         breaks = 100, 
         border = NA,
         main = "Cross-Validation log-mean Error",
         xlab = "Lambda",
         ylab = "Penalty Factor") 
  }
  
  ###################
  ## Optimal parameters:
  min_ind_arr <- which(bic_fit_error == min(bic_fit_error), arr.ind = TRUE)
  
  lambda_min_ind <- min_ind_arr[1,2]
  lambda_opt_val <- lambda.seq[lambda_min_ind]
  
  pf_min_ind <- min_ind_arr[1,1]
  pf_opt_val <- penalty.factor[pf_min_ind]
  
  
  ###################
  ## Fitting model with selected penalty factor:
  alpha_opt <- (p + 1) / (p * pf_opt_val + 1)
  beta_opt  <- (p * pf_opt_val + pf_opt_val) / (p * pf_opt_val + 1)
  names <- colnames(rdata$covariates[-c(1,2)])
  pf.vectorized <- ifelse(str_detect(names, "y0"), alpha_opt, beta_opt)

  glmnet_sparse <- glmnet(
    x = x_train,
    y = y_train, family = "gaussian",
    lambda  = lambda.seq,
    penalty.factor = pf.vectorized, intercept = FALSE)
  
  glmnet_sparse_opt <- glmnet(
    x = x_train,
    y = y_train, family = "gaussian",
    lambda  = lambda_opt_val,
    penalty.factor = pf.vectorized, intercept = FALSE)
  
  glmnet_sparse_opt$beta ## NOT WORKING YET

  ###################
  ## Cleaning output:
  B_tibble <- tibble()
  B_tibble <- process_coeffs(d, p, B_tibble, glmnet_sparse,
                             pf_val = pf_val)
  B_tibble <- B_tibble %>% arrange(lambda1, lambda2, var)
  
  B_tibble_opt <- tibble()
  B_tibble_opt <- process_coeffs(d, p, B_tibble_opt, glmnet_sparse_opt, 
                                 pf_val = pf_val)
  B_tibble_opt <- B_tibble_opt %>% arrange(lambda1, lambda2, var)
  
  output <- list(
    lambda          = lambda.seq,        ## lambda          : Sequence of lambda used.
    penalty_factor  = penalty.factor,    ## penalty_factor  : Sequence of penalty factors used.
    bic_error       = bic_fit_error,     ## cv_error_m      : cross-validation mean error.
    lambda_opt_val  = lambda_opt_val,    ## lambda_opt_val  : optimal lambda according to cross-validation error.
    pf_opt_val      = pf_opt_val,        ## pf_opt_val      : optimal Pen. Fact., according to cross-validation error.
    rvar_coeffs     = B_tibble,          ## rvar_coeffs     : matrix of RVAR coefficients corresponding to optimal Penalty Factors.
    rvar_opt_coeffs = B_tibble_opt,      ## rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
    rvar_glmnet_fit = glmnet_sparse)     ## rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
  
  return(output)
  
}









#################################################
#################################################
#################################################
#################################################
### EXAMPLES:

example <- FALSE

if (example) {
  
  #########################
  ## Generate parameters:
  #########################
  
  source("rvar_project/experiments/regression_var/003_Generating_RvarData.R")
  source("rvar_project/experiments/regression_var/041_rvar_supps.R")
  library(magrittr)
  library(tidyverse)
  
  set.seed(20)
  d         <- 5
  p         <- 2
  prob_phi0 <- 0.35
  prob_phip <- 0.15
  min0      <- 0.3
  max0      <- 0.5
  minp      <- 0.3
  maxp      <- 0.5
  vmin      <- 0.3
  vmax      <- 0.5
  
  output <- generate_rvar_pars(d, p, 
                               prob_phi0, prob_phip, 
                               min0, max0, minp, maxp,
                               vmin, vmax)
  
  par(mfrow = c(3,1), mar = c(5.1, 4.1, 4.1, 4.1))
  col_lims <- seq(-0.6, 0.6, length.out = 10)
  plot(output$phi0, breaks = col_lims)
  plot(output$phip_list[[1]], breaks = col_lims)
  plot(output$phip_list[[2]], breaks = col_lims)
  
  #########################
  ## Generate Data:
  #########################
  
  N <- 8
  sims_data <- simulate_rvar1(output, y_cov = 0.5 * diag(p), N = N, Ti = 10)
  
  X_list  <- sims_data$X_list
  Y       <- sims_data$Y


  lambda.seq      <- 10^(seq(1, -5, length.out = 20))
  penalty.factor  <- 10^(seq(3, -3, length.out = 20))
  verbose <- TRUE



}