
#################################################
#################################################
#################################################
#################################################
## cv.solve_rvar_glmnet: 
##    For data Xlist and Y, perform cross-validation
##    and model selection for the R-VAR model. 
##
##  INPUTS
##    d              : time series dimension.
##    p              : covariate dimension.
##    X_list         : list of length N containing Tk x d matrices.
##    Y              : matrix of exogenous time invariant  covariates
##                      for each individual.
##    lambda.seq     : numeric vector of lambda parameters.
##    penalty.factor : vector that determines the balance of 
##                      lambda1/lambda2 for the model.
##    nfolds         : number of folds in our CV procedure. Must be 
##                      less than the number of subjects in the study.
##    
##  OUTPUT
##    lambda          : Sequence of lambda used.
##    penalty_factor  : Sequence of penalty factors used.
##    cv_error_m      : cross-validation mean error.
##    cv_error_sd     : cross-validation SD of error.
##    lambda_opt_val  : optimal lambda according to cross-validation 
##                        error.
##    pf_opt_val      : optimal Pen. Fact., according to cross-validation 
##                        error.
##    rvar_coeffs     : matrix of RVAR coefficients corresponding to optimal 
##                        choice of Penalty Factors.
##    rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
##    rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
##
cv.solve_rvar_glmnet_vectorized <- function(
    d, p, Y_list, X, 
    lambda.seq, penalty.factor , 
    nfolds = 5, verbose = FALSE, ...) { ## nfolds < N-individuals.
  
  print("step1")
  ###################
  ## Initializing:
  N     <- nrow(X)
  rdata <- vectorize_rvar_data(Y_list, X)
  n_pf  <- length(penalty.factor)
  
  
  print("step2")
  ###################
  ## Preparing variable names:
  ## Setting names of variables for Y
  y_var_names <- NULL
  if (is.null(colnames(Y_list[[1]]))) {
    y_var_names <- paste0("var", 1:ncol(Y_list[[1]]))
  } else {
    y_var_names <- colnames(Y_list[[1]])
  }
  
  
  print("step3")
  ## Setting names of variables for X
  x_var_names <- NULL 
  if (is.null(colnames(X))) {
    x_var_names <- paste0("x", 0:ncol(X))
  } else {
    x_var_names <- c("x0", colnames(X))
  }

  
  print("step4")
  ###################
  ## Consistency check:
  response_ids   <- rdata$response_vectorized[, c(1,2)]
  covariates_ids <- rdata$covariates[, c(1,2)]
  if(sum(response_ids != covariates_ids) > 0) {
    stop("Error: Compatibility problem between covariates and response. Explore issue.")
  }

  
  print("step5")
  ###################
  ## CV-fold selection/setup:
  subject_folds <- sample(1:nfolds, size = N, replace = TRUE)
  foldid <- subject_folds[rdata$covariates$subject]
  
  cv_fit_error_m <- matrix(NA, 
                           nrow = length(penalty.factor), 
                           ncol = length(lambda.seq))
  cv_fit_error_sd <- matrix(NA, 
                            nrow = length(penalty.factor), 
                            ncol = length(lambda.seq))
  
  
  print("step6")
  ###################
  ## CV runs:
  for (pf_val_ind in seq_along(penalty.factor)) {

    print(paste0("step6.", 1))
    pf_val <- penalty.factor[pf_val_ind]
    alpha <- (p + 1) / (p * pf_val + 1)
    beta <- (p * pf_val + pf_val) / (p * pf_val + 1)
        
    pf_vec <- ifelse(
        str_detect(
            string = colnames(rdata$covariates_vectorized[, -c(1:2)]),
            pattern = "x0"), 
        alpha, beta)
        


    print(paste0("step6.", 2))
    #print((rdata$response_vectorized[,4])[[1]])
    #print(dim(as.matrix(rdata$covariates_vectorized[, -c(1:2)])))

    glmnet_sparse <- cv.glmnet(
      x = as.matrix(rdata$covariates_vectorized[, -c(1:2)]),
      y = as.matrix(rdata$response_vectorized[,4][[1]]),
      family = "gaussian",
      lambda  = lambda.seq,
      foldid = foldid,
      
      penalty.factor = pf_vec,
      intercept = FALSE, standardize = FALSE)
    
    print(paste0("step6.", 3))
    cv_fit_error_m[pf_val_ind, ]  <- glmnet_sparse$cvm
    cv_fit_error_sd[pf_val_ind, ] <- glmnet_sparse$cvsd
    
    if(verbose) print(pf_val_ind)

    print(paste0("step6.", 4))
  }
  
  print("step7")
  if(verbose) {
    plot(log(cv_fit_error_m  - min(cv_fit_error_m) + 1e-4) , 
         breaks = 100, 
         border = NA,
         main = "Cross-Validation log-mean Error",
         xlab = "Lambda",
         ylab = "Penalty Factor") 
  }
  
  print("step8")
  ###################
  ## Optimal parameters:
  min_ind_arr <- which(cv_fit_error_m == min(cv_fit_error_m), arr.ind = TRUE)
  
  lambda_min_ind <- min_ind_arr[1,2]
  lambda_opt_val <- lambda.seq[lambda_min_ind]
  
  pf_min_ind <- min_ind_arr[1,1]
  pf_opt_val <- penalty.factor[pf_min_ind]
  
  print("step9")
  ###################
  ## Fitting model with selected penalty factor:
  alpha_opt <- (p + 1) / (p * pf_opt_val + 1)
  beta_opt  <- (p * pf_opt_val + pf_opt_val) / (p * pf_opt_val + 1)
  pf_vec_opt <- ifelse(
        str_detect(
            string = colnames(rdata$covariates_vectorized[, -c(1:2)]),
            pattern = "x0"), 
        alpha_opt, beta_opt)

  print("step10")
  glmnet_sparse <- glmnet(
    x = as.matrix(rdata$covariates_vectorized[, -c(1:2)]),
    y = (rdata$response_vectorized[,4])[[1]], family = "gaussian",
    lambda  = lambda.seq,
    penalty.factor = pf_vec_opt,
    intercept = FALSE, standardize = FALSE)
  
  glmnet_sparse_opt <- glmnet(
    x = as.matrix(rdata$covariates_vectorized[, -c(1:2)]),
    y = (rdata$response_vectorized[,4])[[1]], family = "gaussian",
    lambda  = lambda_opt_val,
    penalty.factor = pf_vec_opt,
    intercept = FALSE, standardize = FALSE)
  
  print("step11")
  ###################
  ## Cleaning output:
  B_tibble <- tibble()
  B_tibble <- process_coeffs(d, p, B_tibble, glmnet_sparse,
                             pf_val = pf_val, 
                             y_var_names, x_var_names)
  B_tibble <- B_tibble %>% arrange(lambda1, lambda2, var)
  
  B_tibble_opt <- tibble()
  B_tibble_opt <- process_coeffs(d, p, B_tibble_opt, glmnet_sparse_opt, 
                                 pf_val = pf_val, 
                                 y_var_names, x_var_names)
  B_tibble_opt <- B_tibble_opt %>% arrange(lambda1, lambda2, var)
  
  print("step12")
  output <- list(
    lambda          = lambda.seq,        ## lambda          : Sequence of lambda used.
    penalty_factor  = penalty.factor,    ## penalty_factor  : Sequence of penalty factors used.
    cv_error_m      = cv_fit_error_m,    ## cv_error_m      : cross-validation mean error.
    cv_error_sd     = cv_fit_error_sd,   ## cv_error_sd     : cross-validation SD of error.
    lambda_opt_val  = lambda_opt_val,    ## lambda_opt_val  : optimal lambda according to cross-validation error.
    pf_opt_val      = pf_opt_val,        ## pf_opt_val      : optimal Pen. Fact., according to cross-validation error.
    rvar_coeffs     = B_tibble,          ## rvar_coeffs     : matrix of RVAR coefficients corresponding to optimal Penalty Factors.
    rvar_opt_coeffs = B_tibble_opt,      ## rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
    rvar_glmnet_fit = glmnet_sparse)     ## rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
  
  return(output)
  
}





#########################
## EXAMPLES:
#########################

example <- FALSE
if (example) {
  
  source("003_Generating_RvarData.R")
  source("041_rvar_supps_vect.R")

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
  u_min     <- 0.2
  u_max     <- 1
  signed    <- TRUE
  nz_x_prob   <- 0.9

  output <- generate_rvar_pars(d, p, 
                               prob_phi0, prob_phip, 
                               min0, max0, minp, maxp,
                               vmin, vmax)
  
  par(mfrow = c(3,1), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0, breaks = 10)
  plot(output$phip_list[[1]], breaks = 10)
  plot(output$phip_list[[2]], breaks = 10)
  
  #########################
  ## Generate Data:
  #########################
  N <- 100
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
  text(x = sims_data$X[,1], y = 0.3 + sims_data$X[,2],  # Fine-tune the position
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


  
  #########################
  ## SOLVING RVAR WITH GLMNET
  #########################

  lambda.seq      <- 10^(seq(1, -5, length.out = 20))
  penalty.factor  <- 10^(seq(3, -3, length.out = 20))
  verbose <- TRUE
  cv_model <- cv.solve_rvar_glmnet_vectorized(
    d = d, p = p, sims_data$Y_list, sims_data$X, sims_data$p, 
    lambda.seq = lambda.seq, nfolds = 10,
    penalty.factor = penalty.factor, verbose = verbose)
  
  ## PLOT OF RESULTS:
  cv_model$rvar_opt_coeffs
  layout(
    matrix(c(1,2,3,
             4,4,4,
             5,5,5), 
           byrow = T, ncol = 3))
  
  plot(sims_data$B_dcmp$phi0, main = "Joint Effect", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[1]], main = "Individual Effect Y1", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[2]], main = "Individual Effect Y2", breaks = col_lims)
  
  cv_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    plot(breaks = 100, border = NA)
  cv_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {abs(.) > 1e-10} %>%
    plot(border = NA)
    
    ## Visualization of non-zeros in log-scale
    par(mfrow = c(1,1))
    cv_model$rvar_opt_coeffs %>% 
        select(-lambda1, -lambda2, -var) %>%
        {log(abs(.) + 1e-10, base = 10)} %>%
        as.matrix() %>%
        hist(breaks = 100)
  
  
  B_data %>% select(lambda1, lambda2) %>%
    as.matrix() %>%
    {.[1:100,]} %>%
    log() %>%
    plot(border = NA, breaks = 30)
  
  
  par(mfrow = c(1,1))  
  B_data %>% 
    filter(lambda1 == lambda.seq[6]) %>%
    select(-lambda1, -lambda2, -var) %>%
    
    as.matrix() %>%
    
    plot(border = NA, breaks = 11)
  
}

rm(example)
