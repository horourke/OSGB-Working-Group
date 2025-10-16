
#################################################
#################################################
#################################################
#################################################
## bic.solve_rvar_glmnet_vectorized: 
##    For data Xlist and Y, perform cross-validation
##    and model selection for the R-VAR model. 
##
##  INPUTS
##    d              : time series dimension.
##    p              : covariate dimension.
##    Y_list         : list of length N containing Tk x d matrices.
##    X              : matrix of exogenous time invariant  covariates
##                      for each individual.
##    lambda.seq     : numeric vector of lambda parameters.
##    penalty.factor : vector that determines the balance of 
##                      lambda1/lambda2 for the model.
##    
##  OUTPUT
##    lambda          : Sequence of lambda used.
##    penalty_factor  : Sequence of penalty factors used.
##    bic_error       : BIC across the multiple tuning parameters.
##    lambda_opt_val  : optimal lambda according to BIC.
##    pf_opt_val      : optimal Pen. Fact., according to BIC.
##    rvar_coeffs     : matrix of RVAR coefficients corresponding to optimal 
##                        choice of Penalty Factors.
##    rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
##    rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
##    var_ind_coeffs  : individual VAR models for each subject. 
##
bic.solve_rvar_glmnet_vectorized <- function(
    d, p, Y_list, X, 
    lambda.seq, penalty.factor , 
    verbose = FALSE, ...) { ## nfolds < N-individuals.
  
  ###################
  ## Preparing variable names:
  ## Setting names of variables for Y
  y_var_names <- NULL
  if (is.null(colnames(Y_list[[1]]))) {
    y_var_names <- paste0("var", 1:ncol(Y_list[[1]]))
  } else {
    y_var_names <- colnames(Y_list[[1]])
  }
  
  ## Setting names of variables for X
  x_var_names <- NULL 
  if (is.null(colnames(X))) {
    x_var_names <- paste0("x", 0:ncol(X))
  } else {
    x_var_names <- c("x0", colnames(X))
  }

  ###################
  ## Initializing:
  N     <- nrow(X)
  rdata <- vectorize_rvar_data(Y_list, X, y_var_names)
  n_pf  <- length(penalty.factor)

  ###################
  ## Consistency check:
  response_ids   <- rdata$response_vectorized[, c(1,2)]
  covariates_ids <- rdata$covariates[, c(1,2)]
  if(sum(response_ids != covariates_ids) > 0) {
    stop("Error: Compatibility problem between covariates and response. Explore issue.")
  }

  ###################
  ## CV-fold selection/setup:
  bic_fit_error <- matrix(NA, 
                          ncol = length(penalty.factor),
                          nrow = length(lambda.seq))
  
  if (verbose) {
    memory_in_bytes <- mem_used()
    print(paste0("BIC_B:", memory_in_bytes / (1024^3), "GB"))
  }

  ###################
  ## CV runs:
  for (pf_val_ind in seq_along(penalty.factor)) {
    pf_val <- penalty.factor[pf_val_ind]
    alpha <- (p + 1) / (p * pf_val + 1)
    beta <- (p * pf_val + pf_val) / (p * pf_val + 1)
        
    pf_vec <- ifelse(
        str_detect(
            string = colnames(rdata$covariates_vectorized[, -c(1:2)]),
            pattern = "x0"), 
        alpha, beta)
        


    glmnet_sparse <- glmnet(
      x = as.matrix(rdata$covariates_vectorized[, -c(1:2)]),
      y = (rdata$response_vectorized[,4])[[1]], family = "gaussian",
      lambda  = lambda.seq,
      penalty.factor = pf_vec,
      intercept = FALSE, standardize = FALSE)
    
    n_obs <- glmnet_sparse$nobs
    dfs   <- glmnet_sparse$df
    dev   <- glmnet_sparse$nulldev * glmnet_sparse$dev.ratio

    bic_fit_error[pf_val_ind, ] <- dfs * log(n_obs) - dev
           
    if(verbose){
      print(pf_val_ind)

      memory_in_bytes <- mem_used()
      print(paste0("BIC_", pf_val_ind, ":", memory_in_bytes / (1024^3), "GB"))
    }

    rm(glmnet_sparse)
  }
  
  if(verbose) {
    plot(log(bic_fit_error  - min(bic_fit_error) + 1e-4) , 
         breaks = 100, 
         border = NA,
         main = "log-BIC",
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
  pf_vec_opt <- ifelse(
        str_detect(
            string = colnames(rdata$covariates_vectorized[, -c(1:2)]),
            pattern = "x0"), 
        alpha_opt, beta_opt)


  #glmnet_sparse <- glmnet(
  #  x = as.matrix(rdata$covariates_vectorized[, -c(1:2)]),
  #  y = (rdata$response_vectorized[,4])[[1]], family = "gaussian",
  #  lambda  = lambda.seq,
  #  penalty.factor = pf_vec_opt,
  #  intercept = FALSE, standardize = FALSE)
  
  glmnet_sparse_opt <- glmnet(
    x = as.matrix(rdata$covariates_vectorized[, -c(1:2)]),
    y = (rdata$response_vectorized[,4])[[1]], family = "gaussian",
    lambda  = lambda_opt_val,
    penalty.factor = pf_vec_opt,
    intercept = FALSE, standardize = FALSE)
  
  if (verbose) {
    memory_in_bytes <- mem_used()
    print(paste0("BIC_F1:", memory_in_bytes / (1024^3), "GB"))
  }
  
  ###################
  ## Cleaning output:
  #B_tibble <- tibble()
  #B_tibble <- process_coeffs(d, p, B_tibble, glmnet_sparse,
  #                           pf_val = pf_val, 
  #                           y_var_names, x_var_names)
  #B_tibble <- B_tibble %>% arrange(lambda1, lambda2, var)
  
  B_tibble_opt <- tibble()
  B_tibble_opt <- process_coeffs(d, p, B_tibble_opt, glmnet_sparse_opt, 
                                 pf_val = pf_opt_val, 
                                 y_var_names, x_var_names) %>%
                    mutate(var = factor(var, levels = y_var_names, ordered = TRUE)) %>%     
                    arrange(lambda1, lambda2, var) %>%
                    mutate(var = as.character(var))

  memory_in_bytes <- mem_used()
  print(paste0("BIC_F2:", memory_in_bytes / (1024^3), "GB"))

  ###################
  ## Creating individual VAR coefficient matrices
  Bmat <- B_tibble_opt %>%
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {(.) * (abs(.) > 1e-10)}

  B_total <- list()
  for (ind in 1:N) {
    x <- X[ind, ]
    prod_mat <- diag(d)
    for (ent_ind in 1:p) 
      prod_mat <- rbind(prod_mat, x[ent_ind] * diag(d))

    B_total[[ind]] <- Bmat %*% prod_mat
  }


  
  output <- list(
    lambda          = lambda.seq,        ## lambda          : Sequence of lambda used.
    penalty_factor  = penalty.factor,    ## penalty_factor  : Sequence of penalty factors used.
    bic_error       = bic_fit_error,     ## cv_error_m      : cross-validation mean error.
    lambda_opt_val  = lambda_opt_val,    ## lambda_opt_val  : optimal lambda according to cross-validation error.
    pf_opt_val      = pf_opt_val,        ## pf_opt_val      : optimal Pen. Fact., according to cross-validation error.
    #rvar_coeffs     = B_tibble,         ## rvar_coeffs     : matrix of RVAR coefficients corresponding to optimal Penalty Factors.
    rvar_opt_coeffs = B_tibble_opt,      ## rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
    rvar_glmnet_fit = glmnet_sparse_opt, ## rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
    var_ind_coeffs  = B_total)           ## var_ind_coeffs  : individual VAR models for each subject. 
  
  memory_in_bytes <- mem_used()
  print(paste0("BIC_F3:", memory_in_bytes / (1024^3), "GB"))


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
  wd <- getwd()
  .libPaths(paste0(wd, "/req_lib/"))
  library(plot.matrix)
  library(tidyverse)
  library(magrittr)
  library(mvtnorm)
  
  set.seed(20)
  ## R-VAR parameters:
  d         <- 15
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
  
  par(mfrow = c(3,1), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0, breaks = 10)
  plot(output$phip_list[[1]], breaks = 10)
  plot(output$phip_list[[2]], breaks = 10)
  dev.off()

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
  text(x = sims_data$X[,1], y = sims_data$X[,2],  # Fine-tune the position
       label = c(1:5, rep("", N-5)), col = rep(c("red","black"), c(5, N-5)),
       pos = 3) 
  
  plot(sims_data$B_list[[1]], main = "Sample 1", breaks = col_lims)
  plot(sims_data$B_list[[2]], main = "Sample 2", breaks = col_lims)
  plot(sims_data$B_list[[3]], main = "Sample 3", breaks = col_lims)
  plot(sims_data$B_list[[4]], main = "Sample 4", breaks = col_lims)
  plot(sims_data$B_list[[5]], main = "Sample 5", breaks = col_lims)
  dev.off()  
  
  
  
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
    dev.off()

  
  #########################
  ## SOLVING RVAR WITH GLMNET
  #########################

  Y_list          <- sims_data$Y_list
  X               <- sims_data$X
  lambda.seq      <- 10^(seq(1, -5, length.out = 20))
  penalty.factor  <- 10^(seq(3, -3, length.out = 20))
  verbose <- TRUE
  cv_model <- bic.solve_rvar_glmnet_vectorized(
    d = d, p = p, sims_data$Y_list, sims_data$X, sims_data$p,
    lambda.seq = lambda.seq, nfolds = 10,
    penalty.factor = penalty.factor, verbose = verbose)
  
  phiest_list <- cv_model$rvar_opt_coeffs %>%
    select(-lambda1, -lambda2, -var) %>% 
    {lapply(0:p, function(x_val, data) {
        #print(colnames(data))
        print(paste0("x", x_val))
        which_keep <- str_detect(colnames(data), paste0("x", x_val))
        mat_filtered <- data %>% select_if(which_keep) %>% 
        as.matrix() 
        print(colnames(mat_filtered))
        return(mat_filtered)
    }, data = .)}

  ## PLOT OF RESULTS:
  cv_model$rvar_opt_coeffs
  par(mfcol = c(3,3))
  
  plot(sims_data$B_dcmp$phi0, main = "Joint Effect", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[1]], main = "Individual Effect Y1", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[2]], main = "Individual Effect Y2", breaks = col_lims)
  
  lapply(phiest_list, plot, breaks = 10)
  lapply(phiest_list, function(x) {
    plot(x!= 0)
  })
  dev.off()

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

  
}

rm(example)
