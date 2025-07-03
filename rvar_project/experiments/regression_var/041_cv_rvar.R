######################################################################
######################################################################
## rearrange rVAR: 
##
## Function to rearrange data to the regression format:
rearrange_rvar_data <- function(X_list, Y, p) {
  
  d <- nrow(X_list[[1]])
  p <- ncol(Y)
  N <- length(X_list)
  
  ## We generate the covariates
  covariate_mat <- lapply(
    1:N, 
    function(k, X_list, Y) {
      p   <- ncol(Y)
      yk  <- Y[k,]
      xk  <- t(X_list[[k]])
      
      xk_y <- lapply(
        1:(p+1),
        function(i, x, y) {
          xy <- y[i] * x 
          name <- c(0, 1:length(y))
          colnames(xy) <- paste0("var",1:ncol(xy), "_",name[i])
          return(xy)
        }, 
        x = xk, y = c(1,yk)) %>%
        {Reduce(cbind, .)}
      
      xk_y <- xk_y %>%
        as_tibble() %>%
        {mutate(., 
                subject = k,
                time = 1:nrow(.),
                .before = 1)} %>%
        filter(time != max(time))
      
      #par(mfrow = c(1,1))
      #xk_y %>% select(-subject, -time) %>% cor() %>% abs() %>% plot()
      #plot(abs(cor(xk_y)))
      
      return(xk_y)
    },
    X_list, Y) %>%
    {Reduce(rbind, .)}
  
  
  response_mat <- lapply(
    1:N, 
    function(k, X_list) {
      xk <- t(X_list[[k]])
      
      colnames(xk) <- paste0("var", 1:ncol(xk))
      
      x_clean <- xk %>%
        as_tibble() %>%
        {mutate(., 
                subject = k,
                time = 1:nrow(.),
                .before = 1)} %>%
        filter(time != min(time))
      
      return(x_clean)
    }, X_list) %>%
    {Reduce(rbind, .)}
  
  output <- list(response = response_mat, covariates = covariate_mat)
  
  return(output)
}


######################################################################
######################################################################
## Solve rvar:
##
## Function to calculate the solutions with LM (no regularization)
solve_rvar_lm <- function(X_list, Y, p) {
  rdata <- rearrange_rvar_data(X_list, Y, p)
  d <- nrow(X_list[[1]])
  
  B <- NULL
  for(var in 1:d) {
    lm <- lm(unlist(rdata$response[,2 + var]) ~ as.matrix(rdata$covariates[, -c(1,2)]))  
    B <- rbind(B, lm$coefficients)
  }
  
  return(B)
  
}



## Function to calculate the solutions with GLMNET (regularization)
solve_rvar_glmnet <- function(d, p, X_list, Y, 
                              lambda.seq, penalty.factor , ...) {
  rdata <- rearrange_rvar_data(X_list, Y, p)
  rdata$response
  n_pf <- length(penalty.factor)
  
  B_tibble <- tibble()
  for (pf_val in penalty.factor) {
    alpha <- (p + 1) / (p * pf_val + 1)
    beta <- (p * pf_val + pf_val) / (p * pf_val + 1)
    
    glmnet_sparse <- glmnet(
      x = as.matrix(rdata$covariates[, -c(1,2)]),
      y = rdata$response[,-c(1:2)], family = "mgaussian",
      lambda  = lambda.seq,
      penalty.factor = rep(c(alpha, beta), c(d, d * p)), intercept = FALSE)#,...)
    
    B_tibble <- process_coeffs(d, p, B_tibble, glmnet_sparse, pf_val = pf_val)
    
  }
  B_tibble <- B_tibble %>% arrange(lambda1, lambda2, var)
  return(B_tibble)
  
}


process_coeffs <- function(d, p, B_tibble, glmnet_sparse, pf_val) {
  
  nlambda <- length(glmnet_sparse$lambda)
  
  B_update <- B_tibble
  
  for(lambda_ind in 1:nlambda) {
    
    ## Calculating penalty values:
    alpha <- (p + 1) / (p * pf_val + 1)
    beta <- (p * pf_val + pf_val) / (p * pf_val + 1)
    lambda1 <- glmnet_sparse$lambda[lambda_ind]
    lambda2 <- lambda1 * pf_val
    
    ## Merging B
    B_update <- sapply(glmnet_sparse$beta, 
                       function(x, ind) {return(x[, ind])}, ind = lambda_ind) %>%
      t() %>%
      as_tibble() %>%
      
      mutate(lambda1 = lambda1,
             lambda2 = lambda2,
             var     = 1:d,
             .before = 1) %>%
      
      rbind(B_update)
  }
  
  return(B_update)
}













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
##    X_list         : list of length N containing d x T matrices.
##    Y              : matrix of covariates for each individual.
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
cv.solve_rvar_glmnet <- function(d, p, X_list, Y, 
                                 lambda.seq, penalty.factor , 
                                 nfolds = 5, verbose = FALSE, ...) { ## nfolds < N-individuals.
  
  ###################
  ## Initializing:
  N     <- nrow(Y)
  rdata <- rearrange_rvar_data(X_list, Y, p)
  n_pf  <- length(penalty.factor)
  
  ###################
  ## CV-fold selection/setup:
  subject_folds <- sample(1:nfolds, size = N, replace = TRUE, )
  foldid <- subject_folds[rdata$covariates$subject]
  
  cv_fit_error_m <- matrix(NA, 
                           nrow = length(penalty.factor), 
                           ncol = length(lambda.seq))
  cv_fit_error_sd <- matrix(NA, 
                            nrow = length(penalty.factor), 
                            ncol = length(lambda.seq))
  
  ###################
  ## CV runs:
  for (pf_val_ind in seq_along(penalty.factor)) {
    pf_val <- penalty.factor[pf_val_ind]
    alpha <- (p + 1) / (p * pf_val + 1)
    beta <- (p * pf_val + pf_val) / (p * pf_val + 1)
    
    glmnet_sparse <- cv.glmnet(
      x = as.matrix(rdata$covariates[, -c(1,2)]),
      y = as.matrix(rdata$response[,-c(1:2)]), family = "mgaussian",
      lambda  = lambda.seq,
      foldid = foldid,
      
      penalty.factor = rep(c(alpha, beta), c(d, d * p)), intercept = FALSE)#,...)
    
    cv_fit_error_m[pf_val_ind, ]  <- glmnet_sparse$cvm
    cv_fit_error_sd[pf_val_ind, ] <- glmnet_sparse$cvsd
    
    if(verbose) print(pf_val_ind)

  }
  
  if(verbose) {
    plot(log(cv_fit_error_m  - min(cv_fit_error_m) + 1e-4) , 
         breaks = 100, 
         border = NA,
         main = "Cross-Validation log-mean Error",
         xlab = "Lambda",
         ylab = "Penalty Factor") 
  }
  
  ###################
  ## Optimal parameters:
  min_ind_arr <- which(cv_fit_error_m == min(cv_fit_error_m), arr.ind = TRUE)
  
  lambda_min_ind <- min_ind_arr[1,2]
  lambda_opt_val <- lambda.seq[lambda_min_ind]
  
  pf_min_ind <- min_ind_arr[1,1]
  pf_opt_val <- penalty.factor[pf_min_ind]
  
  
  ###################
  ## Fitting model with selected penalty factor:
  alpha_opt <- (p + 1) / (p * pf_opt_val + 1)
  beta_opt  <- (p * pf_opt_val + pf_opt_val) / (p * pf_opt_val + 1)
  
  glmnet_sparse <- glmnet(
    x = as.matrix(rdata$covariates[, -c(1,2)]),
    y = rdata$response[,-c(1:2)], family = "mgaussian",
    lambda  = lambda.seq,
    penalty.factor = rep(c(alpha_opt, beta_opt), c(d, d * p)), intercept = FALSE)
  
  glmnet_sparse_opt <- glmnet(
    x = as.matrix(rdata$covariates[, -c(1,2)]),
    y = rdata$response[,-c(1:2)], family = "mgaussian",
    lambda  = lambda_opt_val,
    penalty.factor = rep(c(alpha_opt, beta_opt), c(d, d * p)), intercept = FALSE)
  
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
    cv_error_m      = cv_fit_error_m,    ## cv_error_m      : cross-validation mean error.
    cv_error_sd     = cv_fit_error_sd,   ## cv_error_sd     : cross-validation SD of error.
    lambda_opt_val  = lambda_opt_val,    ## lambda_opt_val  : optimal lambda according to cross-validation error.
    pf_opt_val      = pf_opt_val,        ## pf_opt_val      : optimal Pen. Fact., according to cross-validation error.
    rvar_coeffs     = B_tibble,          ## rvar_coeffs     : matrix of RVAR coefficients corresponding to optimal Penalty Factors.
    rvar_opt_coeffs = B_tibble_opt,      ## rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
    rvar_glmnet_fit = glmnet_sparse)     ## rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
  
  return(output)
  
}
  


if (example) {
  
  #########################
  ## Generate parameters:
  #########################
  
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
  
  N <- 50
  sims_data <- simulate_rvar1(output, y_cov = 0.5 * diag(p), N = N, Ti = 100)
  
  
  #########################
  ## Visualizing parameters:
  #########################
  
  par(mfrow = c(3,3))
  col_lims <- seq(-0.6, 0.6, length.out = 10)
  plot(sims_data$B_dcmp$phi0, main = "Joint Effect", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[1]], main = "Individual Effect Y1", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[2]], main = "Individual Effect Y2", breaks = col_lims)
  
  yrange <- c(-max(abs(sims_data$Y)) - 0.1, 
              max(abs(sims_data$Y))+ 0.1)
  
  plot(sims_data$Y[,1], sims_data$Y[,2], 
       xlab = "Y1", ylab = "Y2", main = "Covariates Y",
       xlim = yrange, ylim = yrange, col = rep(c("red","black"), c(5, N-5)))
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
  
  sims_data$X_mat <- lapply(
    1:length(sims_data$X_list), 
    function(k, data) {
      x <- data[[k]]
      colnames(x) <- paste0("t", 1:ncol(x))
      x <- as_tibble(x) %>%
        mutate(subject = k, 
               var = 1:nrow(x), 
               .before = 1)
      return(x)},
    data = sims_data$X_list) %>%
    
    {Reduce(rbind, .)}
  
  sims_data$X_mat %>%
    
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
  ## SOLVING RVAR AND VISUALIZING
  #########################
  
  X_list <- sims_data$X_list
  Y      <- sims_data$Y
  
  B <- solve_rvar_lm(sims_data$X_list, sims_data$Y, sims_data$p)
  B0 <- B[, 1 + (1:5)]
  B1 <- B[, 6 + (1:5)]
  B2 <- B[, 11 + (1:5)]
  
  col_lims <- seq(-0.6, 0.6, length.out = 10)
  
  par(mfrow = c(2,3))
  plot(sims_data$B_dcmp$phi0, main = "Joint Effect", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[1]], main = "Individual Effect Y1", breaks = col_lims)
  plot(sims_data$B_dcmp$phip_list[[2]], main = "Individual Effect Y2", breaks = col_lims)
  
  plot(B0, main = "Estimated Joint Effect", breaks = col_lims)
  plot(B1, main = "Estimated Individual Effect Y1", breaks = col_lims)
  plot(B2, main = "Estimated Individual Effect Y2", breaks = col_lims)
  
  
  #########################
  ## SOLVING RVAR WITH GLMNET
  #########################
  lambda.seq      <- 10^(seq(1, -5, length.out = 20))
  penalty.factor  <- 10^(seq(3, -3, length.out = 20))
  verbose <- TRUE
  cv_model <- cv.solve_rvar_glmnet(d = d, p = p, sims_data$X_list, sims_data$Y, sims_data$p, 
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
    {abs(.) > 1e-2} %>%
    plot(border = NA)
  
  
  
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
