
#################################################
#################################################
#################################################
#################################################
## solve_rvar_glmnet_adaptive: 
##    For data Xlist and Y, perform cross-validation
##    and model selection for the R-VAR model. 
##
##  INPUTS
##    d              : time series dimension.
##    p              : covariate dimension.
##    Y_list         : list of length N containing Tk x d matrices.
##    X              : matrix of exogenous time invariant  covariates
##                      for each individual.
##    gamma          : Adaptive parameter. 
##    rvar.fit       : the output of either
##                      "cv.solve_rvar_glmnet_vectorized", or
##                      "bic.solve_rvar_glmnet_vectorized".
##    lambda.seq     : numeric vector of lambda parameters.
##    nfolds         : number of folds in our CV procedure. Must be 
##                      less than the number of subjects in the study.
##    
##  OUTPUT
##    lambda          : Sequence of lambda used.
##    adaptive_weight : Sequence of penalty factors used.
##    cv_error_m      : cross-validation mean error.
##    cv_error_sd     : cross-validation SD of error.
##    lambda_opt_val  : optimal lambda according to cross-validation error.
##    rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
##    rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
##    var_ind_coeffs  : individual VAR models for each subject. 
##
solve_rvar_glmnet_adaptive <- function(
  d, p, Y_list, X, gamma = 1, rvar.fit = NULL,
  lambda.seq,
  nfolds = 5, verbose = FALSE) { ## nfolds < N-individuals.
  

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

  ###################
  ## Consistency check:
  response_ids   <- rdata$response_vectorized[, c(1,2)]
  covariates_ids <- rdata$covariates[, c(1,2)]
  if(sum(response_ids != covariates_ids) > 0) {
    stop("Error: Compatibility problem between covariates and response. Explore issue.")
  }

  ###################
  ## CV-fold selection/setup:
  res     <- N %% nfolds 
  cv_ind  <- rep(1:nfolds, floor(N / nfolds))
  if (res > 0) cv_ind <- c(cv_ind, 1:res)
  subject_folds <- cv_ind[sample(1:N, N, replace = FALSE)]
  foldid <- subject_folds[rdata$covariates$subject]

  ###################
  ## Train Adaptive RVAR
  pre_pf_vec <- ifelse(
    abs(rvar.fit$rvar_glmnet_fit$beta) == 0,
    -1,
    abs(rvar.fit$rvar_glmnet_fit$beta))

  pf_vec <- ifelse(
    pre_pf_vec == -1,
    Inf,
    1 / abs(pre_pf_vec) ^ gamma)

  if(sum(pre_pf_vec != -1) == 0) {

    lseq.length <- length(lambda.seq)
    lopt <- rvar.fit$lambda_opt_val
    error <- rvar.fit[[3]][rvar.fit$lambda == lopt]

    output <- list(
      lambda          = rvar.fit$lambda,          ## lambda          : Sequence of lambda used.
      adaptive_weight = pf_vec,                   ## adaptive_weight : Sequence of penalty factors used.
      cv_error_m      = rep(error, lseq.length),  ## cv_error_m      : cross-validation mean error.
      cv_error_sd     = rep(0, lseq.length),      ## cv_error_sd     : cross-validation SD of error.
      lambda_opt_val  = rvar.fit$lambda_opt_val,  ## lambda_opt_val  : optimal lambda according to cross-validation error.
      rvar_opt_coeffs = rvar.fit$rvar_opt_coeffs, ## rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
      rvar_glmnet_fit = rvar.fit$rvar_glmnet_fit, ## rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
      var_ind_coeffs  = rvar.fit$var_ind_coeffs)  ## var_ind_coeffs  : individual VAR models for each subject. 

    return(output)
  }


  ###################
  ## Train Adaptive RVAR

  glmnet_sparse <- cv.glmnet(
    x = as.matrix(rdata$covariates_vectorized[, -c(1:2)]),
    y = unlist(rdata$response_vectorized[,4][[1]]),
    family = "gaussian",
    lambda  = lambda.seq,
    foldid = foldid,
    penalty.factor = pf_vec,
    intercept = FALSE, standardize = FALSE)
    
  if (verbose) {
    print(pf_val_ind)
    memory_in_bytes <- mem_used()
    print(paste0("CV:", memory_in_bytes / (1024^3), "GB"))

    plot(glmnet_sparse$cvm, 
         main = "Cross-Validation mean Error",
         xlab = "CVM",
         ylab = "Lambda") 
  }
  
  ###################
  ## Fitting model with selected penalty factor:
  lambda_opt_val <- glmnet_sparse$lambda.min
    
  glmnet_sparse_opt <- glmnet(
    x = as.matrix(rdata$covariates_vectorized[, -c(1:2)]),
    y = (rdata$response_vectorized[,4])[[1]], family = "gaussian",
    lambda  = lambda_opt_val,
    penalty.factor = pf_vec,
    intercept = FALSE, standardize = FALSE)
  
  B_tibble_opt <- tibble()
  B_tibble_opt <- process_coeffs(d, p, B_tibble_opt, glmnet_sparse_opt, 
                                 pf_val = 1, 
                                 y_var_names, x_var_names) %>%
                    mutate(var = factor(var, levels = y_var_names, ordered = TRUE)) %>%     
                    arrange(lambda1, lambda2, var) %>%
                    mutate(var = as.character(var))
  
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
    B_total[[ind]] <- Bmat[, 1:d]
    for (ent_ind in 1:p) {
      B_total[[ind]] <- B_total[[ind]] + x[ent_ind] * Bmat[,ent_ind * d + (1:d)]
    }
  }

  output <- list(
    lambda          = lambda.seq,         ## lambda          : Sequence of lambda used.
    adaptive_weight = pf_vec,             ## adaptive_weight : Sequence of penalty factors used.
    cv_error_m      = glmnet_sparse$cvm,  ## cv_error_m      : cross-validation mean error.
    cv_error_sd     = glmnet_sparse$cvsd, ## cv_error_sd     : cross-validation SD of error.
    lambda_opt_val  = lambda_opt_val,     ## lambda_opt_val  : optimal lambda according to cross-validation error.
    rvar_opt_coeffs = B_tibble_opt,       ## rvar_opt_coeffs : matrix of optimal RVAR coefficients for PF and lambda.
    rvar_glmnet_fit = glmnet_sparse_opt,  ## rvar_glmnet_fit : unprocessed glmnet output for the rvar fit.
    var_ind_coeffs  = B_total)            ## var_ind_coeffs  : individual VAR models for each subject. 

  if (verbose) {
    memory_in_bytes <- mem_used()
    print(paste0("CV_F3:", memory_in_bytes / (1024^3), "GB"))  
  }

  return(output)
  
}





#########################
## EXAMPLES:
#########################

example <- FALSE
if (example) {
  
  #########################
  ## Import requirements:
  #########################
  
  source("002_requirements_lite.R")
  source("003_Generating_RvarData.R")
  source("041_rvar_supps_vect.R")
  source("045_cv.glmnet_rvar_vectorized.r")
  source("046_bic.glmnet_rvar_vectorized.r")

  
  #########################
  ## Generate parameters:
  #########################
  
  set.seed(20)
  ## R-VAR parameters:
  d         <- 5
  p         <- 2
  prob_phi0 <- 0.2
  prob_phip <- 0.1
  min0      <- 0.3
  max0      <- 0.5
  minp      <- 0.2
  maxp      <- 0.3
  vmin      <- 0.3
  vmax      <- 0.5
  signed    <- TRUE

  ## Exogenous parameters:
  type      <- "unif"
  u_min     <- 0.2
  u_max     <- 1
  signed    <- TRUE
  nz_x_prob   <- 0.9

  output <- generate_rvar_pars(d, p, 
                               prob_phi0, prob_phip, 
                               min0, max0, minp, maxp,
                               vmin, vmax, signed)
  
  pdf("plot.pdf", width = 7, height = 3)
  par(mfrow = c(1,3), mar = c(5.1, 4.1, 4.1, 4.1))
  plot(output$phi0, breaks = 10)
  plot(output$phip_list[[1]], breaks = 10)
  plot(output$phip_list[[2]], breaks = 10)
  dev.off()
  
  #########################
  ## Generate Data:
  #########################
  N <- 100
  Ti <- 100
  X         <- simulate_exogenous_vars(p = p, N = N, type = "unif",
                                       u_min = u_min, u_max = u_max,
                                       signed = signed, nz_x_prob = nz_x_prob)
  sims_data <- simulate_rvar1(output, X = X, N = N, Ti = rep(Ti, N))
  
  
  lapply(sims_data$Y_list, dim)
  
  lapply(sims_data$B_list, dim)
  lapply(sims_data$B_list, function(x) {sum(x != 0)})
  


  #########################
  ## Visualizing parameters:
  #########################
  
  pdf("plot.pdf", width = 7, height = 7)
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
  

  pdf("plot.pdf", width = 7, height = 7)
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

  lambda.seq      <- 10^(seq(1, -6, length.out = 20))
  penalty.factor  <- 10^(seq(4, -4, length.out = 20))
  verbose <- TRUE
  cv_model <- cv.solve_rvar_glmnet_vectorized(
    d = d, p = p, sims_data$Y_list, sims_data$X, sims_data$p, 
    lambda.seq = lambda.seq, nfolds = 10,
    penalty.factor = penalty.factor, verbose = FALSE)
  
  bic_model <- bic.solve_rvar_glmnet_vectorized(
    d = d, p = p, Y_list = sims_data$Y_list, X = sims_data$X, 
    lambda.seq = lambda.seq, penalty.factor = penalty.factor, 
    verbose = FALSE) 

  cv_model_ad <- solve_rvar_glmnet_adaptive(
    d = d, p = p,  Y_list = sims_data$Y_list, X = sims_data$X, 
    gamma = 3, rvar.fit = cv_model,
    lambda.seq = lambda.seq, nfolds = 10,
    verbose = FALSE)

  bic_model_ad <- solve_rvar_glmnet_adaptive(
    d = d, p = p,  Y_list = sims_data$Y_list, X = sims_data$X, 
    gamma = 4, rvar.fit = bic_model,
    lambda.seq = lambda.seq, nfolds = 10,
    verbose = FALSE)

  ######################
  ## CV: PLOT OF RESULTS
  cv_model$rvar_opt_coeffs
  
  pdf("plot1.pdf", width = 10, height = 12)
  layout(
    matrix(c(1,1,1,2,2,2,
             3,3,3,8,8,8,
             4,4,4,9,9,9,
             5,5,5,10,10,10,
             6,6,7,7,11,11), 
           byrow = T, ncol = 6))
  par(mar = c(5.1, 4.1, 4.1, 3.1))

  Bmat <- cbind(
    sims_data$B_dcmp$phi0,
    sims_data$B_dcmp$phip_list[[1]],
    sims_data$B_dcmp$phip_list[[2]])
  plot(Bmat, breaks = col_lims, border = NA)
  Bmat %>%
    {.[abs(.) > 1e-10]} %>%
    density() %>%
    plot()

  
  cv_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    plot(breaks = 100, border = NA)
  cv_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {abs(.) > 1e-10} %>%
    plot(border = NA)
  cv_model$rvar_opt_coeffs %>%
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {.[abs(.) > 1e-10 & abs(.) < 0.1]} %>%
    density(bw = 0.001) %>%
    plot()

  plot(log(cv_model$cv_error_m  - min(cv_model$cv_error_m) + 1e-4) , 
         breaks = 100, 
         border = NA,
         main = "Cross-Validation log-mean Error",
         xlab = "Lambda",
         ylab = "Penalty Factor") 
  plot(cv_model$cv_error_m, 
         breaks = 100, 
         border = NA,
         main = "Cross-Validation mean Error",
         xlab = "Lambda",
         ylab = "Penalty Factor") 

  cv_model_ad$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    plot(breaks = 100, border = NA)
  cv_model_ad$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {abs(.) > 1e-10} %>%
    plot(border = NA)
  cv_model_ad$rvar_opt_coeffs %>%
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {.[abs(.) > 1e-10 & abs(.) < 0.1]} %>%
    density(bw = 0.001) %>%
    plot()

  plot(cv_model_ad$cv_error_m)

  dev.off()





  ######################
  ## BIC: PLOT OF RESULTS
  pdf("plot2.pdf", width = 10, height = 12)
  layout(
    matrix(c(1,1,1,2,2,2,
             3,3,3,8,8,8,
             4,4,4,9,9,9,
             5,5,5,10,10,10,
             6,6,7,7,11,11), 
           byrow = T, ncol = 6))
  par(mar = c(5.1, 4.1, 4.1, 3.1))

  Bmat <- cbind(
    sims_data$B_dcmp$phi0,
    sims_data$B_dcmp$phip_list[[1]],
    sims_data$B_dcmp$phip_list[[2]])
  plot(Bmat, breaks = col_lims, border = NA)
  Bmat %>%
    {.[abs(.) > 1e-10]} %>%
    density(bw = 0.01) %>%
    plot()


  bic_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    plot(breaks = 100, border = NA)
  bic_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {abs(.) > 1e-10} %>%
    plot(border = NA)
  bic_model$rvar_opt_coeffs %>%
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {.[abs(.) > 1e-10]} %>%
    density(bw = 0.01) %>%
    plot()


  plot(log(bic_model$bic_error  - min(bic_model$bic_error) + 1e-4), 
         breaks = 100, 
         border = NA,
         main = "Cross-Validation log-mean Error",
         xlab = "Lambda",
         ylab = "Penalty Factor") 
  plot(bic_model$bic_error, 
         breaks = 100, 
         border = NA,
         main = "Cross-Validation log-mean Error",
         xlab = "Lambda",
         ylab = "Penalty Factor") 

  bic_model_ad$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    plot(breaks = 100, border = NA)
  bic_model_ad$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {abs(.) > 1e-10} %>%
    plot(border = NA)
  bic_model_ad$rvar_opt_coeffs %>%
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {.[abs(.) > 1e-10]} %>%
    density(bw = 0.01) %>%
    plot()


  plot(bic_model_ad$cv_error_m)

  dev.off()



  ######################
  ## Looking at individual RVAR results:
  cv_model_ad$rvar_opt_coeffs
  cv_model_ad$var_ind_coeffs[[1]]
  bic_model_ad$var_ind_coeffs[[1]]
  

}

rm(example)
