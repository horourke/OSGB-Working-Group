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
          colnames(xy) <- paste0("var",1:ncol(xy), "_y",name[i])
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

