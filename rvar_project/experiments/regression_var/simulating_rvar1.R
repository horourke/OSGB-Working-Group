
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
## Simulate rVAR.

simulate_rvar1 <- function(rvar_pars1, y_cov, N, Ti) {
  
  ## Simulate y values.
  d <- ncol(rvar_pars1$phi0)
  p <- ncol(y_cov)
  Y <- rmvnorm(n = N, sigma = y_cov)
  
  B_list <- list()
  for (n_ind in 1:nrow(Y)) {
    yn_ind <- Y[n_ind,]
    phip_list_mod <- rvar_pars1$phip_list
    
    for (l in 1:p) phip_list_mod[[l]] <- yn_ind[l] * phip_list_mod[[l]] 
    
    B_list[[n_ind]] <- rvar_pars1$phi0 + Reduce("+", phip_list_mod)
    
  }
  
  X_list <- list()
  
  for (n_ind in 1:N) {
    
    X <- matrix(0, ncol = Ti + 400, nrow = d)
    
    E <- t(rmvnorm(Ti + 400, sigma = rvar_pars1$covD))
    
    for (i in 2:(Ti + 400)) {
      X[, i] <- B_list[[n_ind]] %*% X[, i-1] + E[, i]
    }
    
    X_list[[n_ind]] <- X[, -(1:400)]
    rm(X)
  }
  
  output <- list(
    Y      = Y,
    X_list = X_list,
    B_list = B_list,
    B_dcmp = rvar_pars1[1:2],
    y_cov  = y_cov,
    x_cov  = rvar_pars1$covD,
    d      = d,
    p      = p)
    
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
  sims_data <- simulate_rvar1(output, y_cov = 0.5 * diag(p), N = N, Ti = 100)
  
  
  lapply(sims_data$X_list, dim)
  
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
  

    
  
}


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




#########################
## EXAMPLES:
#########################

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
  lambda.seq      <- 10^(seq(1, -3, length.out = 10))
  penalty.factor  <- 10^(seq(1, -1, length.out = 10))
  B_data <- solve_rvar_glmnet(d = d, p = p, sims_data$X_list, sims_data$Y, sims_data$p, 
                              lambda.seq = lambda.seq,
                              penalty.factor = penalty.factor)
  dim(B_data)
  
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
