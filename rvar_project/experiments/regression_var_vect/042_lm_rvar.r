
source("041_rvar_supps.R")

#################################################
#################################################
#################################################
#################################################
## solve_rvar_lm:
##  Function to calculate the solutions with 
##  LM (no regularization)
##
##  INPUTS:
##    X_list  : list of length N containing d x T matrices.
##    Y       : N x p matrix of individual-level covariates.
##    p       : number of time-invariant covariates.
##
##  OUTPUTS:
##    B       : d x (d * (p + 1)) matrix of R-VAR coefficients.
##
solve_rvar_lm <- function(Y_list, X, p) {
  rdata <- rearrange_rvar_data(Y_list, X, p)
  d <- ncol(Y_list[[1]])
  
  B <- NULL
  for(var in 1:d) {
    lm <- lm(unlist(rdata$response[,2 + var]) ~ as.matrix(rdata$covariates[, -c(1,2)]))  
    B <- rbind(B, lm$coefficients)
  }
  
  return(B)
  
}



example <- FALSE 

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


  #########################
  ## SOLVING RVAR AND VISUALIZING
  #########################

  B <- solve_rvar_lm(sims_data$Y_list, sims_data$X, sims_data$p)
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

  ## Things look correct! Great!  
  
}

rm(example)