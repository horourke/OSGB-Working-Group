
#################################################
#################################################
## idiographic_models:
##  Function that transforms population and regression
##  coefficients of an R-VAR model into a list of 
##  N matrices of individualized VAR coefficients.
##
##  INPUT:
##    Y                 : Nxp matrix of individual covariates.     
##    pop_coeffs        : dxd matrix of population level coefficients.
##    rvar_coeffs_list  : list of length p, with dxd matrices of 
##                          covariate-level effects in the R-VAR.
##
##  OUTPUT:
##    List of length N, where each element is the VAR model for
##      each individual.
##
idiographic_models <- function(Y, pop_coeffs, rvar_coeffs_list) {
    
  N <- nrow(Y)
  p <- ncol(Y)
  OUTPUT <- lapply(
    1:N,
    function(k, Y, pop_mat, ind_list) {
      coeff_mat <- pop_mat
           
      for(i in 1:p) {
        coeff_mat <- coeff_mat + Y[k,i] * ind_list[[i]] 
      }
      return(coeff_mat)
    },
    Y = Y, pop_mat = pop_coeffs, ind_list = rvar_coeffs_list)
  
  return(OUTPUT)
  
}

#################################################
#################################################

agg_absbias <- function (phi_N_est, phi_N_true) {
  
  
  N <- length(phi_N_est)
  
  ind_absbias_list <- mapply(
    function(x, y) {
      sum(abs(x-y))
    },
    x = phi_N_est, y = phi_N_true)
  
  
  return(sum(unlist(ind_absbias_list)) / N)
  
}


agg_absbias <- function (phi_N_est, phi_N_true) {
  
  
  N <- length(phi_N_est)
  
  ind_absbias_list <- mapply(
    function(x, y) {
      sum(abs(x-y))
    },
    x = phi_N_est, y = phi_N_true)
  
  
  return(sum(unlist(ind_absbias_list)) / N)
  
}






#################################################
#################################################
## Examples:

example <- FALSE

if (example) {
  
  #########################
  ## Generate parameters:
  #########################
  
  source("003_Generating_RvarData.R")
  source("041_rvar_supps.R")
  source("042_cv_rvar.R")
  source("043_bic_rvar.R")
  
  set.seed(20)
  set.seed(30)
  #set.seed(40)
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
  
  phi_bind_true <- cbind(output$phi0, 
                         do.call(cbind, output$phip_list))
  
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
  X_list <- sims_data$X_list
  Y      <- sims_data$Y
  verbose <- TRUE
  
  cv_model <- biccv.solve_rvar_glmnet(d = d, p = p, 
                                   X_list = X_list, Y = Y, 
                                   lambda.seq = lambda.seq, nfolds = 10,
                                   penalty.factor = penalty.factor, verbose = verbose)
  
  ## PLOT OF RESULTS:
  cv_model$rvar_opt_coeffs
  layout(
    matrix(c(1,2,3,
             4,4,4,
             5,5,5), 
           byrow = T, ncol = 3))
  
  sims_data$B_dcmp$phi0 %>% 
    #abs() %>%
    plot(main = "Joint Effect", breaks = 30)
  sims_data$B_dcmp$phip_list[[1]]%>% 
    #abs() %>%
    plot(main = "Individual Effect Y1", breaks = 30)
  sims_data$B_dcmp$phip_list[[2]] %>% 
    #abs() %>%
    plot(main = "Individual Effect Y2", breaks = 30)

  
  ## PLOT OF RESULTS:
  ## PLOT OF RESULTS:
  ## PLOT OF RESULTS:
  layout(
    matrix(c(1, 2, 3), 
           byrow = T, ncol = 1))
  
  phi_bind_true  %>% 
    plot(breaks = 30, border = NA)
  
  cv_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    #abs() %>%
    plot(breaks = 30, border = NA)
  
  cv_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    {abs(.) > 5e-5} %>%
    plot(border = NA)

  
  
  #########################
  ## Performance measures:
  #########################
  
  cv_model$rvar_opt_coeffs %>% 
    select(-lambda1, -lambda2, -var) %>%
    as.matrix() %>%
    agg_absbias(phi_N_true = list(phi_bind_true))

  ## Splitting estimate matrix into parts:
  pop_mat_est  <- cv_model$rvar_opt_coeffs[, 3 + (1:d)]
  ind_list_est <- lapply(
    1:p,
    function(k, mat, p, d){
      return(mat[, d * (k-1) + (1:d)])
    }, 
    mat = cv_model$rvar_opt_coeffs[, -(1:(d + 3))],
    p = p, d = d)
  
  ## Idiographic-format
  igr_est <- idiographic_models(
    Y = Y, 
    pop_coeffs = pop_mat_est, 
    rvar_coeffs_list = ind_list_est)
  
  igr_true <- idiographic_models(
    Y = Y, 
    pop_coeffs = output$phi0, 
    rvar_coeffs_list = output$phip_list)
  
  #############################
  ## Aggregated Bias:
  agg_absbias(igr_est, igr_true)
  
  
}

rm(example)
