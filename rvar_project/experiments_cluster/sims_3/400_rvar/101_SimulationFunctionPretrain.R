

#####################################################
#####################################################
#####################################################


#################################################
#################################################
## FullSimulation:
##    Function that, given a set of simulation parameters
##     args, performs args$nsim simulation replicates of 
##     IPC-HD methods for hub estimation with parameters in
##     args.
##
##  INPUTS
##    args    : object output from function CreateParameters.
##                provides all info for simulation run.
##    
##  OUTPUT:
##    output  : data-frame of simulation outputs. 
##
FullSimulation <- function(args) {

  main_folder <- "400_rvar/"


  #####################################################
  ## Setup output-saving objects
  #####################################################
  
  ## What methods do you want to run?
  method_names  <- c(
    "rvar_bic", "rvar_cv")
  n_methods     <- length(method_names)

  ## Output of simulation:
  ## Dataframe with 
  output <- data.frame(
    ID      = 1:(args$nsim * n_methods),
    method  = character(args$nsim * n_methods),
    sim     = numeric(args$nsim * n_methods),
    time    = numeric(args$nsim * n_methods),
    l0_est  = numeric(args$nsim * n_methods),
    l0_true = numeric(args$nsim * n_methods),
    TPR     = numeric(args$nsim * n_methods),
    FPR     = numeric(args$nsim * n_methods),
    l1      = numeric(args$nsim * n_methods),
    Fr      = numeric(args$nsim * n_methods))
  attach(output)

  #################################################
  #################################################
  ## Cycle:
  sim_ind           <- 1
  loop_start_time   <- Sys.time()
  count             <- 1

  gather <- list()
  while (sim_ind < args$nsim + 1) {
    
    ############################
    ## Generate data:
    {
      print(paste("Generating Model and Data (",
                  round(100 * (sim_ind - 1) / args$nsim, 2),
                  "%", ")"))
      
      # Generate Phi0, Phi1, ... Phip parameters:
      rvar_model_pars <- generate_rvar_pars(
        args$d, args$p, args$prob_phi0, args$prob_phip,
        args$phi0_min, args$phi0_max,
        args$phip_min, args$phip_max,
        args$vmin, args$vmax)
      phi_list <- c(
        list(rvar_model_pars$phi0), 
        rvar_model_pars$phip_list) ## Concatenated phi0 and phip's.

      # Generate exogenous data.
      X  <-  simulate_exogenous_vars(
        p = args$p, N = args$N, type = args$type,
        u_min = args$u_min, u_max = args$u_max, g_sd = args$g_sd,
        signed = args$signed, nz_x_prob = args$nz_x_prob) 

      # Generate R-VAR Y data.
      data <- simulate_rvar1(
        rvar_pars1 = rvar_model_pars, 
        X = X, 
        N = args$N, 
        Ti = rep(args$T, args$N))

      Y_list <- data$Y_list

      gather$model_pars <- rvar_model_pars
      gather$phi_list   <- phi_list
      gather$X          <- X
      gather$data       <- data
  
    }
    ############################
    ############################
    ############################
    ######## Sparse Multi-VAR: adaptive
    {
      ## 
      print(paste0("Step ", sim_ind,": bic.solve_rvar"))
      start_time                  <- Sys.time()
      lambda.seq      <- 10^(seq(1, -5, length.out = 20))
      penalty.factor  <- 10^(seq(3, -3, length.out = 20))
      bic.model <- bic.solve_rvar_glmnet_vectorized(
        d = args$d, p = args$p,
        Y_list = Y_list, X = X,
        lambda.seq = lambda.seq,
        penalty.factor = penalty.factor,
        verbose = FALSE)
      print(bic.model$rvar_opt_coeffs)
      end_time                    <- Sys.time()

      ## Saving things! bla bla bla
      output[count, 2]      <- "rvar_bic"
      output[count, 3]      <- sim_ind
      output[count, 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[count, -(1:4)] <- eval_msr(data$B_list, bic.model$var_ind_coeffs)
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used BIC:", memory_in_bytes / (1024^3), "GB"))
      
      count <- count + 1

      gather$BICmodel <- bic.model

    }
    ############################
    ############################
    ############################
    ######## Sparse Multi-VAR: adaptive
    {
      print(paste0("Step ", sim_ind,": cv.solve_rvar"))
      start_time                  <- Sys.time()
      lambda.seq      <- 10^(seq(1, -5, length.out = 20))
      penalty.factor  <- 10^(seq(3, -3, length.out = 20))

      cv.model <- cv.solve_rvar_glmnet_vectorized(
        d = args$d, p = args$p,
        Y_list = Y_list, X = X,
        lambda.seq = lambda.seq,
        penalty.factor = penalty.factor,
        nfolds = 5,
        verbose = FALSE)
      print(cv.model$rvar_opt_coeffs)
      end_time                    <- Sys.time()

      ## Saving things! bla bla bla
      output[count, 2]      <- "rvar_cv"
      output[count, 3]      <- sim_ind
      output[count, 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[count, -(1:4)] <- eval_msr(data$B_list, cv.model$var_ind_coeffs)
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used CV:", memory_in_bytes / (1024^3), "GB"))

      count <- count + 1

      gather$CVmodel <- cv.model

    }
    ############################
    ######## Time analysis:
    {
      print(paste0("Step ", sim_ind,": Time Analysis."))
      
      ## If it will take more than X days to run,
      ## save results and leave.
      time_stamp <- Sys.time()
      current.rt.hour   <- 
        difftime(time_stamp, loop_start_time, units = "hours") %>%
        as.numeric()
      current.rt.days   <- 
        difftime(time_stamp, loop_start_time, units = "days") %>%
        as.numeric()
      mean.rt.days      <- current.rt.days / sim_ind
      expected.rt.days  <- current.rt.days + 1.5 * mean.rt.days
      ncompleted        <- sim_ind
      
      if (expected.rt.days >= args$running_days) { ## days.
        print(paste("---> Expected running time (+1):", 
                    round(expected.rt.days, digits = 4),
                    "days."))
        print("---> Stopping process...")
        sim_ind = args$nsim + 1
      }
      sim_ind <- sim_ind + 1
    }
  }

  #################################################
  ## Return output:

  print(output)

  gather$output <- output

  return(gather)

}
