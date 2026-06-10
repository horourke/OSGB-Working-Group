

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
FullSimulation <- function(args, microrun) {

  main_folder <- "400_modvar/"


  #####################################################
  ## Setup output-saving objects
  #####################################################
  
  ## What methods do you want to run?
  method_names  <- c(
    "modvar_bic", "modvar_cv_roll",
    "modvar_ada.bic", "modvar_ada.cv")
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
    Fr      = numeric(args$nsim * n_methods),
    sens    = numeric(args$nsim * n_methods),
    spec    = numeric(args$nsim * n_methods),
    fcast1  = numeric(args$nsim * n_methods),
    fcast2  = numeric(args$nsim * n_methods),
    fcast3  = numeric(args$nsim * n_methods))
  attach(output)

  #################################################
  #################################################
  ## Cycle:
  sim_ind           <- 1
  loop_start_time   <- Sys.time()
  count             <- 1

  while (sim_ind < args$nsim + 1) {
    
    ############################
    ## Generate data:
    {
      print(paste("Generating Model and Data (",
                  round(100 * (sim_ind - 1) / args$nsim, 2),
                  "%", ")"))
      
      # Setting common seed for each id_task*microrun*sim_ind
      seed_val <- 10000 * args$id_task  + 100 * microrun + sim_ind
      set.seed(seed_val)
      print(paste0("Seed val: ", seed_val))

      # Generate Phi0, Phi1, ... Phip parameters:
      modvar_model_pars <- generate_modvar_pars(
        args$d, args$p, args$n,
        args$prob_phi0, args$prob_phip, args$prob_delta, 
        args$phi0_min, args$phi0_max,
        args$phip_min, args$phip_max,
        args$delta_min, args$delta_max,
        args$vmin, args$vmax, signed = args$signed)
      phi_list <- c(
        modvar_model_pars$phi0, 
        modvar_model_pars$phip_list) ## Concatenated phi0 and phip's.

      # Generate exogenous data.
      X  <-  simulate_exogenous_vars(
        p = args$p, n = args$n, type = args$type,
        u_min = args$u_min, u_max = args$u_max, g_sd = args$g_sd,
        signed = args$signed, nz_x_prob = args$nz_x_prob) 

      # Generate R-VAR Y data.
      data <- simulate_modvar1(
        modvar_pars1 = modvar_model_pars, 
        X = X, 
        n = args$n, 
        Ti = rep(args$T +  2 * (args$range + args$horizon), args$n))

      Y_list <- lapply(
        data$Y_list,
        function(Y, Ti) {return(Y[1:Ti, ])},
         Ti = args$T)
      
      Y_forecast <- lapply(
        data$Y_list,
        function(Y, Ti) {return(Y[-(1:Ti), ])},
        Ti = args$T)
  
    }
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## BIC.MOD-VAR
    {
      ## 
      print(paste0("Step ", sim_ind,".1: BIC.MOD-VAR"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(1, -3, length.out = 20))
      ratios    <- 10^(seq(3, 0, length.out = 20))
      bic.model <- bic.modvar(
        Ylist = Y_list,
        X = X, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        multi = TRUE,
        alpha = 0.90)
      end_time  <- Sys.time()

      ## Saving things! bla bla bla
      output[count, 2]      <- "modvar_bic"
      output[count, 3]      <- sim_ind
      output[count, 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[count, -(1:4)] <- eval_msr(data$B_list, bic.model$bysubject_coeffs, Y_forecast, args$range, args$horizon)
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used BIC:", memory_in_bytes / (1024^3), "GB"))
      
      count <- count + 1
    }
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## CV.MOD-VAR rolling window
    {
      print(paste0("Step ", sim_ind,".2: CV.MOD-VAR rolling window"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(1, -3, length.out = 20))
      ratios    <- 10^(seq(3, 0, length.out = 20))

      cv.model <- cv.modvar(
        Ylist = Y_list, 
        X = X, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        multi = TRUE,
        cv.type = "rolling")
      end_time                    <- Sys.time()

      ## Saving things! bla bla bla
      cv_time <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[count, 2]      <- "modvar_cv_roll"
      output[count, 3]      <- sim_ind
      output[count, 4]      <- cv_time
      output[count, -(1:4)] <- eval_msr(data$B_list, cv.model$bysubject_coeffs, Y_forecast, args$range, args$horizon)
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used CV:", memory_in_bytes / (1024^3), "GB"))

      count <- count + 1
    }
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## adaBIC.MOD-VAR
    {
      ## 
      print(paste0("Step ", sim_ind,".1: adaBIC.MOD-VAR"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(1, -5, length.out = 10))
      ratios    <- 10^(seq(3, -3, length.out = 10))

      W.ada <- adaweights(
        bic.model, 
        p = args$p,
        multi = TRUE, 
        alpha = 1, 
        thr = 1e-5)

      ada.model <- cv.modvar(
        Ylist = Y_list, 
        X = X, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = W.ada,
        multi = TRUE, 
        cv.type = "rolling")
      end_time  <- Sys.time()

      ## Saving things! bla bla bla
      output[count, 2]      <- "modvar_ada.bic"
      output[count, 3]      <- sim_ind
      output[count, 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric() %>% {. + cv_time}
      output[count, -(1:4)] <- eval_msr(data$B_list, ada.model$bysubject_coeffs, Y_forecast, args$range, args$horizon)
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used BIC:", memory_in_bytes / (1024^3), "GB"))
      
      count <- count + 1
    }
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## CV.MOD-VAR by subject
    {
      print(paste0("Step ", sim_ind,".3: CV.MOD-VAR rolling"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(1, -5, length.out = 10))
      ratios    <- 10^(seq(3, -3, length.out = 10))

      W.ada <- adaweights(
        cv.model, 
        p = args$p,
        multi = TRUE, 
        alpha = 1, 
        thr = 1e-5)

      ada.model  <- cv.modvar(
        Ylist = Y_list, 
        X = X, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = W.ada,
        multi = TRUE, 
        cv.type = "rolling")
      end_time  <- Sys.time()

      ## Saving things! bla bla bla
      output[count, 2]      <- "modvar_ada.cv"
      output[count, 3]      <- sim_ind
      output[count, 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric() %>% {. + cv_time}
      output[count, -(1:4)] <- eval_msr(data$B_list, ada.model$bysubject_coeffs, Y_forecast, args$range, args$horizon)
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used CV:", memory_in_bytes / (1024^3), "GB"))

      count <- count + 1
    }
    ############################
    ######## Time analysis:
    {
      print(paste0("Step ", sim_ind,"5: Time Analysis."))
      
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

  return(output)

}
