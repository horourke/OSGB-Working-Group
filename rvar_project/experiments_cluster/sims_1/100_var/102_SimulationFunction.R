

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

  main_folder <- "100_var/"


  #####################################################
  ## Setup output-saving objects
  #####################################################
  
  ## What methods do you want to run?
  method_names  <- c(
    "var_standard", "var_adaptive")
  n_methods     <- length(method_names)

  ## Output of simulation:
  ## Dataframe with 
  output <- data.frame(
    ID      = 1:(args$nsim * n_methods),
    method  = character(args$nsim * n_methods),
    measure = rep("inf_meas", args$nsim * n_methods),
    sim     = numeric(args$nsim * n_methods),
    time    = numeric(args$nsim * n_methods))
  
  for (var in 1:args$p) {
    output[[var + 5]] <- numeric(args$nsim * n_methods)
    colnames(output)[var + 5] <- paste0("var", var)
  }
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
      
      # Generate PM.
      model        <- generate_rvar_models_1(
        d = args$d, p = args$p, 
        prob_tot = args$prob_tot, prob_c = args$prob_c,
        entry_min = args$entry_min, entry_max = args$entry_max, 
        signed = args$signed,
        N = args$N, max_iter = 1000)

      Y <-        model$Y
      Phi_list <- model$Phi_list
      
      # Generate data.
      data      <- simulate_rvar_data_1(
        Phi_list = Phi_list, X_covar =  args$sigma2 * diag(args$d), 
        Y = Y, N = args$N, T = args$T)

      Xlist <- lapply(data$X_list, t)
  
    }
    ############################
    ############################
    ############################
    ######## Sparse VAR: STANDARD
    {
      ## 
      start_time                  <- Sys.time()
      results <- list()
      for (n_ind in 1:args$N) {
        model <- BigVAR::constructModel(
          Y = Xlist[[n_ind]], 
          p = 1,
          gran = c(50,10),
          struct = "Basic",
          cv = "Rolling",
          verbose = TRUE,
          ownlambdas = FALSE,
          model.controls=list(intercept=FALSE),
          linear = FALSE)
        
        results[[n_ind]] <- BigVAR::cv.BigVAR(model)
      }
      end_time                    <- Sys.time()

      ## Saving things! bla bla bla
      #output[count, 2]      <- "COR_Scr_IPCHD"
      #output[count, 3]      <- "inf_meas"
      #output[count, 4]      <- sim_ind
      #output[count, 5]      <- difftime(
      #    time1 = end_time, time2 = start_time, units = "s") %>%
      #    as.numeric()
      #output[count, -(1:5)] <- result_cor
      
      count <- count + 1
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

  return(output)

}
