

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

  main_folder <- "300_mvar/"


  #####################################################
  ## Setup output-saving objects
  #####################################################
  
  ## What methods do you want to run?
  method_names  <- c(
    "mvar_standard", "mvar_adaptive")
  
  n_methods     <- length(method_names)
  args$nsim     <- 1 

  ## Output of simulation:
  ## Dataframe with 
  output <- data.frame(
    ID      = 1:(args$nsim * n_methods),
    method  = character(args$nsim * n_methods),
    sim     = numeric(args$nsim * n_methods),
    time    = numeric(args$nsim * n_methods),
    rel_msfe1  = numeric(args$nsim * n_methods),
    rel_msfe2  = numeric(args$nsim * n_methods),
    rel_msfe3  = numeric(args$nsim * n_methods),
    rel_msfe4  = numeric(args$nsim * n_methods),
    rel_msfe5  = numeric(args$nsim * n_methods),
    rel_msfe6  = numeric(args$nsim * n_methods),
    rel_msfe7  = numeric(args$nsim * n_methods),
    rel_msfe8  = numeric(args$nsim * n_methods),
    rel_msfe9  = numeric(args$nsim * n_methods),
    rel_msfe10 = numeric(args$nsim * n_methods))
  attach(output)

  #################################################
  #################################################
  ## Cycle:
  sim_ind           <- 1
  loop_start_time   <- Sys.time()
  count             <- 1

    
    ############################
    ## Generate data:
    {
      print(paste("Generating Model and Data (",
                  round(100 * (sim_ind - 1) / args$nsim, 2),
                  "%", ")"))
      


      ## Clean data:
      load("103_Data.RData")

      dfn_list <- df_list %>% 
        lapply(
          function(x) {
            x %>%
            dplyr::select(-t,-sub) %>%
            mutate_if(is.numeric, scale) %>%
            return()
          }
        ) 

      ## Separate into forecasting and estimation sets:
      Y_list <- dfn_list %>% 
        tibble() %>%
        dplyr::select(-sub,-t) %>% 
        as.matrix() %>%
        {(.)[1:args$T]}
      
      Y_forecast <- dfn_list %>% 
        tibble() %>%
        dplyr::select(-sub,-t) %>% 
        as.matrix() %>%
        {(.)[-(1:args$T)]}
  
    }
    ############################
    ############################
    ############################
    ######## Sparse Multi-VAR: adaptive
    {
      ## 
      start_time                  <- Sys.time()
      model <- multivar::constructModel(data = Y_list, lassotype = "adaptive")
      fit <- multivar::cv.multivar(model)
      fit$mats %>% str()
      end_time                    <- Sys.time()

      ## Saving things! bla bla bla
      output[count, 2]      <- "mvar_adaptive"
      output[count, 3]      <- sim_ind
      output[count, 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[count, -(1:4)] <- eval_msr((fit$mats)$total, Y_forecast, args$range, args$horizon)
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used ADAPTIVE:", memory_in_bytes / (1024^3), "GB"))

      count <- count + 1
    }
    ############################
    ############################
    ############################
    ######## Sparse Multi-VAR: standard
    {
      ## 
      start_time                  <- Sys.time()
      model <- multivar::constructModel(data = Y_list, lassotype = "standard")
      fit <- multivar::cv.multivar(model)
      fit$mats %>% str()
      end_time                    <- Sys.time()

      ## Saving things! bla bla bla
      output[count, 2]      <- "mvar_standard"
      output[count, 3]      <- sim_ind
      output[count, 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[count, -(1:4)] <- eval_msr((fit$mats)$total, Y_forecast, args$range, args$horizon)
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used STANDARD:", memory_in_bytes / (1024^3), "GB"))

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

  #################################################
  ## Return output:
  print(output)

  return(output)

}
