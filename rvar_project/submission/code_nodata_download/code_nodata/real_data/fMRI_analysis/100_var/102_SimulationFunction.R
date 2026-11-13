

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
  
  ## Clean data:
  load("110_DataForSims.RData")
  print(ls())

  ## What methods do you want to run?
  method_names  <- c(
    "var_standard")
  n_methods     <- length(method_names)
  n_sub          <- length(df_list)

  ## Output of simulation:
  ## Dataframe with 
  output <- data.frame(
    method  = character(n_sub * n_methods),
    subject = numeric(n_sub * n_methods),
    T       = rep(args$T, n_sub * n_methods),
    time    = numeric(n_sub * n_methods),
    rel_msfe1  = numeric(n_sub * n_methods),
    rel_msfe2  = numeric(n_sub * n_methods),
    rel_msfe3  = numeric(n_sub * n_methods),
    rel_msfe4  = numeric(n_sub * n_methods),
    rel_msfe5  = numeric(n_sub * n_methods),
    rel_msfe6  = numeric(n_sub * n_methods),
    rel_msfe7  = numeric(n_sub * n_methods),
    rel_msfe8  = numeric(n_sub * n_methods),
    rel_msfe9  = numeric(n_sub * n_methods),
    rel_msfe10 = numeric(n_sub * n_methods))
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
      print(paste("Generating Model and Data"))
      
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
        lapply(function(x) {
          x %>% 
          tibble() %>%
          as.matrix() %>%
          {(.)[1:args$T,]} %>%
          return()
        })
      
      Y_forecast <- dfn_list %>% 
        lapply(function(x) {
          x %>% 
          tibble() %>%
          as.matrix() %>%
          {(.)[-(1:args$T), ]} %>%
          return()
        })

        args$n <- length(Y_list)
    }
    ############################
    ############################
    ############################
    ######## Sparse VAR: STANDARD
    {

      start_time                  <- Sys.time()
      results <- list()
      for (n_ind in 1:args$n) {
        
        model <- BigVAR::constructModel(
          Y = Y_list[[n_ind]], 
          p = 1,
          gran = c(50,10),
          struct = "Basic",
          cv = "Rolling",
          verbose = TRUE,
          ownlambdas = FALSE,
          model.controls=list(intercept=FALSE),
          linear = FALSE)
        
        results[[n_ind]] <- BigVAR::cv.BigVAR(model)
        print(coef(results[[n_ind]]))

      }
      B_est_list <- 
        lapply(results, coef) %>%
        lapply(as.matrix) 
      end_time                    <- Sys.time()

      ## Saving things! bla bla bla
      output[(count - 1) * n_sub + (1:n_sub), 1]      <- "var_standard"
      output[(count - 1) * n_sub + (1:n_sub), 2]      <- (1:n_sub)
      output[(count - 1) * n_sub + (1:n_sub), 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()

      print(length(Y_forecast))
      print(length(B_est_list))
      output[(count - 1) * n_sub + (1:n_sub), -(1:4)] <- t(eval_forecast(Y_forecast, B_est_list, args$range, args$horizon))
      
      print(dim(
        output[(count - 1) * n_sub + (1:n_sub), -(1:4)]
      ))
      print(dim(
        t(eval_forecast(Y_forecast, B_est_list, args$range, args$horizon))
      ))

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
        sim_ind = n_sub + 1
      }
      sim_ind <- sim_ind + 1
    }

  #################################################
  ## Return output:
  print(output)

  return(output)

}
