

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

  main_folder <- "500_modvar1sd/"


  #####################################################
  ## Setup output-saving objects
  #####################################################
  
  ## Clean data:
  load("110_DataForSims.RData")

  ## What methods do you want to run?
  method_names  <- c(
    "modvar1sd_MO", 
    "modvar1sd_FI", 
    "modvar1sd_MI",
    "ada.modvar1sd_MO", 
    "ada.var1sd_FI", 
    "ada.modvar1sd_MI")

  n_methods     <- length(method_names)
  n_sub         <- length(df_list)

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

  print(dim(output))

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
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## "modvar_MO"
    {
      ## 
      print(paste0("Step ", sim_ind,".1: MOD-VAR MOD ONLY"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(2, -2, length.out = 30))
      ratios    <- 10^(seq(2, 0, length.out = 30))
      MO.model <- cv.modvar(
        Ylist = Y_list, 
        X = lcdf_norm, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        multi = FALSE,
        cv.type = "bysubject",
        nfolds = 6,
        sparse1sd = TRUE)
      end_time  <- Sys.time()

      ## Saving things! bla bla bla
      MO_time <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[(count - 1) * n_sub + (1:n_sub), 1]      <- "modvar1sd_MO"
      output[(count - 1) * n_sub + (1:n_sub), 2]      <- (1:n_sub)
      output[(count - 1) * n_sub + (1:n_sub), 4]      <- MO_time
      output[(count - 1) * n_sub + (1:n_sub), -(1:4)] <- t(eval_forecast(Y_forecast, MO.model$bysubject_coeffs, args$range, args$horizon))
      
      print(dim(
        output[(count - 1) * n_sub + (1:n_sub), -(1:4)]
      ))
      print(dim(
        t(eval_forecast(Y_forecast, MO.model$bysubject_coeffs, args$range, args$horizon))
      ))

      memory_in_bytes <- mem_used()
      print(paste0("Memory used MO:", memory_in_bytes / (1024^3), "GB"))
      
      count <- count + 1
    }
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## "modvar_FI"
    {
      print(paste0("Step ", sim_ind,".2: CV.MOD-VAR rolling window"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(2, -2, length.out = 30))
      ratios    <- 10^(seq(2, 0, length.out = 30))

      FI.model <- cv.modvar(
        Ylist = Y_list, 
        X = NULL, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        multi = TRUE,
        cv.type = "rolling",
        sparse1sd = TRUE)
      end_time                    <- Sys.time()


      ## Saving things! bla bla bla
      FI_time <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[(count - 1) * n_sub + (1:n_sub), 1]  <- "modvar1sd_FI"
      output[(count - 1) * n_sub + (1:n_sub), 2]  <- (1:n_sub)
      output[(count - 1) * n_sub + (1:n_sub), 4]  <- FI_time
      output[(count - 1) * n_sub + (1:n_sub), -(1:4)] <- t(eval_forecast(Y_forecast, FI.model$bysubject_coeffs, args$range, args$horizon))
      
      print(dim(
        output[(count - 1) * n_sub + (1:n_sub), -(1:4)]
      ))
      print(dim(
        t(eval_forecast(Y_forecast, FI.model$bysubject_coeffs, args$range, args$horizon))
      ))

      memory_in_bytes <- mem_used()
      print(paste0("Memory used FI:", memory_in_bytes / (1024^3), "GB"))

      count <- count + 1
    }
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## "modvar_MI"
    {
      print(paste0("Step ", sim_ind,".2: CV.MOD-VAR rolling window"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(2, -2, length.out = 30))
      ratios    <- 10^(seq(2, 0, length.out = 30))

      MI.model <- cv.modvar(
        Ylist = Y_list, 
        X = lcdf_norm, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        multi = TRUE,
        cv.type = "rolling",
        sparse1sd = TRUE)
      end_time                    <- Sys.time()


      ## Saving things! bla bla bla
      MI_time <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric()
      output[(count - 1) * n_sub + (1:n_sub), 1]  <- "modvar1sd_MI"
      output[(count - 1) * n_sub + (1:n_sub), 2]      <- (1:n_sub)
      output[(count - 1) * n_sub + (1:n_sub), 4]  <- MI_time
      output[(count - 1) * n_sub + (1:n_sub), -(1:4)] <- t(eval_forecast(Y_forecast, MI.model$bysubject_coeffs, args$range, args$horizon))
      
      print(dim(
        output[(count - 1) * n_sub + (1:n_sub), -(1:4)]
      ))
      print(dim(
        t(eval_forecast(Y_forecast, MI.model$bysubject_coeffs, args$range, args$horizon))
      ))

      memory_in_bytes <- mem_used()
      print(paste0("Memory used MI:", memory_in_bytes / (1024^3), "GB"))

      count <- count + 1
    }
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## "ada.modvar_MO"
    {
      ## 
      print(paste0("Step ", sim_ind,".1: adaBIC.MOD-VAR"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(2, -2, length.out = 30))
      ratios    <- 10^(seq(2, 0, length.out = 30))


      W.ada <- adaweights(
        MO.model, 
        p = ncol(lcdf_norm),
        multi = FALSE, 
        alpha = 1, 
        thr = 1e-5)

      adaMO.model <- cv.modvar(
        Ylist = Y_list,
        X = lcdf_norm, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = W.ada,
        multi = FALSE,
        cv.type = "bysubject", 
        nfolds = 6,
        sparse1sd = TRUE)
      end_time  <- Sys.time()

      ## Saving things! bla bla bla
      output[(count - 1) * n_sub + (1:n_sub), 1]      <- "ada.modvar1sd_MO"
      output[(count - 1) * n_sub + (1:n_sub), 2]      <- (1:n_sub)
      output[(count - 1) * n_sub + (1:n_sub), 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric() %>% {. + MO_time}
      output[(count - 1) * n_sub + (1:n_sub), -(1:4)] <- t(eval_forecast(Y_forecast, adaMO.model$bysubject_coeffs, args$range, args$horizon))
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used adaMO:", memory_in_bytes / (1024^3), "GB"))
      
      count <- count + 1
    }
    
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## "ada.modvar_FI"
    {
      ## 
      print(paste0("Step ", sim_ind,".1: adaBIC.MOD-VAR"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(2, -2, length.out = 30))
      ratios    <- 10^(seq(2, 0, length.out = 30))


      W.ada <- adaweights(
        FI.model, 
        p = 0,
        multi = TRUE, 
        alpha = 1, 
        thr = 1e-5)

      adaFI.model <- cv.modvar(
        Ylist = Y_list,
        X = NULL, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = W.ada,
        multi = TRUE,
        cv.type = "rolling",
        sparse1sd = TRUE)
      end_time  <- Sys.time()

      ## Saving things! bla bla bla
      output[(count - 1) * n_sub + (1:n_sub), 1]      <- "ada.modvar1sd_FI"
      output[(count - 1) * n_sub + (1:n_sub), 2]      <- (1:n_sub)
      output[(count - 1) * n_sub + (1:n_sub), 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric() %>% {. + FI_time}
      output[(count - 1) * n_sub + (1:n_sub), -(1:4)] <- t(eval_forecast(Y_forecast, adaFI.model$bysubject_coeffs, args$range, args$horizon))
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used adaFI:", memory_in_bytes / (1024^3), "GB"))
      
      count <- count + 1
    }
    ####################################################################################
    ####################################################################################
    ####################################################################################
    ######## "ada.modvar_MI"
    {
      ## 
      print(paste0("Step ", sim_ind,".1: adaBIC.MOD-VAR"))
      start_time                  <- Sys.time()
      lambdas1  <- 10^(seq(2, -2, length.out = 30))
      ratios    <- 10^(seq(2, 0, length.out = 30))


      W.ada <- adaweights(
        MI.model, 
        p = ncol(lcdf_norm),
        multi = TRUE, 
        alpha = 1, 
        thr = 1e-5)

      adaMI.model <- cv.modvar(
        Ylist = Y_list,
        X = lcdf_norm, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = W.ada,
        multi = TRUE,
        cv.type = "rolling",
        sparse1sd = TRUE)
      end_time  <- Sys.time()

      ## Saving things! bla bla bla
      output[(count - 1) * n_sub + (1:n_sub), 1]      <- "ada.modvar1sd_MI"
      output[(count - 1) * n_sub + (1:n_sub), 2]      <- (1:n_sub)
      output[(count - 1) * n_sub + (1:n_sub), 4]      <- difftime(
          time1 = end_time, time2 = start_time, units = "s") %>%
          as.numeric() %>% {. + MI_time}
      output[(count - 1) * n_sub + (1:n_sub), -(1:4)] <- t(eval_forecast(Y_forecast, adaMI.model$bysubject_coeffs, args$range, args$horizon))
      
      memory_in_bytes <- mem_used()
      print(paste0("Memory used adaMI:", memory_in_bytes / (1024^3), "GB"))
      
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
        sim_ind = n_sub + 1
      }
      sim_ind <- sim_ind + 1
    }
  

  #################################################
  ## Return output:

  print(output)

  return(output)

}
