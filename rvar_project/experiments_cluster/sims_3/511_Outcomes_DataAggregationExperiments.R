
############################################
############################################
## SETTING PARAMETERS:
wd <- getwd()
.libPaths(paste0(wd,"/req_lib"))
library(readr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(stringr)
library(lazyeval)


###################### Parameter table:
runtype       <- 2 # FOR EXPERIMENT RUNS
#runtype       <- 3 # FOR FULL RUNS
index_old     <- 1 # run index to use
sim_par_table <- expand.grid(
  running_days  = 2,
  entry_min     = 0.3,                            ## entry_min : minimum B0,B1,...,Bp entry magnitude
  entry_max     = 0.9,                            ## entry_max : maximum B0,B1,...,Bp entry magnitude
    
  g_sd          = 0.3,                            ## g_sd  : 
  u_min         = 0.5,                            ## u_min : minimum entry of exogenous X for type = "unif".
  u_max         = 1,                              ## u_max : maximum entry of exogenous X for type = "unif".
  type          = "unif",                         ## type  : distribution of exogenous variables.
  nz_x_prob     = c(0.75),                        ## nz_x_prob : proportion of entries in exogenous
                                                    ##              data matrix X with non-zero values.
  signed        = c(TRUE),                        ## signed    : are entries signed or all positive?


  prob_c        = c(2/3, 1/3),                    ## prob_c   : proportion of common entries.
  prob_tot      = 0.1,                            ## prob_tot : total proportion of non-zero entries.

  nsim          = ifelse(
                    runtype == 1, 
                    1, 
                    ifelse(runtype == 2, 2, 10)),## nsim     : no of simulation repetitions.
  sigma2        = c(0.05, 0.1),                  ## sigma2   : variance o VAR error term.
  N             = c(20, 50),                    ## N        : No. of individuals
  T             = c(50, 100),                    ## T        : timepoints per individual.
  p             = c(2, 5),                       ## p        : covariate dimension
  d             = c(10, 20))                     ## d        : Time series dimension
## Add specific parameters to table:
sim_par_table$prob_phi0 <- sim_par_table$prob_tot * sim_par_table$prob_c
sim_par_table$prob_phip <- sim_par_table$prob_tot * (1 - sim_par_table$prob_c)
  
sim_par_table$phi0_min <- sim_par_table$entry_min
sim_par_table$phi0_max <- sim_par_table$entry_max
sim_par_table$phip_min <- sim_par_table$entry_min
sim_par_table$phip_max <- sim_par_table$entry_max
  
sim_par_table$vmax <- sim_par_table$sigma2
sim_par_table$vmin <- 0.5 * sim_par_table$sigma2
attach(sim_par_table)


###################### Creating folders:
subfolder_new        <- paste0("500_AggregatedDataExperiments/")
subfolder_data_new   <- paste0(subfolder_new, "data_all/")
subfolder_plots_new  <- paste0(subfolder_new, "plots_all/")

if (!dir.exists(subfolder_new)) {
       dir.create(subfolder_new)
}
if (!dir.exists(subfolder_data_new)) {
       dir.create(subfolder_data_new)
}
if (!dir.exists(subfolder_plots_new)) {
       dir.create(subfolder_plots_new)
}

##################################################################
##################################################################
## LOAD + MERGE RESULTS:
##################################################################
##################################################################
# Which simulations are we using?
##  RUNTYPE = 2: experiment partial runs.
##  RUNTYPE = 3: full simulation runs. 
run_info <- list(
  list(
    main_dir       = "100_var/",
    run_index      = 1,
    runtype        = runtype,
    abrev_name     = "VAR"),

  list(
    main_dir       = "200_gimme/",
    run_index      = 1,
    runtype        = runtype,
    abrev_name     = "GIMME"),
    
  list(
    main_dir       = "300_mvar/",
    run_index      = 1,
    runtype        = runtype,
    abrev_name     = "M-VAR"),

  list(
    main_dir       = "400_rvar/",
    run_index      = 1,
    runtype        = runtype,
    abrev_name     = "RVAR")
  )

###########################
## LOOP OVER ALL 288 SIMULATION PARAMETER COMBINATIONS
for (id_task in 1:64) {
  print(paste("XXXXXXXXXXXXXXXX ID-TASK", id_task))

  output <- NULL

  for(id_microrun in 0:9) {

    print(paste("XXXXXXXX MICRO-RUN", id_microrun))

    ###########################
    ## Merge all data, from all methods into output.
    for (method_ind in c(1,3,4)) {
      main_dir     <- run_info[[method_ind]]$main_dir
      run_index    <- run_info[[method_ind]]$run_index
      runtype      <- run_info[[method_ind]]$runtype
      runtype_name <- c("pretrainings", "experiments", "outputs")[runtype]
      abrev_name   <- run_info[[method_ind]]$abrev_name

      load(paste0(
        main_dir, runtype_name, run_index,
        "/data/output", id_task, "_", id_microrun, ".RData"))

      ## Add id-task info:
      output_temp <- get(paste0("output", id_task, "_", id_microrun))
      args_temp   <- get(paste0("args", id_task))
      output_temp <- output_temp %>%
        mutate(
          d = args_temp$d,
          p = args_temp$p,
          T = args_temp$T,
          N = args_temp$N,
          sigma2 = args_temp$sigma2,
          prob_c = args_temp$prob_c,
          signed = args_temp$signed,
          .before = method)
      
      ## Clean names:
      colnames(output_temp) <- gsub(" ", "", colnames(output_temp))
      
      ## Visualize current data:
      print(paste(
        abrev_name, nrow(output_temp), ncol(output_temp)))

      ## Add data to output:
      output <- rbind(
        output,
        output_temp)
      rm(output_temp, args_temp)
      rm(list = paste0("output", id_task, "_", id_microrun))
    }
  }
  


  ###########################
  ## Saving data:
  print(paste("Saving data: Task-ID", id_task))
  print(paste0(subfolder_data_new, "output", id_task, ".RData"))
  
  args <- get(paste0("args", id_task))
  print(paste(
      "TOTAL",
      nrow(output),
      ncol(output)))

  tmp.env <- new.env()
  assign(
    paste0("output", id_task),
    output,
    pos = tmp.env)
  assign(
    paste0("args", id_task),
    args,
    pos = tmp.env)
  save(
    list = ls(all.names = TRUE, pos = tmp.env), envir = tmp.env, 
    file = paste0(subfolder_data_new, "output", id_task, ".RData"))
  
  rm(output)
  rm(args)
  rm(list = paste0("args", id_task))
}

ls()
rm(list = ls())
