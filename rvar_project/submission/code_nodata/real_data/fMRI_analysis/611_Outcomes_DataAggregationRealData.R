
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


###################### Creating folders:
subfolder_new        <- paste0("600_AggregatedDataReal_v3/")
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
    abrev_name     = "VAR"),

  list(
    main_dir       = "200_gimme/",
    run_index      = 1,
    abrev_name     = "GIMME"),
    
  list(
    main_dir       = "300_mvar/",
    run_index      = 1,
    abrev_name     = "M-VAR"),

  list(
    main_dir       = "400_modvar/",
    run_index      = 3,
    abrev_name     = "MOD-VAR"),

  list(
    main_dir       = "500_modvar1sd/",
    run_index      = 2,
    abrev_name     = "MOD-VAR (1SD)")
  )

###########################
## LOOP OVER ALL 9 sample sizes
for (id_task in 1:9) {
  print(paste("XXXXXXXXXXXXXXXX ID-TASK", id_task))

  output <- NULL


  ###########################
  ## Merge all data, from all methods into output.
  for (method_ind in c(1,3,4,5)) {
    main_dir     <- run_info[[method_ind]]$main_dir
    run_index    <- run_info[[method_ind]]$run_index
    runtype_name <- "experiments"
    abrev_name   <- run_info[[method_ind]]$abrev_name

    load(paste0(
      main_dir, runtype_name, run_index,
      "/data/output", id_task, ".RData"))

    ## Add id-task info:
    output_temp <- get(paste0("output", id_task))
    args_temp   <- get(paste0("args", id_task))
      
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
    rm(list = paste0("output", id_task))
  }
  
  ###########################
  ## Saving data:
  print(paste("Saving data: Task-ID", id_task))
  print(paste0(subfolder_data_new, "output", id_task, ".RData"))
  
  tmp.env <- new.env()
  assign(
    paste0("output", id_task),
    output,
    pos = tmp.env)
  save(
    list = ls(all.names = TRUE, pos = tmp.env), envir = tmp.env, 
    file = paste0(subfolder_data_new, "output", id_task, ".RData"))
  
  rm(output)
}

ls()
rm(list = ls())
