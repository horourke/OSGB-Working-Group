############################################
############################################
## SETTING PARAMETERS:
wd <- getwd()
.libPaths(paste0(wd,"/req_lib"))

library(readr)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(stringr)
library(lazyeval)
library(pROC)

input <- commandArgs(trailingOnly = TRUE)
diag_shift_val <- as.numeric(input[1])
T0_prop_val    <- as.numeric(input[2])


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


  prob_c        = c(0.25, 0.75),                  ## prob_c   : proportion of common entries.
  prob_tot      = 0.1,                            ## prob_tot : total proportion of non-zero entries.

  nsim          = ifelse(
                    runtype == 1, 
                    1, 
                    ifelse(runtype == 2, 2, 10)), ## nsim     : no of simulation repetitions.
  sigma2        = c(0.1),                         ## sigma2   : variance o VAR error term.
  n             = c(10, 20),                      ## n        : No. of individuals
  T             = c(30, 50, 100),                 ## T        : timepoints per individual.
  p             = c(2, 5),                        ## p        : covariate dimension
  d             = c(10, 20, 30))                  ## d        : Time series dimension

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
###################### Processing Data
##################################################################
##################################################################
method_names <- c(
    "var_standard",
    "mvar_standard", "mvar_adaptive",
    "modvar_bic", "modvar_cv_roll", "modvar_cv_bsubj")
method_names_clean <- c(
    "VAR",          
    "M-VAR: Standard", "M-VAR: Adaptive",
    "MOD-VAR: BIC", "MOD-VAR: RWCV", "MOD-VAR: BSCV")


for (d_val in c(10, 20, 30)) {
  
  ##############################
  ##############################
  ## LOADING ALL DATA WITH T0 = P.
  sim_ind_load    <- which(d == d_val) #p == p_val)
  type            <- "all"
  results_dir     <- paste0(subfolder_new, "plots_", type, "/")

  for (sim_ind in sim_ind_load) {
    load(paste0(
      subfolder_new, "data_all/",
      "output", sim_ind, ".RData"))
  }

  ##############################
  ##############################
  ## TRUE POSITIVE RATE OF JIC-HD METHODS:

  ## Merge dataset of JIC-HD-derived data.
  output_merged  <- distinct(bind_rows(mget(ls(pattern = '^output\\d+')))) %>%
    filter(method %in% method_names) %>%
    #mutate(method = str_replace_all(method, setNames(method_names_clean, method_names))) %>%
    mutate(total_N = n * T) %>% 
    select(-n, -T) %>%
    group_by(d, p, total_N, sigma2, prob_c, signed, method) %>%

    summarise(
      meanTPR = mean(TPR),
      sdTPR = sd(TPR),
      meanFPR = mean(FPR),
      sdFPR = sd(FPR),
      meanTime = mean(time),
      sdTime = sd(time),
      meanL1 = mean(l1),
      sdL1 = sd(l1),
      meanFr = mean(Fr),
      sdFr = sd(Fr))

  file_name <- paste0(
  subfolder_plots_new, 
  "521_MainTPR",
  "_d", d_val, 
  ".pdf")
  pdf(file_name, width = 8, height = 5)
  p1 <-  output_merged %>% ggplot(aes(x = total_N, y = meanTPR)) +
    geom_line(aes(col = method, linetype = method), linewidth = 1) +
    geom_ribbon(aes(ymin = meanTPR - sdTPR, ymax = meanTPR + sdTPR, fill = method), alpha = 0.1) +
    geom_hline(yintercept = c(0,1), linetype = 2) +
    facet_grid(p ~ sigma2 + prob_c, scales = "free_x", labeller = label_parsed) +
    theme(legend.position="bottom")
  print(p1)

  p2 <-  output_merged %>% ggplot(aes(x = total_N, y = meanFPR)) +
    geom_line(aes(col = method, linetype = method), linewidth = 1) +
    geom_ribbon(aes(ymin = meanFPR - sdFPR, ymax = meanFPR + sdFPR, fill = method), alpha = 0.1) +
    geom_hline(yintercept = c(0), linetype = 2) +
    facet_grid(p ~ sigma2 + prob_c, scales = "free_x", labeller = label_parsed) +
    theme(legend.position="bottom")
  print(p2)

  p3 <-  output_merged %>% 
    mutate(logmt = log(meanTime, base = 60)) %>%
    
    ggplot(aes(x = total_N, y = logmt)) + 
      geom_line(aes(col = method, linetype = method), linewidth = 1) +
      #geom_ribbon(aes(ymin = meanTime - sdTime, ymax = meanTime + sdTime, fill = method), alpha = 0.1) +
      geom_hline(yintercept = c(0,1,2), linetype = 3) +
      facet_grid(p ~ sigma2 + prob_c, scales = "free_x", labeller = label_parsed) +
      theme(legend.position="bottom")
  print(p3)
    
  dev.off()

  rm(list = grep('^output\\d+', ls(), value = TRUE))
  rm(list = grep("args", ls(), value = TRUE))
  rm(
    output_merged, file_name, p1, p2, p3)

}

