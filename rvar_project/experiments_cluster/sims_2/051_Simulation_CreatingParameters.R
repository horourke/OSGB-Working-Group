#################################################
#################################################
#################################################
##
## In the following document, we introduce the
## functions that run the simulations and
## estimate and turn them into manageable outputs
##
#################################################
#################################################

#################################################
#################################################
#################################################
#################################################
#################################################
#################################################

#################################################
#################################################
## CreateParameters:
##    Function that, given a particular simulation
##    ID id_task, returns the simulation parameters
##    corresponding to such task ID. 
##
##  INPUTS
##    id_task  : ID of the current simulation task.
##                If id_task==0, is a small debugging example.
##                If id_task %in% 1:108, simulations with
##                  different choices of p,T,n, ph, pnh, etc.
##    runtype  : numeric. determines if the run is for
##                runtype== 1: reduced size experiment (2 repetitions)
##                runtype== 1: full size simulation (10 repetitions).
##    
##  OUTPUT:
##    entry_min : minimum entry magnitude
##    entry_max : maximum entry magnitude
##    prob_c    : proportion of common entries.
##    prob_tot  : total proportion of non-zero entries.
##    nsim      : no of simulation repetitions.
##    sigma2    : variance o VAR error term.
##    N         : No. of individuals
##    T         : timepoints per individual.
##    p         : covariate dimension
##    d         : Time series dimension
##
##
CreateParameters <- function(id_task, runtype = c(1, 2, 3)) {
  ## id_task = 0 
  ##  corresponds to a reduced experiment that is useful
  ##  for debugging the code.
  if(id_task == 0) {
    args <- list(
    running_days  = 2,
    entry_min     = 0.1,        ## entry_min : minimum entry magnitude
    entry_max     = 0.9,        ## entry_max : maximum entry magnitude
    
    g_sd          = 0.3,        ## g_sd  : 
    u_min         = 0.5,        ## u_min : minimum entry of exogenous X for type = "unif".
    u_max         = 1,          ## u_max : maximum entry of exogenous X for type = "unif".
    type          = "unif",     ## type  : distribution of exogenous variables.
    nz_x_prob     = 075,        ## nz_x_prob : proportion of entries in exogenous
                                ##              data matrix X with non-zero values.
    signed        = TRUE,       ## signed    : are entries signed or all positive?


    signed        = TRUE,       ## signed   : are entries signed or all positive?
    prob_c        = 0.5,        ## prob_c   : proportion of common entries.
    prob_tot      = 0.1,       ## prob_tot : total proportion of non-zero entries.

    nsim          = 2,          ## nsim     : no of simulation repetitions.
    sigma2        = 0.05,       ## sigma2   : variance o VAR error term.
    N             = 5,          ## N        : No. of individuals
    T             = 30,         ## T        : timepoints per individual.
    p             = 2,          ## p        : covariate dimension
    d             = 5)          ## d        : Time series dimension

    ## Add specific parameters to table:
    args$prob_phi0 <- args$prob_tot * args$prob_c
    args$prob_phip <- args$prob_tot * (1 - args$prob_c)
  
    args$phi0_min <- args$entry_min
    args$phi0_max <- args$entry_max
    args$phip_min <- args$entry_min
    args$phip_max <- args$entry_max
  
    args$vmax <- args$sigma2
    args$vmin <- 0.5 * args$sigma2

    args$id_task <- 0
    
    return(args)
  }
  ## id_task in 1-288
  ##  corresponds to the parameters used for our 
  ##  systematic simulations.

  ## TABLE OF ALL PARAMETER COMBINATIONS.
  sim_par_table <- expand.grid(
    running_days  = 2,
    entry_min     = 0.1,                            ## entry_min : minimum B0,B1,...,Bp entry magnitude
    entry_max     = 0.9,                            ## entry_max : maximum B0,B1,...,Bp entry magnitude
    
    g_sd          = 0.3,                            ## g_sd  : 
    u_min         = 0.4,                            ## u_min : minimum entry of exogenous X for type = "unif".
    u_max         = 1,                              ## u_max : maximum entry of exogenous X for type = "unif".
    type          = "unif",                         ## type  : distribution of exogenous variables.
    nz_x_prob     = c(0.75),                        ## nz_x_prob : proportion of entries in exogenous
                                                    ##              data matrix X with non-zero values.
    signed        = c(TRUE),                        ## signed    : are entries signed or all positive?


    prob_c        = c(2/3, 1/3),                    ## prob_c   : proportion of common entries.
    prob_tot      = 0.05,                           ## prob_tot : total proportion of non-zero entries.

    nsim          = ifelse(runtype <= 2, 2, 5),    ## nsim     : no of simulation repetitions.
    sigma2        = c(0.05, 0.1),                  ## sigma2   : variance o VAR error term.
    N             = c(50, 100),                    ## N        : No. of individuals
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
  
  
  ## Function returns row of index id_task.
  args <- sim_par_table[id_task, ]
  args_list <- list()
  for (i in 1:ncol(args)) {
    if (class(args[, i]) == "factor") {
      args_list[[i]] <- as.character(args[1,i])
    } else args_list[[i]] <- args[1,i]
  }
  names(args_list) <- colnames(args)

  args_list$id_task <- id_task
  
  return(args_list)

}

