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
    signed        = TRUE,       ## signed   : are entries signed or all positive?
    
    prob_c        = 0.5,        ## prob_c   : proportion of common entries.
    prob_tot      = 0.05,       ## prob_tot : total proportion of non-zero entries.

    nsim          = 2,          ## nsim     : no of simulation repetitions.
    sigma2        = 0.05,       ## sigma2   : variance o VAR error term.
    N             = 5,          ## N        : No. of individuals
    T             = 30,         ## T        : timepoints per individual.
    p             = 2,          ## p        : covariate dimension
    d             = 5)          ## d        : Time series dimension

    args$id_task <- 0
    
    return(args)
  }
  ## id_task in 1-288
  ##  corresponds to the parameters used for our 
  ##  systematic simulations.

  ## TABLE OF ALL PARAMETER COMBINATIONS.
  sim_par_table <- expand.grid(
    running_days  = 2,
    entry_min     = 0.1,                            ## entry_min : minimum entry magnitude
    entry_max     = 0.9,                            ## entry_max : maximum entry magnitude
    signed        = c(TRUE, FALSE),                 ## signed   : are entries signed or all positive?
    
    prob_c        = c(2/3, 1/2, 1/3),               ## prob_c   : proportion of common entries.
    prob_tot      = 0.05,                           ## prob_tot : total proportion of non-zero entries.

    nsim          = ifelse(runtype <= 2, 2, 10),    ## nsim     : no of simulation repetitions.
    sigma2        = c(0.05, 0.1, 0.2),              ## sigma2   : variance o VAR error term.
    N             = c(10, 20),                      ## N        : No. of individuals
    T             = c(30, 50, 100),                 ## T        : timepoints per individual.
    p             = c(2, 5),                        ## p        : covariate dimension
    d             = c(10, 20, 50))                  ## d        : Time series dimension

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

