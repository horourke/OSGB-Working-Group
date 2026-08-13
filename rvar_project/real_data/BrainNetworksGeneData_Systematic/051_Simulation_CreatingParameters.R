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
##                If id_task %in% 1:8, simulations with
##                  different choices of p,T,n, ph, pnh, etc.
##    
##  OUTPUT:
##    running_days :
##    range        :
##    horizon      :
##    n            :
##    args$id_task :
##
CreateParameters <- function(id_task) {
  ## id_task = 0 
  ##  corresponds to a reduced experiment that is useful
  ##  for debugging the code.
  if(id_task == 0) {
    args <- list(
      running_days  = 2,
      range         = 200,
      horizon       = 10,
      n             = 25)

    args$id_task <- 0
    return(args)
  }

  ## TABLE OF ALL PARAMETER COMBINATIONS.
  T_vec <- c(20,40,60,80,100,125,150,175,200)
  args <- list(
    running_days  = 2,
    range         = 150,
    horizon       = 10,
    T             = T_vec[id_task])                 ## T        : No. of time points

  args$id_task <- id_task
  
  return(args)

}


