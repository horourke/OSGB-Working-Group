#################################################
#################################################
#################################################
##
## In the following document, we introduce the
## script that has to be run to get the 
## simulation outputs.
##
#################################################
#################################################
#################################################


#################################################
## Sourcing:
#################################################
source("002_requirements_lite.R"); search()
source("003_Generating_ModVarData.R")

Rcpp::sourceCpp("040_modvar/matrix_fista.cpp")
source("040_modvar/auxfunct.r")
source("040_modvar/adaweights.r")
source("040_modvar/bic.modvar.r")
source("040_modvar/cv.modvar.r")
source("040_modvar/ada.modvar.r")

source("051_Simulation_CreatingParameters.R")
source("052_Simulation_EvaluatingMeasures.r")

index  <- 1
main_folder <- "400_modvar/"

#################################################
#################################################
## Step 1: Read imputs:
print("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("#-----------------------Reading inputs")

input                <- commandArgs(trailingOnly = TRUE)
id_original          <- as.numeric(input[1])
id_task              <- id_original %/% 10
id_microrun          <- id_original %% 10
runtype              <- as.numeric(input[2])

args <- CreateParameters(id_task, runtype)
print("#----------------------Verifying Inputs")
print(paste("id_original:", id_original))
print(paste("id_task:", id_task))
print(paste("id_microrun:", id_microrun))
print(args)

## Load Sim-Function:
if (runtype == 1) {
       source(paste0(main_folder, index, "01_SimulationFunctionPretrain.R"))
} else source(paste0(main_folder, index, "02_SimulationFunction.R"))



#################################################
#################################################
## Step 2: Setup run numbers and create folders
print("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("#-----------------------Creating Folders")

## runtype = 1 (experiment), runtype = 2 (full simulation)
runtype_name         <- c("pretrainings","experiments", "outputs")[runtype]
subfolder            <- paste0(main_folder, runtype_name, index, "/")
subfolder_data       <- paste0(subfolder, "data/")
subfolder_logs       <- paste0(subfolder, "logs/")
subfolder_plots       <- paste0(subfolder, "plots/")

if (!dir.exists(subfolder_data)) {
       dir.create(subfolder_data)
}
if (!dir.exists(subfolder_logs)) {
       dir.create(subfolder_logs)
}
if (!dir.exists(subfolder_plots)) {
       dir.create(subfolder_plots)
}

#################################################
#################################################
## Step 2: Run simulation scenarios:
print("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("#------------------Running simulations")
set.seed(id_original)
assign(paste0("output", id_task), FullSimulation(args = args, microrun = id_microrun))

#################################################
#################################################
## Step 3: Save SIM outputs
print("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("#-----------------------Saving Outputs")
tmp_env_sim <- new.env()
assign(paste0("output", id_task,"_", id_microrun),
       get(paste0("output", id_task)),
       pos = tmp_env_sim)
assign(paste0("args", id_task),
       args,
       pos = tmp_env_sim)

print(paste0(subfolder_data, "output", id_task,"_", id_microrun, ".RData"))
save(
       list = ls(all.names = TRUE, pos = tmp_env_sim), envir = tmp_env_sim,
       file = paste0(subfolder_data, "output", id_task,"_", id_microrun, ".RData"))

warnings()