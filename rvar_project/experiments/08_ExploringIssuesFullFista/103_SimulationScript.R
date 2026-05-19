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


#################################################
#################################################
## Step 1: Generating data:
id_task <- 10
runtype <- 2
args <- CreateParameters(id_task, runtype)

## Generate PHI
modvar_model_pars <- generate_modvar_pars(
       args$d, args$p, args$n,
       args$prob_phi0, args$prob_phip, args$prob_delta, 
       args$phi0_min, args$phi0_max,
       args$phip_min, args$phip_max,
       args$delta_min, args$delta_max,
       args$vmin, args$vmax, signed = args$signed)

# Generate exogenous data X.
X  <-  simulate_exogenous_vars(
       p = args$p, n = args$n, type = args$type,
       u_min = args$u_min, u_max = args$u_max, g_sd = args$g_sd,
       signed = args$signed, nz_x_prob = args$nz_x_prob) 

# Generate R-VAR Y data.
data <- simulate_modvar1(
       modvar_pars1 = modvar_model_pars, 
       X = X, 
       n = args$n, 
       Ti = rep(args$T +  2 * (args$range + args$horizon), args$n))

# Setup correct format for Y data.
Y_list <- lapply(
       data$Y_list,
       function(Y, Ti) {return(Y[1:Ti, ])},
       Ti = args$T)
      
Y_forecast <- lapply(
       data$Y_list,
       function(Y, Ti) {return(Y[-(1:Ti), ])},
       Ti = args$T)


#########################
#########################
### Visualizing:
str(args)
str(modvar_model_pars)

options(device = "windows")
layout(matrix(c(1,2,3), ncol = 1, byrow = FALSE))
plot(modvar_model_pars$phi0)
modvar_model_pars$phip_list %>% lapply(plot)






#################################################
#################################################
## Evaluating BIC
lambdas1  <- 10^(seq(3, -3, length.out = 20))
ratios    <- 10^(seq(3, -3, length.out = 20))
bic.model <- bic.modvar(
       Ylist = Y_list,
       X = X, 
       lambdas1 = lambdas1, 
       ratios = ratios, 
       multi = FALSE,
       alpha = 0.90)

str(bic.model)

layout(matrix(1, ncol = 1))
plot(log(bic.model$eval.mat - min(bic.model$eval.mat) + 0.01))


#########################
#########################
### Visualizing:
options(device = "windows")
layout(matrix(c(1:12), ncol = 4, byrow = FALSE))
plot(modvar_model_pars$phi0 != 0)
modvar_model_pars$phip_list %>% lapply(function(x) {plot(x != 0)})
plot(bic.model$joint_coeffs != 0)
bic.model$moderator_coeffs %>% lapply(function(x) {plot(x != 0)})

modvar_model_pars$phi0 %>% abs() %>% plot( )
modvar_model_pars$phip_list %>% lapply(abs) %>% lapply(plot)
bic.model$joint_coeffs %>% abs() %>% plot()
bic.model$moderator_coeffs %>% lapply(abs) %>% lapply(plot)

bic.model$joint_coeffs <- bic.model$joint_coeffs %>%
       {. * (abs(.) > 0.01)}


eval_msr(data$B_list, bic.model$bysubject_coeffs, Y_forecast, args$range, args$horizon)





#################################################
#################################################
## Evaluating BIC
lambdas1  <- 10^(seq(3, -3, length.out = 20))
ratios    <- 10^(seq(3, -3, length.out = 20))
cv.model <- cv.modvar(
       Ylist = Y_list, 
       X = X, 
       lambdas1 = lambdas1, 
       ratios = ratios, 
       multi = FALSE,
       cv.type = "bysubject",
        nfolds = 5)


str(cv.model)

layout(matrix(1:2, ncol = 1))
par(mar = c(5.1, 4.1, 4.1, 4.1))
plot(log(bic.model$eval.mat - min(bic.model$eval.mat) + 0.01), breaks = 10)
plot(log(cv.model$eval.mat - min(cv.model$eval.mat) + 0.01), breaks = 10)
#plot(cv.model$eval.mat <= 1.1* min(cv.model$eval.mat))


#########################
#########################
### Visualizing:
options(device = "windows")
layout(matrix(c(1:12), ncol = 4, byrow = FALSE))
plot(modvar_model_pars$phi0 != 0)
modvar_model_pars$phip_list %>% lapply(function(x) {plot(x != 0)})
plot(cv.model$joint_coeffs != 0)
cv.model$moderator_coeffs %>% lapply(function(x) {plot(x != 0)})

modvar_model_pars$phi0 %>% abs() %>% plot( )
modvar_model_pars$phip_list %>% lapply(abs) %>% lapply(plot)
cv.model$joint_coeffs %>% abs() %>% plot()
cv.model$moderator_coeffs %>% lapply(abs) %>% lapply(plot)

cv.model$joint_coeffs <- cv.model$joint_coeffs %>%
       {. * (abs(.) > 0.05)}




#################################################
#################################################
## Evaluating BIC
W.ada <- adaweights(
       cv.model, 
       p = 2,
       multi = FALSE, 
       alpha = 1, 
       inf = 1e3, 
       thr = 1e-4)

## Exploring ADA weights:
layout(matrix(1:4, ncol = 1))
par(mar = c(5.1, 4.1, 4.1, 8.1))
plot(cv.model$opt_coeffs, breaks = 15)
plot(t(log(W.ada)))
plot(t(W.ada))
plot(abs(cv.model$opt_coeffs) > 1e-5)


ada.model <- bic.modvar(
       Ylist = Y_list,
       X = X, 
       lambdas1 = lambdas1, 
       ratios = 1, 
       weights = W.ada,
       multi = FALSE,
       alpha = 0.90)
dim(t(W.ada))

#########################
#########################
### Visualizing:
options(device = "windows")
layout(matrix(c(1:12), ncol = 4, byrow = FALSE))
plot(modvar_model_pars$phi0 != 0)
modvar_model_pars$phip_list %>% lapply(function(x) {plot(x != 0)})
plot(ada.model$joint_coeffs != 0)
ada.model$moderator_coeffs %>% lapply(function(x) {plot(x != 0)})

modvar_model_pars$phi0 %>% abs() %>% plot( )
modvar_model_pars$phip_list %>% lapply(abs) %>% lapply(plot)
ada.model$joint_coeffs %>% abs() %>% plot()
ada.model$moderator_coeffs %>% lapply(abs) %>% lapply(plot)

