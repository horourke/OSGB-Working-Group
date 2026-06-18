source("002_requirements_lite.R"); search()
source("003_Generating_ModVarData.R")

Rcpp::sourceCpp("../040_modvar/matrix_fista.cpp")
source("../040_modvar/auxfunct.r")
source("../040_modvar/adaweights.r")
source("../040_modvar/bic.modvar.r")
source("../040_modvar/cv.modvar.r")
source("../040_modvar/ada.modvar.r")

source("051_Simulation_CreatingParameters.R")

#################################################
#################################################
## Step 1: Read imputs:
id_task <- 5
runtype <- 2
args <- CreateParameters(id_task, runtype)


seed_val <- 1000* id_task + runtype + 15
set.seed(seed_val)

args$d <- 5
args$n <- 5
args$prob_phi0 <- 0.1
args$prob_phip <- 0
args$prob_delta <- 0.1
args$T <- 50

# Generate Phi0, Phi1, ... Phip parameters:
modvar_model_pars <- generate_modvar_pars(
    args$d, args$p, args$n,
    args$prob_phi0, args$prob_phip, args$prob_delta, 
    args$phi0_min, args$phi0_max,
    args$phip_min, args$phip_max,
    args$delta_min, args$delta_max,
    args$vmin, args$vmax, signed = args$signed)
phi_list <- c(
    modvar_model_pars$phi0, 
    modvar_model_pars$phip_list) ## Concatenated phi0 and phip's.

# Generate exogenous data.
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

Y_list <- lapply(
    data$Y_list,
    function(Y, Ti) {return(Y[1:Ti, ])},
    Ti = args$T)
      
Y_forecast <- lapply(
    data$Y_list,
    function(Y, Ti) {return(Y[-(1:Ti), ])},
    Ti = args$T)


pdf("Coeff_pop.pdf", 8, 8)
layout(matrix(c(1:12), ncol = 4, byrow = FALSE))
modvar_model_pars$phi0 %>% abs() %>% plot( )
modvar_model_pars$delta_list %>% lapply(abs) %>% lapply(plot)
dev.off()

####################################################################################
####################################################################################
####################################################################################
######## BIC.MOD-VAR
{

    W <- build_design(Y_list, X = NULL, multi = TRUE)
    
    pdf("build_mat.pdf", 8, 8)
    plot(W$Wmat, breaks = 15)
    dev.off()

    
    ## 
    lambdas1  <- 10^(seq(3, -3, length.out = 20))
    ratios    <- 10^(seq(3, -3, length.out = 20))
    bic.model <- bic.modvar(
        Ylist = Y_list,
        X = NULL, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        multi = TRUE,
        alpha = 0.90)



    pdf("BIC_fit.pdf", 5, 5)
    layout(matrix(1, ncol = 1))
    plot(log(bic.model$eval.mat - min(bic.model$eval.mat) + 0.01),
         breaks = 15)
    dev.off()


    #########################
    #########################
    ### Visualizing:

    pdf("Coeff_fit.pdf", 8, 8)
    layout(matrix(c(1:12), ncol = 4, byrow = FALSE))
    plot(modvar_model_pars$phi0 != 0)
    modvar_model_pars$delta_list[1:2] %>% lapply(function(x) {plot(x != 0)})
    plot(bic.model$joint_coeffs != 0)
    bic.model$idiographic_coeffs[1:2] %>% lapply(function(x) {plot(x != 0)})

    modvar_model_pars$phi0 %>% abs() %>% plot( )
    modvar_model_pars$delta_list[1:2] %>% lapply(abs) %>% lapply(plot)
    bic.model$joint_coeffs %>% abs() %>% plot()
    bic.model$idiographic_coeffs[1:2] %>% lapply(abs) %>% lapply(plot)
    dev.off()


}

{
    ## 
    lambdas1  <- 10^(seq(3, -3, length.out = 20))
    ratios    <- 10^(seq(3, -3, length.out = 20))
    cv.model <- cv.modvar(
       Ylist = Y_list, 
       X = NULL, 
       lambdas1 = lambdas1, 
       ratios = ratios, 
       multi = TRUE,
       cv.type = "rolling")


    pdf("CV_fit.pdf", 5, 5)
    layout(matrix(1, ncol = 1))
    plot(log(cv.model$eval.mat - min(cv.model$eval.mat) + 0.01),
         breaks = 15)
    dev.off()

    #########################
    #########################
    ### Visualizing:

    pdf("CV_Supp_fit.pdf", width = 10,  height = 10)
    layout(matrix(c(1:24), ncol = 6, byrow = TRUE))
    plot(modvar_model_pars$phi0 != 0)
    modvar_model_pars$delta_list %>% lapply(function(x) {plot(x != 0)})
    plot(cv.model$joint_coeffs != 0)
    cv.model$idiographic_coeffs %>% lapply(function(x) {plot(x != 0)})

    modvar_model_pars$phi0 %>% abs() %>% plot( )
    modvar_model_pars$delta_list %>% lapply(abs) %>% lapply(plot)
    cv.model$joint_coeffs %>% abs() %>% plot()
    cv.model$idiographic_coeffs %>% lapply(abs) %>% lapply(plot)
    dev.off()



}


####################################################################################
####################################################################################
####################################################################################
######## adaBIC.MOD-VAR
{
    ## 
    print(paste0("Step ", sim_ind,".1: adaBIC.MOD-VAR"))
    start_time                  <- Sys.time()
    lambdas1  <- 10^(seq(1, -5, length.out = 20))
    ratios    <- 10^(seq(0, -3, length.out = 20))

    W.ada <- adaweights(
        cv.model, 
        p = 0,
        multi = TRUE, 
        alpha = 1.5, 
        thr = 1e-5)

    ada.model <- cv.modvar(
        Ylist = Y_list,
        X = NULL, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = W.ada,
        multi = TRUE,
        cv.type = "rolling")
    end_time  <- Sys.time()


    pdf("ADA_fit.pdf", 5, 5)
    layout(matrix(1, ncol = 1))
    plot(log(ada.model$eval.mat - min(ada.model$eval.mat) + 0.01),
         breaks = 15)
    dev.off()

    #########################
    #########################
    ### Visualizing:

    pdf("ADA_Supp_fit.pdf", width = 10,  height = 10)
    layout(matrix(c(1:24), ncol = 6, byrow = TRUE))
    plot(modvar_model_pars$phi0 != 0)
    modvar_model_pars$delta_list %>% lapply(function(x) {plot(x != 0)})
    plot(ada.model$joint_coeffs != 0)
    ada.model$idiographic_coeffs %>% lapply(function(x) {plot(x != 0)})

    modvar_model_pars$phi0 %>% abs() %>% plot( )
    modvar_model_pars$delta_list %>% lapply(abs) %>% lapply(plot, breaks = 15)
    ada.model$joint_coeffs %>% abs() %>% plot()
    ada.model$idiographic_coeffs %>% lapply(abs) %>% lapply(plot, breaks = 15)
    dev.off()

}

####################################################################################
####################################################################################
####################################################################################
######## adaBIC.MOD-VAR
{
    ## 
    print(paste0("Step ", sim_ind,".1: adaBIC.MOD-VAR"))
    start_time                  <- Sys.time()
    lambdas1  <- 10^(seq(1, -5, length.out = 20))
    ratios    <- 10^(seq(3, -3, length.out = 20))

    W.ada <- adaweights(
        bic.model, 
        p = 0,
        multi = TRUE, 
        alpha = 1.5, 
        thr = 1e-5)

    ## Exploring ADA weights:
    pdf("ADA_weights.pdf", 5, 8)
    layout(matrix(1:3, ncol = 1))
    plot(abs(bic.model$opt_coeffs), breaks = 15)
    plot(t(log(W.ada)), breaks = 15)
    plot(t(W.ada), breaks = 15)
    dev.off()
    
    ada.model <- cv.modvar(
        Ylist = Y_list,
        X = NULL, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = W.ada,
        multi = TRUE,
        cv.type = "rolling")
    end_time  <- Sys.time()


    pdf("ADA_fit.pdf", 5, 5)
    layout(matrix(1, ncol = 1))
    plot(log(ada.model$eval.mat - min(ada.model$eval.mat) + 0.01),
         breaks = 15)
    dev.off()

    #########################
    #########################
    ### Visualizing:

    pdf("ADA_Supp_fit.pdf", width = 10,  height = 10)
    layout(matrix(c(1:24), ncol = 6, byrow = TRUE))
    plot(modvar_model_pars$phi0 != 0)
    modvar_model_pars$delta_list %>% lapply(function(x) {plot(x != 0)})
    plot(ada.model$joint_coeffs != 0)
    ada.model$idiographic_coeffs %>% lapply(function(x) {plot(x != 0)})

    modvar_model_pars$phi0 %>% abs() %>% plot( )
    modvar_model_pars$delta_list %>% lapply(abs) %>% lapply(plot, breaks = 15)
    ada.model$joint_coeffs %>% abs() %>% plot()
    ada.model$idiographic_coeffs %>% lapply(abs) %>% lapply(plot, breaks = 15)
    dev.off()

}


####################################################################################
####################################################################################
####################################################################################
######## adaBIC.MOD-VAR
{
    ## 
    print(paste0("Step ", sim_ind,".1: adaBIC.MOD-VAR"))
    start_time                  <- Sys.time()
    lambdas1  <- 10^(seq(1, -5, length.out = 20))
    ratios    <- 10^(seq(0, -3, length.out = 20))

    W.ada <- adaweights(
        ada.model, 
        p = 0,
        multi = TRUE, 
        alpha = 1.5, 
        thr = 1e-5)

    ada.model <- cv.modvar(
        Ylist = Y_list,
        X = NULL, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = W.ada,
        multi = TRUE,
        cv.type = "rolling")
    end_time  <- Sys.time()


    pdf("ADA_fit.pdf", 5, 5)
    layout(matrix(1, ncol = 1))
    plot(log(ada.model$eval.mat - min(ada.model$eval.mat) + 0.01),
         breaks = 15)
    dev.off()

    #########################
    #########################
    ### Visualizing:

    pdf("ADA_Supp_fit.pdf", width = 10,  height = 10)
    layout(matrix(c(1:24), ncol = 6, byrow = TRUE))
    plot(modvar_model_pars$phi0 != 0)
    modvar_model_pars$delta_list %>% lapply(function(x) {plot(x != 0)})
    plot(ada.model$joint_coeffs != 0)
    ada.model$idiographic_coeffs %>% lapply(function(x) {plot(x != 0)})

    modvar_model_pars$phi0 %>% abs() %>% plot( )
    modvar_model_pars$delta_list %>% lapply(abs) %>% lapply(plot, breaks = 15)
    ada.model$joint_coeffs %>% abs() %>% plot()
    ada.model$idiographic_coeffs %>% lapply(abs) %>% lapply(plot, breaks = 15)
    dev.off()

}
