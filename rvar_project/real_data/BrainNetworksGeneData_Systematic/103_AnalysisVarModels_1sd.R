
library(tseries)
library(MTS)
library(ppcor)
library(TSA)
library(forecast)
library(ggplot2)

library(magrittr)
library(readxl)
library(tidyverse)
library(plot.matrix)

## Model fitting packages:
library(glmnet)
library(mvtnorm)
library(multivar)
library(BigVAR)
library(expm)
library(gridExtra)

Rcpp::sourceCpp("042_modvar/matrix_fista.cpp")
source("042_modvar/auxfunct.r")
source("042_modvar/adaweights.r")
source("042_modvar/bic.modvar.r")
source("042_modvar/cv.modvar.r")
source("042_modvar/ada.modvar.r")

load("BrainNetworksGeneData_Systematic/101_Data.RData")

##################################################
##################################################
##################################################
eval_forecast <- function(Y_forecast, B_est_list, range, horizon) {
  
  n <- length(Y_forecast)
  
  #msfe <- numeric(horizon)
  msfe <- matrix(0, ncol = n, nrow = horizon)
  for (h in 1:horizon) {
    
    bysubj_msfe <- numeric(n)
    for (k in 1:n) {
      Y <- Y_forecast[[k]]
      B <- B_est_list[[k]] %^% h
      
      Z <- Y[h + 1:(range - h),]
      X <- Y[1:(range - h),]
      
      err_mat <- Z - X %*% t(B)
      
      bysubj_msfe[k] <- sqrt(mean(err_mat^2))
      
    }
    msfe[h,] <- bysubj_msfe
    
  }
  
  rownames(msfe) <- paste0("msfe_step", 1:h)
  
  return(msfe)
  
}
##################################################
##################################################
##################################################

Tsamp <- 50
dfn_list <- df_list %>% 
  lapply(function(x) {
    x %>%
      dplyr::select(-t,-sub) %>%
      mutate_if(is.numeric, scale) %>%
      return()
  })
dfn_list[[1]] %>% apply(2,mean)  
dfn_list[[1]] %>% apply(2,var)  



models <- list()
results <- list()
Y_forecast <- list()
for(i in 1:length(df_list)) {
  models[[i]] <- df_list[[i]][1:Tsamp,] %>% 
    tibble() %>%
    dplyr::select(-sub,-t) %>% 
    as.matrix() %>%
    BigVAR::constructModel( 
      p = 1,
      gran = c(50,10),
      struct = "Basic",
      cv = "Rolling",
      verbose = TRUE,
      ownlambdas = FALSE,
      model.controls=list(intercept=FALSE),
      linear = FALSE)
  results[[i]] <- BigVAR::cv.BigVAR(models[[i]])
  print(coef(results[[i]]))
  
  Y_forecast[[i]] <- df_list[[i]][-(1:Tsamp),] %>% 
    tibble() %>%
    dplyr::select(-sub,-t) %>% 
    as.matrix() 

}


str(results[[1]])
apply(results[[1]]@InSampMSFE, 2, mean)
apply(results[[1]]@OOSMSFE, 2, mean)
results[[1]]
results[[1]]@LambdaGrid
results[[1]]@InSampMSFE

B_est_list <- lapply(results,coef) %>% 
  lapply(as.matrix) 
eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) %>%
  apply(1, mean)
##RW Out of Sample Loss
## [1] 0.294


B_est_list2 <- rep(ncol(df) - 2, length(df_list)) %>%
  lapply(function(x) {return(diag(x))})
results %>% lapply(function(x) {})



par(mar = c(5.1, 4.1, 4.1, 4.1))
((eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  plot(breaks = 10)


((eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  apply(1, mean) %>%
  plot()


## Benefits compared to persistent forecasting.
((eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
    geom_boxplot() +
    geom_hline(yintercept = 0)
  
## We observe that for the 1st step prediction, the model performs
## somewhat poorly (~1% worse than persistence forecasting), but
## as we increase the time horizon, the estimation improves
## significantly (on average, 20% relative improvement)


##################################################
##################################################
##################################################
##################################################
##################################################
##################################################

dfn_list <- df_list %>% 
  lapply(function(x) {
    x %>%
      dplyr::filter(t %in% 1:50) %>%
      dplyr::select(-t,-sub) %>%
      mutate_if(is.numeric, scale) %>%
      return()
  })

dfn_list[[1]] %>% dim()
dfn_list[[1]] %>% apply(2,mean)  
dfn_list[[1]] %>% apply(2,var)  


modelM50 <- multivar::constructModel(data = dfn_list, lassotype = "standard")
modelAdaM50 <- multivar::constructModel(data = dfn_list, lassotype = "adaptive")
fitM50 <- multivar::cv.multivar(modelM50)
fitAdaM50 <- multivar::cv.multivar(modelAdaM50)

fitM50$mats %>% str()
B_estM50_list <- (fitM50$mats)$total
B_estAdaM50_list <- (fitAdaM50$mats)$total
str(fitM50)

dim(fitM50$MSFE)
str(fitM50)

((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  plot(breaks = 10)

((eval_forecast(Y_forecast, B_estAdaM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  plot(breaks = 10)


((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  apply(1, mean) %>%
  plot()
((eval_forecast(Y_forecast, B_estAdaM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  apply(1, mean) %>%
  plot()



## Benefits compared to persistent forecasting.
rfeM50 <- ((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeAdaM50 <- ((eval_forecast(Y_forecast, B_estAdaM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
p1 <- rfeM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  ggtitle("Multi VAR: Standard")

p2 <- rfeAdaM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  ggtitle("Multi VAR: Adaptive")

## Benefits compared to persistent forecasting.
rfeVAR <- (eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)
p3 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  ggtitle("VAR")



grid.arrange(p1, p2, p3, ncol = 3)

ls()

rm(i, p1, p2, Tsamp)

save.image("BrainNetworksGeneData_Systematic/103_Data.RData")
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
## MOD-VAR with only individual structure:


#lambdas1  <- 10^(seq(2, -2, length.out = 10)) # Best so far!
#ratios    <- 10^(seq(2, 0, length.out = 10)) 

lambdas1  <- 10^(seq(2, -5, length.out = 20)) # Best so far!
ratios    <- 10^(seq(2, -2, length.out = 20)) 

lcdf_norm <- lcdf_norm %>% 
  as.matrix()
colnames(lcdf_norm) <- NULL 
MO.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = lcdf_norm,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = FALSE,
    cv.type = "rolling",
    sparse1sd = TRUE)


IO.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = NULL,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = TRUE,
    cv.type = "rolling",
    sparse1sd = TRUE)


MI.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = lcdf_norm,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = TRUE,
    cv.type = "rolling",
    sparse1sd = TRUE)


save.image("BrainNetworksGeneData_Systematic/103_Data.RData")
rm(list = ls())
##################################################
##################################################
## Visualizations before comparison:
load("BrainNetworksGeneData_Systematic/103_Data.RData")

## Verifying correct fit conditions:
str(MO.model$idiographic_coeffs)
str(MO.model$moderator_coeffs)
str(MO.model$joint_coeffs)

str(IO.model$idiographic_coeffs)
str(IO.model$moderator_coeffs)
str(IO.model$joint_coeffs)

str(MI.model$idiographic_coeffs)
str(MI.model$moderator_coeffs)
str(MI.model$joint_coeffs)

## Are the models having the same sparsity pattern?
par(mfcol = c(3,2))
plot(log(abs(MI.model$joint_coeffs) + 1e-10, 10) > -2)
plot(log(abs(MO.model$joint_coeffs) + 1e-10, 10) > -2)
plot(log(abs(IO.model$joint_coeffs) + 1e-10, 10) > -2)

plot(MI.model$joint_coeffs == MO.model$joint_coeffs)
plot(MI.model$joint_coeffs == IO.model$joint_coeffs)


## How does the evaluation matrix look?
par(mfcol = c(3,2))
plot(log(MO.model$eval.mat, base = 10))
plot(log(IO.model$eval.mat, base = 10))
plot(log(MI.model$eval.mat, base = 10))
which(MO.model$lambda1 == MO.model$lambda1_opt)
which(MO.model$ratios == MO.model$ratio_opt)

plot(MO.model$eval.mat == min(MO.model$eval.mat))
plot(IO.model$eval.mat == min(IO.model$eval.mat))
plot(MI.model$eval.mat == min(MI.model$eval.mat))


MO.model$lambda1_opt
MO.model$ratio_opt


##################################################
##################################################
## Visualizing results:
B_estMO50_list <- MO.model$bysubject_coeffs
B_estIO50_list <- IO.model$bysubject_coeffs
B_estMI50_list <- MI.model$bysubject_coeffs

rfeMO50 <- ((eval_forecast(Y_forecast, B_estMO50_list, range = 100, horizon = 10) -
               eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeIO50 <- ((eval_forecast(Y_forecast, B_estIO50_list, range = 100, horizon = 10) -
               eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeMI50 <- ((eval_forecast(Y_forecast, B_estMI50_list, range = 100, horizon = 10) -
               eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 


rfeM50 <- ((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 

rfeAdaM50 <- ((eval_forecast(Y_forecast, B_estAdaM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 

rfeVAR <- (eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)

#######################
#######################
## Generating plots:
dMO50 <- rfeMO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: MOD")

dIO50 <- rfeIO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: FI")

dMI50 <- rfeMI50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: MOD+FI")

dM50 <- rfeM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "Multi-VAR: Standard")

dAdaM50 <- rfeAdaM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "Multi-VAR: Adaptive")

dV50 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "VAR")

perf_data <- rbind(dMO50, 
                   dIO50, 
                   dMI50,
                   dM50, 
                   dAdaM50, 
                   dV50) %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE))

head(perf_data)

ggplot(perf_data, aes(x = factor(StepNo), y = error, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    x = "Step Number",
    y = "Error",
    fill = "Type"
  ) +
  theme_minimal()





##################################################
##################################################
##################################################
##################################################
##################################################
##################################################


#lambdas1  <- 10^(seq(2, -2, length.out = 10)) # Best so far!
#ratios    <- 10^(seq(2, 0, length.out = 10)) 

lcdf_norm <- lcdf_norm %>% 
  as.matrix()
colnames(lcdf_norm) <- NULL


lambdas1  <- 10^(seq(2, -2, length.out = 20)) # Best so far!
ratios    <- 10^(seq(2, -2, length.out = 20)) 
W.ada <- adaweights(
  modvar_fit = MO.model,
  p = 4,
  multi = FALSE, 
  alpha = 1, 
  thr = 1e-5)
adaMO.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = lcdf_norm, 
    lambdas1 = lambdas1,
    ratios = ratios, 
    multi = FALSE, 
    weights = W.ada,
    cv.type = "rolling",
    sparse1sd = TRUE)




lambdas1  <- 10^(seq(2, -2, length.out = 10)) # Best so far!
ratios    <- 10^(seq(2, -2, length.out = 10)) 
W.ada <- adaweights(
  IO.model, 
  p = 0,
  multi = TRUE, 
  alpha = 1, 
  thr = 1e-5)
adaIO.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = NULL,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = TRUE,
    weights = W.ada,
    cv.type = "rolling",
    sparse1sd = TRUE)


lambdas1  <- 10^(seq(2, -2, length.out = 10)) # Best so far!
ratios    <- 10^(seq(2, -2, length.out = 10))
W.ada <- adaweights(
  MI.model, 
  p = 4,
  multi = TRUE, 
  alpha = 1, 
  thr = 1e-5)
adaMI.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = lcdf_norm,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = TRUE,
    weights = W.ada,
    cv.type = "rolling",
    sparse1sd = TRUE)




save.image("BrainNetworksGeneData_Systematic/103_Data.RData")



## Verifying correct fit conditions:
str(MO.model$idiographic_coeffs)
str(MO.model$moderator_coeffs)
str(MO.model$joint_coeffs)

str(IO.model$idiographic_coeffs)
str(IO.model$moderator_coeffs)
str(IO.model$joint_coeffs)

str(MI.model$idiographic_coeffs)
str(MI.model$moderator_coeffs)
str(MI.model$joint_coeffs)

## Are the models having the same sparsity pattern?
par(mfcol = c(3,2))
plot(log(abs(MI.model$joint_coeffs) + 1e-10, 10) > -2)
plot(log(abs(MO.model$joint_coeffs) + 1e-10, 10) > -2)
plot(log(abs(IO.model$joint_coeffs) + 1e-10, 10) > -2)

plot(MI.model$joint_coeffs == MO.model$joint_coeffs)
plot(MI.model$joint_coeffs == IO.model$joint_coeffs)


## How does the evaluation matrix look?
par(mfcol = c(3,2))
plot(log(MO.model$eval.mat, base = 10))
plot(log(IO.model$eval.mat, base = 10))
plot(log(MI.model$eval.mat, base = 10))

plot(MO.model$eval.mat == min(MO.model$eval.mat))
plot(IO.model$eval.mat == min(IO.model$eval.mat))
plot(MI.model$eval.mat == min(MI.model$eval.mat))


MO.model$lambda1_opt
MO.model$ratio_opt


##################################################
##################################################
## Visualizing results:
B_estMO50_list <- MO.model$bysubject_coeffs
B_estIO50_list <- IO.model$bysubject_coeffs
B_estMI50_list <- MI.model$bysubject_coeffs

B_estAdaMO50_list <- adaMO.model$bysubject_coeffs
B_estAdaIO50_list <- adaIO.model$bysubject_coeffs
B_estAdaMI50_list <- adaMI.model$bysubject_coeffs

rfeAdaMO50 <- ((eval_forecast(Y_forecast, B_estAdaMO50_list, range = 100, horizon = 10) -
                  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
                 eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeAdaIO50 <- ((eval_forecast(Y_forecast, B_estAdaIO50_list, range = 100, horizon = 10) -
                  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
                 eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeAdaMI50 <- ((eval_forecast(Y_forecast, B_estAdaMI50_list, range = 100, horizon = 10) -
                  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
                 eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 


rfeMO50 <- ((eval_forecast(Y_forecast, B_estMO50_list, range = 100, horizon = 10) -
                  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
                 eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeIO50 <- ((eval_forecast(Y_forecast, B_estIO50_list, range = 100, horizon = 10) -
                  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
                 eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeMI50 <- ((eval_forecast(Y_forecast, B_estMI50_list, range = 100, horizon = 10) -
                  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
                 eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 


rfeM50 <- ((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 

rfeAdaM50 <- ((eval_forecast(Y_forecast, B_estAdaM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 

rfeVAR <- (eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)

#######################
#######################
## Generating plots:
dAdaMO50 <- rfeAdaMO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "adaMOD-VAR: MOD")

dAdaIO50 <- rfeAdaIO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "adaMOD-VAR: FI")

dAdaMI50 <- rfeAdaMI50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "adaMOD-VAR: MOD+FI")

dMO50 <- rfeMO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: MOD")

dIO50 <- rfeIO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: FI")

dMI50 <- rfeMI50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: MOD+FI")

dM50 <- rfeM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "Multi-VAR: Standard")

dAdaM50 <- rfeAdaM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "Multi-VAR: Adaptive")

dV50 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step5,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "VAR")

perf_data <- rbind(
  dAdaMO50, dAdaIO50, dAdaMI50,
  dMO50, dIO50, dMI50,
  dM50, dAdaM50, dV50) %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE))

unique(perf_data$type)
head(perf_data)

ggplot(perf_data, aes(x = factor(StepNo), y = error, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    x = "Step Number",
    y = "Error",
    fill = "Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")




##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
adaMO.model$moderator_coeffs
plot(adaMO.model$joint_coeffs)

plot(adaMO.model$joint_coeffs - diag(diag(adaMO.model$joint_coeffs)))
fitadaM50$mats$unique

par(mfrow = c(2,2))
MO.model$moderator_coeffs %>%
  lapply(plot)


head(lcdf_norm)

plot(MO.model$moderator_coeffs[[3]])
