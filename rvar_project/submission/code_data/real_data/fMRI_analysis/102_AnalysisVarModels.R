
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

load("101_Data.RData")

##################################################
##################################################
##################################################
## Function to evaluate J-step forecast performance with
## Root Mean Squared Forecast Error (RMSFE):
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
## INDIVIDUAL VAR MODEL FITTING:
##################################################
##################################################
##################################################

Tsamp <- 100
dfn_list <- df_list %>% 
  lapply(function(x) {
    x %>%
      dplyr::select(-t,-sub) %>%
      mutate_if(is.numeric, scale) %>%
      return()
  })
dfn_list[[1]] %>% apply(2,mean)  
dfn_list[[1]] %>% apply(2,var)  


## Individual VAR model fitting:
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

## Evaluating forecast performance:
B_est_list <- lapply(results,coef) %>% 
  lapply(as.matrix) 
eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) %>%
  apply(1, mean)

## Generating trivial prediction model (persistence model):
B_est_list2 <- rep(ncol(df) - 2, length(df_list)) %>%
  lapply(function(x) {return(diag(x))})
results %>% lapply(function(x) {})


## Relative error of VAR compared to persistence model:
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
## MULTI-VAR MODEL FITTING:
##################################################
##################################################
##################################################

dfn_list <- df_list %>% 
  lapply(function(x) {
    x %>%
      dplyr::filter(t %in% 1:Tsamp) %>%
      dplyr::select(-t,-sub) %>%
      mutate_if(is.numeric, scale) %>%
      return()
  })

## Model fitting:
modelM50 <- multivar::constructModel(data = dfn_list, lassotype = "standard")
fitM50 <- multivar::cv.multivar(modelM50)
fitM50$mats %>% str()
B_estM50_list <- (fitM50$mats)$total

## Relative error of Multi-VAR compared to persistence model:
((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  plot(breaks = 10)
((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  apply(1, mean) %>%
  plot()



## Benefits compared to persistent forecasting.
rfeM50 <- ((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
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
  ggtitle("Multi VAR")

## Benefits compared to persistent forecasting.
rfeVAR <- (eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)
p2 <- rfeVAR %>% 
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

grid.arrange(p1, p2, ncol = 2)



##################################################
##################################################
##################################################
## MOD-VAR WITH SUBJECT-SPECIFIC COMPONENTS ONLY:
## MODEL FITTING
##################################################
##################################################
##################################################

Rcpp::sourceCpp("041_modvar/matrix_fista.cpp")
source("041_modvar/auxfunct.r")
source("041_modvar/adaweights.r")
source("041_modvar/bic.modvar.r")
source("041_modvar/cv.modvar.r")
source("041_modvar/ada.modvar.r")


## Model fitting:
lambdas1  <- 10^(seq(2, -2, length.out = 30))
ratios    <- 10^(seq(2, 0, length.out = 30))
cv.model <- dfn_list %>% 
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
    cv.type = "rolling")


save.image("BrainNetworksGeneData/102_Data.RData")
rm(list = ls())
##################################################
##################################################
load("BrainNetworksGeneData/102_Data.RData")


## Optimization landscape:
par(mfrow = c(1,2),
    mar = c(5.1, 4.1, 4.1, 4.1))
plot(log(cv.model$eval.mat, 10), breaks = 30)
plot(cv.model$eval.mat == min(cv.model$eval.mat))


## Visualizing GCNs per subject:
## Comparison to Multi-VAR and VAR:
par(mfrow = c(2,2))
cv.model$bysubject_coeffs[1:4] %>%
  lapply(plot)
(fitM50$mats)$total[1:4] %>%
  lapply(plot)
B_est_list [1:4] %>%
  lapply(plot)


B_estMOD50_list <- cv.model$bysubject_coeffs


## Relative error of MOD-VAR compared to persistence model:
((eval_forecast(Y_forecast, B_estMOD50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  plot(breaks = 10)
((eval_forecast(Y_forecast, B_estMOD50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  apply(1, mean) %>%
  plot()




## Benefits compared to persistent forecasting.
rfeMOD50 <- ((eval_forecast(Y_forecast, B_estMOD50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
p1 <- rfeMOD50 %>% 
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
  ylim(-0.5,1) +
  ggtitle("MOD-VAR")

## Benefits compared to persistent forecasting.
rfeM50 <- ((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
p2 <- rfeM50 %>% 
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
  ylim(-0.5,1) +
  ggtitle("Multi VAR")

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
  ylim(-0.5,1) +
  ggtitle("VAR")


grid.arrange(p1, p2, p3, ncol = 3)



##################################################
##################################################
## Visualizing results:
dMOD50 <- rfeMOD50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: FI")

dM50 <- rfeM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "Multi-VAR")

dV50 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "VAR")

perf_data <- rbind(dMOD50, dM50, dV50) %>%
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

## ALL OPTIONS OF MOD-VAR MODELS:
Rcpp::sourceCpp("041_modvar/matrix_fista.cpp")
source("041_modvar/auxfunct.r")
source("041_modvar/adaweights.r")
source("041_modvar/bic.modvar.r")
source("041_modvar/cv.modvar.r")
source("041_modvar/ada.modvar.r")


## Fitting MOD-VAR with
## 1. Moderator-only components
lambdas1  <- 10^(seq(0, -5, length.out = 10)) 
ratios    <- 10^(seq(2, -2, length.out = 10)) 
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
    cv.type = "rolling")


## Fitting MOD-VAR with
## 2. Subject-specific components only
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
    cv.type = "rolling")


## Fitting MOD-VAR with
## 3. Both moderator and subject-specific components
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
    cv.type = "rolling")


save.image("BrainNetworksGeneData/102_Data.RData")
rm(list = ls())
##################################################
##################################################
## Visualizations before comparison:
load("BrainNetworksGeneData/102_Data.RData")

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


## How does the evaluation landscape look?
par(mfcol = c(3,2))
plot(log(MO.model$eval.mat, base = 10))
plot(log(IO.model$eval.mat, base = 10))
plot(log(MI.model$eval.mat, base = 10))

plot(MO.model$eval.mat == min(MO.model$eval.mat))
plot(IO.model$eval.mat == min(IO.model$eval.mat))
plot(MI.model$eval.mat == min(MI.model$eval.mat))



##################################################
##################################################
## Visualizing results:
B_estMO50_list <- MO.model$bysubject_coeffs
B_estIO50_list <- IO.model$bysubject_coeffs
B_estMI50_list <- MI.model$bysubject_coeffs



## Relative error of MOD-VAR compared to persistence model:
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
  mutate(type = "Multi-VAR")

dV50 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "VAR")

perf_data <- rbind(dMO50, dIO50, dMI50,  dM50, dV50) %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE))

head(perf_data)


## Performance comparison of all models for T = 100.
ggplot(perf_data, aes(x = factor(StepNo), y = error, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    x = "Step Number",
    y = "Error",
    fill = "Type"
  ) +
  theme_minimal()


