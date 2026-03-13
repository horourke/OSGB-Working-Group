
#################################################
#################################################
#################################################
#################################################
## bic.solve_rvar_glmnet_vectorized: 
##    For data Xlist and Y, perform cross-validation
##    and model selection for the R-VAR model. 
##
##  INPUTS
##    
##
eval_msr <-function(B_true_list, B_est_list, Y_forecast, range, horizon) {

    l0_est <- mapply(
        function(a, b) {
            a_sp <- (abs(a) > 1e-10)
            b_sp <- (abs(b) > 1e-10)
            return(sum(b_sp != 0))
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()

    l0_true <- mapply(
        function(a, b) {
            a_sp <- (abs(a) > 1e-10)
            b_sp <- (abs(b) > 1e-10)
            return(sum(a_sp != 0))
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()

    TPR <- mapply(
        function(a, b) {
            a_sp <- (abs(a) > 1e-10)
            b_sp <- (abs(b) > 1e-10)
            return(sum(a_sp != 0 & b_sp != 0) / (max(sum(a_sp != 0), 1)))
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()

    FPR <- mapply(
        function(a, b) {
            a_sp <- (abs(a) > 1e-10)
            b_sp <- (abs(b) > 1e-10)
            return(sum(a_sp == 0 & b_sp != 0) / (max(sum(a_sp == 0), 1)))
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()

    l1 <- mapply(
        function(a, b) {
            return(sum(abs(a - b)))
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()

    Fr <- mapply(
        function(a, b) {
            return(sqrt(sum(abs(a - b)^2)))
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()
    

    ## Sensitivity / Specificity
    sens <- mapply(
        function(a, b) {
            a_sp <- (abs(a) > 1e-10)
            b_sp <- (abs(b) > 1e-10)
            return(sum(a_sp * b_sp) / sum(a_sp)  )
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()

    spec <- mapply(
        function(a, b) {
            a_sp <- (abs(a) < 1e-10)
            b_sp <- (abs(b) < 1e-10)
            return(sum(a_sp * b_sp) / sum(a_sp)  )
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()

    ## Forecast performance:
    forecast <- eval_forecast(Y_forecast, B_est_list, range, horizon)

    return(c(
        l0_est = l0_est, l0_true = l0_true, 
        TPR = TPR, FPR = FPR, l1 = l1, Fr = Fr,
        sens = sens, spec = spec, 
        forecast))
}

eval_forecast <- function(Y_forecast, B_est_list, range, horizon) {
    
    n <- length(Y_forecast)
    
    msfe <- numeric(horizon)
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
        msfe[h] <- mean(bysubj_msfe)
        
    }
    
    names(msfe) <- paste0("msfe_step", 1:h)
    
    return(msfe)

}
