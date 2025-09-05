
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
eval_msr <-function(B_true_list, B_est_list) {

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
            return(sum(a_sp != 0 & b_sp != 0) / (max(sum(b_sp != 0), 1)))
        },
        a = B_true_list, b = B_est_list
    ) %>% unlist() %>% mean()

    FPR <- mapply(
        function(a, b) {
            a_sp <- (abs(a) > 1e-10)
            b_sp <- (abs(b) > 1e-10)
            return(sum(a_sp != 0 & b_sp == 0) / (max(sum(b_sp == 0), 1)))
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
    
    return(c(
        l0_est = l0_est, l0_true = l0_true, 
        TPR = TPR, FPR = FPR, l1 = l1, Fr = Fr))
}