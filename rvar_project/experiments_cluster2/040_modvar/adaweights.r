
###########################################################
###########################################################
## ada.weights:
##  Given a previously fit MOD-VAR model, we calculate new
##  weights for our estimator to perform adaptive estimation.
##
##  INPUTS:
##      modvar_fit  : Object output from cv.modvar or bic.modvar.
##      alpha       : Adaptive exponent parameter.
##      inf         : value of infinity to use for zero entries.
##      thr         : threshold that determines numeric zeros.
##      
##  OUTPUTS:
##      weights     : matrix of adaptive weights.
##
adaweights <- function(
    modvar_fit, 
    p = 0,
    multi = FALSE, 
    alpha, 
    inf = NULL, 
    thr = 1e-5) { 

    n <- length(modvar_fit$idiographic_coeffs)
    d <- ncol(modvar_fit$joint_coeffs)
    q <- ncol(modvar_fit$opt_coeffs)
    p <- ifelse(
        is.null(modvar_fit$moderator_coeffs), 
        0, 
        length(modvar_fit$moderator_coeffs))

    min_nz <- modvar_fit$opt_coeffs[, 1:(d * (1 + p))] %>%
        {(.)[abs(.) > thr]} %>%
        abs() %>% min()
    if (is.null(inf)) {
        inf_val <- (min_nz/2)^(-alpha)
    } else {
        inf_val <- inf
    }
    M1 <- modvar_fit$opt_coeffs[, 1:(d * (1 + p))] %>%
        apply(c(1,2), function(x) {
            entry_adaweight(x, alpha = alpha, thr = thr, inf = inf_val)
        })## What if you use median?

    if(multi) { ## TODO: FIX 
        median_idiographic <- modvar_fit$idiographic_coeffs %>%
        unlist() %>% array(c(d,d,n)) %>%
        apply(c(1,2), FUN = stats::median)
        print(median_idiographic)

        M2 <- modvar_fit$idiographic_coeffs %>%
            lapply(function(x, median) {abs(x - median)},
            median = median_idiographic)
        print(str(M2))

        min_nz <- M2 %>% {.[abs(.) > thr]} %>%
            abs() %>% min()
        if (is.null(inf)){ 
            inf_val <- (min_nz / 2)^(-alpha)
        } else {
            inf_val <- inf  
        }
        M2 <- M2 %>%
            lapply(
                function(x) {
                    apply(
                        x, c(1,2), 
                        FUN = function(x) {
                            entry_adaweight(x, alpha = alpha, thr = thr, inf = inf_val)
                        })
                }
            ) 
        print(str(M2))
        M2 <- do.call(what = cbind, M2)

        M1 <- cbind(M1, M2)
    }   

    output <- t(M1) 

    return(output)
}



entry_adaweight <- function(x, alpha = 1, thr = 1e-5, inf = 10^(10)) {
    w <- 0
    if (abs(x) > thr) { 
        w <- exp(-alpha * log(abs(x)))
    } else w <- inf
    return(w)
}
