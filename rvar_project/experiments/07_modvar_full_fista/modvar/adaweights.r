
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
adaweights <- function(modvar_fit, alpha, inf = 1e10, thr = 1e(-5)) { 
    Boptt <- t(modar_fit$opt_coeffs)

    apply(
        Boptt, c(1,2), 
        function(x) {
            w <- 0
            if (abs(x) > thr) { w <- exp(-alpha * log(abs(x)))
            } else inf
            return(w)
        } ) %>% return()
    
}
