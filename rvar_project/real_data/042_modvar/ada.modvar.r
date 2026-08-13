
###########################################################
###########################################################
## cv.modvar:
##  Solver that performs cross validation for the selection
##  of tuning parameters for our proposed MOD-VAR.
##  It estimates the CV-MSFE and selects the combination of
##  tuning parameters that minimizes the error.
##  Allows for rolling or by-subject cross validation.
##
##  INPUTS:
##      Ylist       : list of length n, containing matrices of dimensions 
##                      $Tk \times d$, with the time series data for 
##                      each subject.
##      X           : numeric matrix of dimensions $n x p$. The $k$-th 
##                      row X[k,] must contain the moderator variables of
##                      subject $k$, associated with the $k$-th matrix
##                      in Ylist[[k]].
##      lambdas1    : decreasing numeric vector. Penalizes the joint
##                      structure of the time series.
##      ratios      : decreasing numeric vector. Contains the ratio of 
##                      lambda2 / lambda1. Determines the magnitude of 
##                      of penalty for the moderators and multi-var effects.
##      multi       : if multi = TRUE, we add a fully idiographic
##                      component to each subject, independent from
##                      the moderator effects. 
##      nfolds      : number of folds. Used only if multi = FALSE and
##                      cv.type = "bysubject".
##
##  OUTPUTS:
##      lambda1           : vector with lambda1 values to explore with BIC.
##      ratios            : vector with values of ratio explored with BIC.
##      eval.type         : Type of measure used to perform model selection. Can
##                              take values "cv.rolling", "cv.bysubject" and "BIC".
##      eval.mat          : matrix of model evaluation for each combination of 
##                              lambdas1 and ratios.
##      lambda1_opt       : optimal value of lambda1 according to the BIC procedure.
##      ratio_opt         : optimal value of ratio according to our BIC procedure.
##      lambda2_opt       : optimal value of lambda2 according to the BIC procedure.
##      opt_coeffs        : matrix (d x q) containing optimal coefficients of MOD-VAR
##                          model selected by BIC.
##      joint_coeffs      : (d x d) matrix containing optimal joint autoregressive 
##                              effects across all subjects.
##      moderator_coeffs  : list of length p, containing (d x d) matrices, which
##                              correspond to the autoregressive moderator effects.
##      idiographic_coeffs: list of length n, containing (d x d) matrices, which
##                              correspond to fully idiographic autoregressive effects
##      bysubject_coeffs  : list of length n, containing the aggregregated autoregressive 
##                              effects for each subject. 
##
ada.modvar <- function(
    Ylist, X, 
    lambdas1, 
    ratios, 
    weights = NULL,
    multi = FALSE,
    type = c("bic","cv.rolling", "cv.bysubject"), ## "bysubject" only valid if "multi = FALSE"
    nfolds = NULL,                        ## Only needed if multi = FALSE, and cv.type = "bysubject"
    bic.prop = 0.90,
    ada.alpha = 1
    ) { 
    
    n <- length(Ylist)
    p <- ncol(X)
    d <- ncol(Ylist[[1]])


    if (type == "bic") {
        modvar1 <- bic.modvar(
            Ylist = Ylist, X = X, 
            lambdas1 = lambdas1, 
            ratios = ratios, 
            weights = NULL,
            multi = multi,
            alpha = bic.prop)

        ada.weight <- adaweights(
            modvar_fit = modvar1, 
            p = p,
            multi = FALSE, 
            alpha = ada.alpha, 
            inf = 1e10, 
            thr = 1e-5)

        modvar_fit <- bic.modvar(
            Ylist = Ylist, X = X,
            lambdas1 = lambdas1,
            ratios = ratios,
            weights = ada.weight,
            multi = multi,
            alpha = bic.prop)

        return(modvar_fit)

    } else if (type == "cv.rolling") {
        cv.type = "rolling"
    } else if (type == "cv.bysubject") {
        cv.type = "bysubject"
    }

        
    modvar1 <- cv.modvar(
        Ylist = Ylist, X = X, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = NULL,
        multi = multi,
        cv.type = cv.type,
        nfolds = nfolds)

    ada.weight <- adaweights(
        modvar_fit = modvar1, 
        p = p,
        multi = FALSE, 
        alpha = ada.alpha, 
        inf = 1e10, 
        thr = 1e-5)

    modvar_fit <- cv.modvar(
        Ylist = Ylist, X = X, 
        lambdas1 = lambdas1, 
        ratios = ratios, 
        weights = ada.weight,
        multi = multi,
        cv.type = cv.type,
        nfolds = nfolds)

    return(modvar_fit)


}





