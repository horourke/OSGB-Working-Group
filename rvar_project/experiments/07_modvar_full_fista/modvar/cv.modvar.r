
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
##      eval.type         : measure for evaluating our CV performance.
##      eval.mat          : matrix of cross validation forecasting error. rows correspond
##                              to lambdas1 range and columns to ratios)
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
cv.modvar <- function(
    Ylist, X, 
    lambdas1, 
    ratios, 
    weights = NULL,
    multi = FALSE,
    cv.type = c("rolling", "bysubject"), ## "bysubject" only valid if "multi = FALSE"
    nfolds = NULL                        ## Only needed if multi = FALSE, and cv.type = "bysubject"
    ) { 

    ############################
    ## Build design matrix:
    BD <- build_design(Ylist, X, multi)
    Wmat <- BD$Wmat
    Ymat <- BD$Ymat

    n <- length(Ylist)
    q <- ncol(Wmat)
    N <- nrow(Ymat)
    d <- ncol(Ymat)
    p <- ifelse(is.null(X), 0, ncol(X))

    ############################
    ## Construct weights
    if (is.null(weights)) {
        weights <- matrix(1, nrow = q, ncol = d)
    }

    ############################
    ## Building the cross validation indexes:
    nlam1 <- length(lambdas1)
    nrat  <- length(ratios)
    Bcvprev <- array(0, c(q, d, nlam1 * nrat))

    Ts  <- unlist(lapply(Ylist, nrow))
    cv.windows <- NULL
    if (cv.type == "rolling") {
        cv.windows <- build_rolling_windows(Ts)
    } else if (cv.type == "bysubject") {
        if (multi) {
            print("Error: the cv.type \"bysubject\" invalid for multi = TRUE.")
            return(-1)
        } else if (is.null(nfolds)) {
            print("Error: invalid value of nfolds.") ; return(-1)
        } else {
            cv.windows <- build_bysubject_windows(Ts, nfolds)
        }
    } else {print("Error: invalid value of cv.type.") ; return(-1)}

    ############################
    ## Calculate CV.MSFE:
    MSFE <- array(0, c(nlam1, nrat, length(cv.windows[[1]])))
    pb   <- txtProgressBar(1, length(cv.windows[[1]]), style=3)

    for (ind in seq_along(cv.windows[[1]])) {
        
        setTxtProgressBar(pb, ind)

        ## Build cv indexes:
        cv.train    <- cv.windows$cv.train_list[[ind]]
        cv.test     <- cv.windows$cv.test_list[[ind]]

        ## Select warm start:
        Bcvprev <- NULL
        if (ind == 1) {
            Bcvprev <- array(0, c(q, d, nlam1 * nrat))
        } else {
            Bcvprev <- Bout
        }

        ############################
        ## Fit model:
        Wmat.cv <- Wmat[cv.train, ]
        Ymat.cv <- Ymat[cv.train, ]

        Bout  <- weighted_lasso_path(
            X = Wmat.cv,
            Y = Ymat.cv,
            lambda1vec = lambdas1,
            ratiovec = ratios,
            weights = weights,
            max_iter = 2000,
            Bcvprev = Bcvprev,
            tol = 1e-8)

        ############################
        ## Forecasting error:
        Ymat.test <- Ymat[cv.test, ]
        Wmat.test <- Wmat[cv.test, ]

        MSFE[, , ind] <- forecast_mse_batch(
            Xf = Wmat.test,
            Yf = Ymat.test,
            Bflat = Bout,
            nlambda1 = nlam1,
            nratios = nrat)
 
    }
    CV.MSFE <- apply(MSFE, c(1,2), mean)

    ############################
    ## Find minimal lambda1, lambda2:
    opt.inds <- which(CV.MSFE == min(CV.MSFE), arr.ind = TRUE)
    print(opt.inds)

    lambda1_opt <- lambdas1[opt.inds[1,1]]
    ratio_opt <- ratios[opt.inds[1,2]]

    ############################
    ## Fit optimal model:
    Bopt <- weighted_lasso_path(
        X = Wmat,
        Y = Ymat,
        lambda1vec = lambda1_opt,
        ratiovec = ratio_opt,
        weights = weights,
        max_iter = 2000,
        Bcvprev = array(0, c(q, d, nlam1 * nrat)),
        tol = 1e-8)

    dim(Bopt) <- c(q,d)
    Bopt <- t(as.matrix(Bopt))

    ############################
    ## Build outputs:
    joint_coeffs <- Bopt[, 1:d]
    if(p > 0) {
        moderator_coeffs <- lapply(
            1:p,
            function(k) {Bopt[, d * k + 1:d]})
    } else moderator_coeffs = NULL
    if (multi) {
        idiographic_coeffs <- lapply(
            1:n,
            function(k) {Bopt[, d *(p + k) + 1:d]})
    } else idiographic_coeffs = NULL

    ############################
    ## Build individual graphs:
    B_total <- list()
    for (ind in 1:n) {

        B_total[[ind]] <- Bopt[, 1:d]

        if (p > 0) {## Incorporating moderated effects.
            x <- X[ind, ]
            for (ent_ind in 1:p) {
                B_total[[ind]] <- B_total[[ind]] + x[ent_ind] * Bopt[, ent_ind * d + (1:d)]
            }
        } 
        if (multi) {## Incorporating fully idiographic effects.
            B_total[[ind]] <- B_total[[ind]] + Bopt[, (p + ind) * d + (1:d)]
        }
    }

    OUTPUT <- list(
        lambda1 = lambdas1,
        ratios = ratios,
        eval.type = paste0("cv.", cv.type),
        eval.mat = CV.MSFE,
        lambda1_opt = lambda1_opt,
        ratio_opt = ratio_opt,
        lambda2_opt = lambda1_opt * ratio_opt,
        opt_coeffs = Bopt,
        joint_coeffs = joint_coeffs,
        moderator_coeffs = moderator_coeffs,
        idiographic_coeffs = idiographic_coeffs,
        bysubject_coeffs = B_total)

    return(OUTPUT)
}
