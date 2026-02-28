###########################################################
###########################################################
source("modvar/auxfunct.r")
###########################################################
###########################################################


###########################################################
###########################################################
## bic.modvar:
##  Solver that applies BIC for the selection
##  of tuning parameters for our proposed MOD-VAR.
##  It estimates the CV-MSFE and selects the combination of
##  tuning parameters that minimizes the error.
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
##      alpha       : Balance between in-sample and out-of-sample data. 
##      
##  OUTPUTS:
##      Bopt        : optimal estimation of the model B.
##      lambda1_opt : optimal value of lambda1 according to the cv procedure.
##      lambda2_opt : optimal value of lambda2 according to the cv procedure.
##      BIC.MSFE    : matrix of cross validation forecasting error. rows correspond
##                      to lambdas1 range and columns to ratios)
##
bic.modvar <- function(
    Ylist, X, 
    lambdas1, 
    ratios, 
    weights = NULL,
    multi = FALSE,
    alpha = 0.90) { 

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
    ## Building the training indexes:
    nlam1 <- length(lambdas1)
    nrat  <- length(ratios)
    Bcvprev <- array(0, c(q, d, nlam1 * nrat))

    Ts  <- unlist(lapply(Ylist, nrow))
    bic.windows <- build_bic_windows(Ts, alpha)
    bic.train    <- bic.windows$bic.train
    
    ############################
    ## Fit model:
    Bcvprev <- array(0, c(q, d, nlam1 * nrat))
    n.train <- length(bic.train)
    Wmat.train  <- Wmat[bic.train, ]
    Ymat.train  <- Ymat[bic.train, ]

    Bout  <- weighted_lasso_path(
        X = Wmat.train,
        Y = Ymat.train,
        lambda1vec = lambdas1,
        ratiovec = ratios,
        weights = weights,
        max_iter = 2000,
        Bcvprev = Bcvprev,
        tol = 1e-8)

    ############################
    ## Calculate BIC:
    MSFE_train <- forecast_mse_batch(
        Xf = Wmat.train,
        Yf = Ymat.train,
        Bflat = Bout,
        nlambda1 = nlam1,
        nratios = nrat)
    bic_term2 <- bic_batch(
        Bflat = Bout,
        nf = n.train,
        nlambda1 = nlam1,
        nratios = nrat,
        eps = 1e-5)
    BIC.train <- n.train * log(MSFE_train) + bic_term2

    ############################
    ## Find minimal lambda1, lambda2:
    opt.inds <- which(BIC.train == min(BIC.train), arr.ind = TRUE)
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
        BIC = BIC.train,
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
