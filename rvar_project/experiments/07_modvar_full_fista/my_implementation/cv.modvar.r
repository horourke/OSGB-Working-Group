
###########################################################
###########################################################
## build_design:
##  Creates the matrices $W$ and $Y$ to solve the optimization 
##  problem of MOD-VAR, of the form,
##    $$ \hat{B} = \argmin_{B} \frac{1}{2N} \|Y - WB\|_F^2 + \lambda_1\|B\|_{1,w}.$$
##
##  INPUTS:
##      Ylist   : list of length n, containing matrices of dimensions 
##                  $Tk \times d$, with the time series data for 
##                  each subject.
##      X       : numeric matrix of dimensions $n x p$. The $k$-th 
##                  row X[k,] must contain the moderator variables of
##                  subject $k$, associated with the $k$-th matrix
##                  in Ylist[[k]].
##      multi   : if multi = TRUE, we add a fully idiographic
##                  component to each subject, independent from
##                  the moderator effects. 
##
##  OUTPUTS:
##      Wmat    : Design matrix of dimension $N\times q$, where 
##                  $N = \sum_k (T_k-1)$ and $q = d (p+1)$ if 
##                  multi = FALSE, or $q = d (n+p+1)$ if multi = TRUE.
##      Ymat    : Response matrix of dimension $N \times d$.
## 
build_design <- function(Ylist, X = NULL, multi = FALSE) {

    n <- length(Ylist)
    d <- ncol(Ylist[[1]])
    p <- ncol(X)

    ## Response matrix:
    Yl <- Ylist %>%
        lapply(function(x) x[-1,])
    ## Pre-processed predictor matrix:
    Zl <- Ylist %>%
        lapply(function(x) {T = nrow(x); return(x[-T,])} )
    Wl <- list()

    ## Building MOD-VAR predictor matrix W = [W,   X x W,   I x W]
    
    for (subj_ind1 in 1:n) {
        Zfirst  <- Zl[[subj_ind1]]
        Mzero   <- matrix(0, nrow = nrow(Zfirst), ncol = d)
        Wl[[subj_ind1]] <- Zfirst

        if (!is.null(X)) {
            for (mod_ind in 1:p) {
                Zmod <- X[subj_ind1 , mod_ind] * Zfirst
                Wl[[subj_ind1]] <- cbind(Wl[[subj_ind1]], Zmod)
            }
        }
        if (multi) {
            for (subj_ind2 in 1:n) {
        
                if(subj_ind1 != subj_ind2) {
                    Wl[[subj_ind1]] <- cbind(Wl[[subj_ind1]], Mzero)
                } else {
                    Wl[[subj_ind1]] <- cbind(Wl[[subj_ind1]], Zfirst)
                }
            }
        }
    }
    
    Wmat <- do.call(rbind, Wl)
    Ymat <- do.call(rbind, Yl)

    OUTPUT <- list(Wmat = Wmat, Ymat = Ymat)
    return(OUTPUT)
}







###########################################################
###########################################################
## build_rolling_windows:
##  Creates the indexes necessary to perform rolling cross
##  validation for multi-subject data based on Ymat and Wmat.
##
##  INPUTS:
##      Ts          : numeric vector containing the number of time measurements
##                      per subject.
##      folds_max   : maximum number of windows to use.
##                      aimed at improving computational performance.
##      
##  OUTPUTS:
##      cv.train_list   : List. Contains the rows of Wmat and Ymat
##                          used for training at each rolling window 
##                          CV-fold.
##      cv.test_list    : List. Contains the rows of Wmat and Ymat
##                          used for testing the CV-fold. 
##      nfolds          : Numeric counting the number of 
##                          rolling windows.
##
build_rolling_windows <- function(Ts, folds_max = 20) {
    n <- length(Ts)
    Tms <- Ts - 1
    trange  <- min(c(floor(Ts / 3), folds_max)) 
    T1s     <- Tms - trange - 3
    starts  <- cumsum(c(1, Tms))[-(n+1)]

    cv.train_list <- list()
    cv.test_list <- list()
    for (t1 in 0:trange) {

        ## Build rolling window indexes:
        cv.train_temp <- list()
        cv.test_temp <- rep(0, n)
        for (i in 1:n) {
            cv.train_temp[[i]] <- starts[i]:(starts[i] + T1s[i] + t1)
            cv.test_temp[i] <- (starts[i] + T1s[i] + t1 + 1)
        }
        cv.train_list[[t1+1]] <- unlist(cv.train_temp)
        cv.test_list[[t1+1]] <- cv.test_temp
    }

    OUTPUT <- list(
        cv.train_list = cv.train_list,
        cv.test_list = cv.test_list,
        nfolds = length(cv.train_list))
    return(OUTPUT)

}







###########################################################
###########################################################
## build_bysubject_windows:
##  Creates the indexes necessary to perform by-subject cross
##  validation for multi-subject data based on Ymat and Wmat.
##  This procedure divides the folds by subjects, and thus
##  is only applicable when multi = FALSE.
##
##  INPUTS:
##      Ts      : numeric vector containing the number of time measurements
##                  per subject.  
##      nfolds  : number of folds. 
##      
##  OUTPUTS:
##      cv.train_list   : List. Contains the rows of Wmat and Ymat
##                          used for training at each by-subject
##                          CV-fold.
##      cv.test_list    : List. Contains the rows of Wmat and Ymat
##                          used for testing the CV-fold. 
##      nfolds          : Numeric counting the number of 
##                          rolling windows.
##
build_bysubject_windows <- function(Ts, nfolds) {
    n <- length(Ts)
    Tms <- Ts - 1
    starts  <- cumsum(c(1, Tms))[-(n+1)]
    ends  <- cumsum(c(Tms))
    subj_list <- mapply(
        function(a,b) return(a:b),
        a = starts, b = ends)
    
    nlen_temp   <- ceiling((n / nfolds))
    reordering  <- sample(1:(nfolds*nlen_temp))
    order_temp <- rep(1:nfolds, nlen_temp)[reordering]
    order_folds <- order_temp[1:n]
    
    cv.train_list <- lapply(1:nfolds,
        function(k) return(as.vector(subj_list[,order_folds != k])))
    cv.test_list <- lapply(1:nfolds,
        function(k) return(as.vector(subj_list[,order_folds == k])))
    
    OUTPUT <- list(
        cv.train_list = cv.train_list,
        cv.test_list = cv.test_list,
        nfolds = nfolds)
    return(OUTPUT)

}







###########################################################
###########################################################
## cv_rolling.mod_var:
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
##      Bopt        : optimal estimation of the model B.
##      lambda1_opt : optimal value of lambda1 according to the cv procedure.
##      lambda2_opt : optimal value of lambda2 according to the cv procedure.
##      CV.MSFE     : matrix of cross validation forecasting error. rows correspond
##                      to lambdas1 range and columns to ratios)
##
cv_rolling.mod_var <- function(
    Ylist, X, 
    lambdas1, 
    ratios, 
    multi = FALSE,
    cv.type = c("rolling", "bysubject"), ## "bysubject" only valid if "multi = FALSE"
    nfolds = NULL                        ## Only needed if multi = FALSE, and cv.type = "bysubject"
    ) { 

    ## Build design matrix:
    BD <- build_design(Ylist, X, multi)
    Wmat <- BD$Wmat
    Ymat <- BD$Ymat

    n <- length(Ylist)
    q <- ncol(Wmat)
    N <- nrow(Ymat)
    d <- ncol(Ymat)

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
        
        ## Fit model:
        Wmat.cv <- Wmat[cv.train, ]
        Ymat.cv <- Ymat[cv.train, ]

        Bout  <- weighted_lasso_path(
            X = Wmat.cv,
            Y = Ymat.cv,
            lambda1vec = lambdas1,
            ratiovec = ratios,
            max_iter = 2000,
            Bcvprev = Bcvprev,
            tol = 1e-8)

        ## Forecasting error:
        Ymat.test <- Ymat[cv.test, ]
        Wmat.test <- Wmat[cv.test, ]

        MSFE[, , ind] <- forecast_mse_batch(
            Xf = Wmat.test,
            Yf = Ymat.test,
            Bflat = Bout,
            nlambda1 = nlam1,
            nratios = nrat)
        #print(MSFE[, , ind])
 
    }
    #print(dim(MSFE))
    CV.MSFE <- apply(MSFE, c(1,2), mean)
    #print(dim(CV.MSFE))

    ## Find minimal lambda1, lambda2:
    opt.inds <- which(CV.MSFE == min(CV.MSFE), arr.ind = TRUE)
    print(opt.inds)

    lambda1_opt <- lambdas1[opt.inds[1,1]]
    ratio_opt <- ratios[opt.inds[1,2]]

    ## Fit optimal model:
    Bopt <- weighted_lasso_path(
        X = Wmat,
        Y = Ymat,
        lambda1vec = lambda1_opt,
        ratiovec = ratio_opt,
        max_iter = 2000,
        Bcvprev = array(0, c(q, d, nlam1 * nrat)),
        tol = 1e-8)

    OUTPUT <- list(
        Bopt = Bopt,
        lambda1_opt = lambda1_opt,
        lambda2_opt = lambda1_opt * ratio_opt,
        CV.MSFE = CV.MSFE)
    
    return(OUTPUT)
}
