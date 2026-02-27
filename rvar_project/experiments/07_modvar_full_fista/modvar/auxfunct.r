
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
## build_bic_windows:
##  Creates the indexes necessary to perform by-subject cross
##  validation for multi-subject data based on Ymat and Wmat.
##  This procedure divides the folds by subjects, and thus
##  is only applicable when multi = FALSE.
##
##  INPUTS:
##      Ts      : numeric vector containing the number of time measurements
##                  per subject.  
##      alpha   : Balance between in-sample and out-of-sample data. 
##      
##  OUTPUTS:
##      bic.train_list   : List. Contains the rows of Wmat and Ymat
##                          used for training at each by-subject
##                          CV-fold.
##      bic.test_list    : List. Contains the rows of Wmat and Ymat
##                          used for testing the CV-fold. 
##      nfolds           : Numeric counting the number of 
##                          rolling windows.
##
build_bic_windows <- function(Ts, alpha) {

    n <- length(Ts)
    Tms <- Ts - 1
    starts  <- cumsum(c(1, Tms))[-(n+1)]
    ends    <- cumsum(Tms)
    train.lengths <- floor(Tms * alpha)
    
    bic.train_list <- list()
    bic.test_list <- list()
    for (ind in 1:n) {
        if (train.lengths[ind] + 1 > ends[ind]) {
            stop("Error: no bic.test set. The value of alpha must be strictly between 0 and 1.")
        }
        bic.train_list[[ind]] <- starts[ind]:(starts[ind] + train.lengths[ind])
        bic.test_list[[ind]] <- (starts[ind] + train.lengths[ind] + 1):ends[ind]
    }
    bic.train   <- unlist(bic.train_list)
    bic.test    <- unlist(bic.test_list)

    OUTPUT <- list(
        bic.train = bic.train,
        bic.test = bic.test)
    return(OUTPUT)
}
