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
##      var_fits    : optional list of subject VAR fits for multi = TRUE, p = 0.
##      p           : number of moderator variables.
##      multi       : logical indicating whether fully idiographic components are present.
##
##  OUTPUTS:
##      weights     : matrix of adaptive weights.
##
adaweights <- function(
    modvar_fit, 
    var_fits = NULL,
    p = 0,
    multi = FALSE, 
    alpha, 
    inf = NULL, 
    thr = 1e-5) {

    entry_weight <- function(x, alpha = 1, thr = 1e-5, inf = 1e10) {
        ifelse(abs(x) <= thr, inf, 1 / abs(x)^alpha)
    }

    if (is.null(modvar_fit)) {
        if (is.null(var_fits)) {
            stop("modvar_fit or var_fits must be provided to compute adaptive weights.")
        }
        if (!(p == 0 && multi)) {
            stop("var_fits-only adaptive weights are only supported when p == 0 and multi == TRUE.")
        }

        B_est_list <- lapply(var_fits, coef) %>% lapply(as.matrix)
        n <- length(B_est_list)
        d <- ncol(B_est_list[[1]])
        B_array <- array(unlist(B_est_list), dim = c(d, d, n))

        beta0_hat <- apply(B_array, c(1,2), median)
        min_nz <- abs(beta0_hat)[abs(beta0_hat) > thr] %>% min()
        if (is.null(inf)) {
            inf_val <- (min_nz/2)^(-alpha)
        } else {
            inf_val <- inf
        }

        W_common <- matrix(
            entry_weight(beta0_hat, alpha = alpha, thr = thr, inf = inf_val),
            nrow = d, ncol = d)

        W_indiv <- lapply(
            seq_len(n),
            function(i) {
                mat <- beta0_hat - B_array[,,i]
                matrix(
                    entry_weight(mat, alpha = alpha, thr = thr, inf = inf_val),
                    nrow = d, ncol = d)
            })

        M1 <- cbind(W_common, do.call(cbind, W_indiv))
        return(t(M1))
    }

    if (max(abs(modvar_fit$opt_coeffs)) < thr) {
        mat <- modvar_fit$opt_coeffs
        mat <- mat + 1
        return(t(mat))
    }

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
            entry_weight(x, alpha = alpha, thr = thr, inf = inf_val)
        })

    if (multi && !is.null(var_fits)) {
        B_est_list <- lapply(var_fits, coef) %>% lapply(as.matrix)
        B_array <- array(unlist(B_est_list), dim = c(d, d, n))
        beta0_hat <- apply(B_array, c(1,2), median)

        W_common <- matrix(
            entry_weight(beta0_hat, alpha = alpha, thr = thr, inf = inf_val),
            nrow = d, ncol = d)

        W_indiv <- lapply(
            seq_len(n),
            function(i) {
                mat <- beta0_hat - B_array[,,i]
                matrix(
                    entry_weight(mat, alpha = alpha, thr = thr, inf = inf_val),
                    nrow = d, ncol = d)
            })
        W_ind <- do.call(cbind, W_indiv)

        M1[, 1:d] <- W_common
        M1 <- cbind(M1, W_ind)
    } else if (multi) {
        beta0_hat <- modvar_fit$bysubject_coeffs %>%
            unlist() %>% array(c(d, d, n)) %>%
            apply(c(1,2), FUN = stats::median)

        W_common <- matrix(
            entry_weight(beta0_hat, alpha = alpha, thr = thr, inf = inf_val),
            nrow = d, ncol = d)

        W_indiv <- lapply(
            seq_len(n),
            function(i) {
                mat <- beta0_hat - modvar_fit$bysubject_coeffs[[i]]
                matrix(
                    entry_weight(mat, alpha = alpha, thr = thr, inf = inf_val),
                    nrow = d, ncol = d)
            })
        W_ind <- do.call(cbind, W_indiv)

        M1[, 1:d] <- W_common
        M1 <- cbind(M1, W_ind)
    }

    output <- t(M1)
    return(output)
}
