## Pending:
## 1: Finish testing if ada.modvar works properly and without typos.
## 2: Implement ada.modvar in simulations.
## 3: run and analyze results.


#####################################################
#####################################################
#####################################################
#####################################################
## Example: estimating fully idiographic structure
##          without moderators.

library(Rcpp)
library(RcppArmadillo)
library(magrittr)
library(plot.matrix)
Rcpp::sourceCpp("modvar/matrix_fista.cpp")
source("modvar/auxfunct.r")
source("modvar/cv.modvar.R")
source("modvar/bic.modvar.r")
source("modvar/adaweights.r")
source("modvar/ada.modvar.r")


set.seed(123)
T  <- 100
p  <- 2
d  <- 5
n  <- 4
q  <- d * (n + 1)

X <- matrix(runif(n * p, -1, 1), n, p)
Btrue <- list()
sigma <- 0.1
Ylist <- list()
for(i in 1:n) {
    # True sparse coefficient matrix
    Btrue[[i]] <- 0.5 * diag(d)
    Btrue[[i]][i, i+1] <- -0.5
    
    print(Btrue[[i]])

    Ymat <- matrix(0, nrow = T + 100, ncol = d)
    for(t_ind in 2:(T + 100)) {
        Ymat[t_ind, ] <- Btrue[[i]] %*% as.vector(Ymat[t_ind - 1, ]) + rnorm(d, sd = sigma)
    }
    Ylist[[i]] <- Ymat[101:(T + 100), ]
}
lapply(Ylist, dim)



lambda1vec <- 10^(seq(0, -4, length.out = 25))
ratiovec   <- 10^(seq(5, -5, length.out = 20))
Bcvprev <- array(0, c(q,d, 4*5))

res <- cv.modvar(
    Ylist, X = NULL, lambda1vec, ratiovec, multi = TRUE,
    cv.type = "rolling")
opt_coeffs <- res$opt_coeffs

adaW <- adaweights(
    modvar_fit = res, 
    p = 0,
    multi = TRUE, 
    alpha = 3, 
    inf = 1e10, 
    thr = 1e-4)

res.ada <- cv.modvar(
    Ylist, X = NULL, lambda1vec, c(1), weights = adaW,
    multi = TRUE, cv.type = "rolling")
opt_coeffs_ada <- res.ada$opt_coeffs

par(mar = c(5.1, 4.1, 4.1, 4.1), mfcol = c(4,1))
plot(log(res.ada$eval.mat, 10))
plot(t(adaW), breaks = 30)
plot(opt_coeffs_ada, breaks = 30)
plot(abs(opt_coeffs_ada) > 10^(-1))
res$lambda1_opt
res$lambda2_opt

par(mar = c(5.1, 4.1, 4.1, 4.1), mfrow = c(1,1))
plot(res.ada$joint_coeffs, breaks = seq(-0.6,0.6, length.out = 30))
res.ada$moderator_coeffs

par(mfrow = c(3,1)) 
for(i in 1:3) {
    plot(res.ada$bysubject_coeffs[[i]], breaks = seq(-1,1, length.out = 30))
} 


#####################################################
#####################################################
#####################################################
#####################################################
## Example: estimating MOD-VAR with moderators, but no
##          fully idiographic component.

library(Rcpp)
library(RcppArmadillo)
library(magrittr)
library(plot.matrix)
setwd("rvar_project/experiments/07_modvar_full_fista/")
Rcpp::sourceCpp("modvar/matrix_fista.cpp")
source("modvar/auxfunct.r")
source("modvar/cv.modvar.R")

set.seed(123)

T  <- 30
p  <- 2
d  <- 5
n  <- 20
q  <- d * (p + 1)

X <- matrix(runif(n * p, -1, 1), n, p)
Btrue <- list()
sigma <- 0.1
Ylist <- list()
for(i in 1:n) {
    # True sparse coefficient matrix
    Btrue[[i]] <- 0.5 * diag(d)
    
    
    Bx1 <- matrix(0,d,d)
    Bx2 <- matrix(0,d,d)
    Bx1[3,1] <- -0.5 * X[i, 1] 
    Bx2[2,5] <- -0.5 * X[i, 2]
    
    Btrue[[i]] <- Btrue[[i]] + Bx1 + Bx2
    print(i)
    print(Btrue[[i]])

    Ymat <- matrix(0, nrow = T + 100, ncol = d)
    for(t_ind in 2:(T + 100)) {
        Ymat[t_ind, ] <- Btrue[[i]] %*% as.vector(Ymat[t_ind - 1, ]) + rnorm(d, sd = sigma)
    }
    Ylist[[i]] <- Ymat[101:(T + 100), ]
    
    print(length(Ylist))
}
lapply(Ylist, dim)



lambda1vec <- 10^(seq(0, -6, length.out = 25))
ratiovec   <- 10^(seq(3, -3, length.out = 20))
Bcvprev <- array(0, c(q,d, 4*5))

res <- cv.modvar(
    Ylist, X = X, lambda1vec, ratiovec, multi = FALSE,
    cv.type = "rolling")
opt_coeffs <- res$opt_coeffs

adaW <- adaweights(
    modvar_fit = res, 
    p = p,
    multi = FALSE, 
    alpha = 3, 
    inf = 1e10, 
    thr = 1e-4)

res.ada <- cv.modvar(
    Ylist, X = X, lambda1vec, c(1), weights = adaW,
    multi = FALSE, cv.type = "rolling")
opt_coeffs_ada <- res.ada$opt_coeffs


par(mar = c(5.1, 4.1, 4.1, 4.1), mfcol = c(3,1))
plot(abs(opt_coeffs) > 10^(-4))
plot(adaW)
plot(abs(opt_coeffs_ada) > 10^(-4))

par(mar = c(5.1, 4.1, 4.1, 4.1), mfcol = c(4,1))
plot(log(res.ada$eval.mat, 10))
plot(t(adaW), breaks = 30)
plot(opt_coeffs_ada, breaks = 30)
plot(abs(opt_coeffs_ada) > 10^(-5))
res$lambda1_opt
res$lambda2_opt


par(mar = c(5.1, 4.1, 4.1, 4.1), mfrow = c(1,1))
plot(res.ada$joint_coeffs, breaks = seq(-0.6,0.6, length.out = 30))
res.ada$moderator_coeffs

par(mfrow = c(3,1)) 
for(i in 1:3) {
    plot(res.ada$bysubject_coeffs[[i]], breaks = seq(-1,1, length.out = 30))
} 







#####################################################
#####################################################
#####################################################
#####################################################
## Example: estimating MOD-VAR with moderators, and
##          fully idiographic components.

library(Rcpp)
library(RcppArmadillo)
library(magrittr)
library(plot.matrix)
setwd("rvar_project/experiments/07_modvar_full_fista/")
Rcpp::sourceCpp("modvar/matrix_fista.cpp")
source("modvar/auxfunct.r")
source("modvar/cv.modvar.R")

set.seed(123)

T  <- 70
p  <- 2
d  <- 10
n  <- 9
q <- d * (n + p + 1)

X <- matrix(runif(n * p, -1, 1), n, p)
Btrue <- list()
sigma <- 0.1
Ylist <- list()
for(i in 1:n) {
    # True sparse coefficient matrix
    Btrue[[i]] <- 0.5 * diag(d)
    Btrue[[i]][i, i+1] <- -0.7
    
    Bx1 <- matrix(0,d,d)
    Bx2 <- matrix(0,d,d)
    Bx1[3,1] <- -0.5 * X[i, 1] 
    Bx2[2,5] <- -0.5 * X[i, 2]
    
    Btrue[[i]] <- Btrue[[i]] + Bx1 + Bx2
    print(i)
    print(Btrue[[i]])

    Ymat <- matrix(0, nrow = T + 100, ncol = d)
    for(t_ind in 2:(T + 100)) {
        Ymat[t_ind, ] <- Btrue[[i]] %*% as.vector(Ymat[t_ind - 1, ]) + rnorm(d, sd = sigma)
    }
    Ylist[[i]] <- Ymat[101:(T + 100), ]
    
    print(length(Ylist))
}
lapply(Ylist, dim)



lambda1vec <- 10^(seq(2, -5, length.out = 10))
ratiovec   <- 10^(seq(3, -3, length.out = 10))
Bcvprev <- array(0, c(q,d, 4*5))

res <- cv.modvar(
    Ylist, X = X, lambda1vec, ratiovec, multi = TRUE,
    cv.type = "rolling")
opt_coeffs <- res$opt_coeffs

adaW <- adaweights(
    modvar_fit = res, 
    p = p,
    multi = TRUE, 
    alpha = 3, 
    inf = 1e10, 
    thr = 1e-4)
plot(log(adaW, base = 10))

res.ada <- cv.modvar(
    Ylist, X = X, lambda1vec, c(1), weights = adaW,
    multi = TRUE, cv.type = "rolling")
opt_coeffs_ada <- res.ada$opt_coeffs

par(mar = c(5.1, 4.1, 4.1, 4.1), mfcol = c(3,1))
plot(abs(opt_coeffs) > 10^(-4))
plot(adaW)
plot(abs(opt_coeffs_ada) > 10^(-4))

par(mar = c(5.1, 4.1, 4.1, 4.1), mfcol = c(4,1))
plot(log(t(adaW), base = 10))
plot(t(adaW), breaks = 30)
plot(opt_coeffs_ada, breaks = 30)
plot(abs(opt_coeffs_ada) > 10^(-5))
res$lambda1_opt
res$lambda2_opt


par(mar = c(5.1, 4.1, 4.1, 4.1), mfrow = c(1,1))
plot(res.ada$joint_coeffs, breaks = seq(-0.6,0.6, length.out = 30))
res.ada$moderator_coeffs

par(mfrow = c(3,1)) 
for(i in 1:3) {
    plot(res.ada$bysubject_coeffs[[i]], breaks = seq(-1,1, length.out = 30))
} 









#####################################################
#####################################################
#####################################################
#####################################################
## Example: estimating MOD-VAR with moderators, but no
##          fully idiographic component. Using by-subject
##          CV.

library(Rcpp)
library(RcppArmadillo)
library(magrittr)
library(plot.matrix)
setwd("rvar_project/experiments/07_modvar_full_fista/")
Rcpp::sourceCpp("modvar/matrix_fista.cpp")
source("modvar/auxfunct.r")
source("modvar/cv.modvar.R")

set.seed(123)

T  <- 30
p  <- 2
d  <- 5
n  <- 20
q  <- d * (p + 1)

X <- matrix(runif(n * p, -1, 1), n, p)
Btrue <- list()
sigma <- 0.1
Ylist <- list()
for(i in 1:n) {
    # True sparse coefficient matrix
    Btrue[[i]] <- 0.5 * diag(d)
    
    
    Bx1 <- matrix(0,d,d)
    Bx2 <- matrix(0,d,d)
    Bx1[3,1] <- -0.5 * X[i, 1] 
    Bx2[2,5] <- -0.5 * X[i, 2]
    
    Btrue[[i]] <- Btrue[[i]] + Bx1 + Bx2
    print(i)
    print(Btrue[[i]])

    Ymat <- matrix(0, nrow = T + 100, ncol = d)
    for(t_ind in 2:(T + 100)) {
        Ymat[t_ind, ] <- Btrue[[i]] %*% as.vector(Ymat[t_ind - 1, ]) + rnorm(d, sd = sigma)
    }
    Ylist[[i]] <- Ymat[101:(T + 100), ]
    
    print(length(Ylist))
}
lapply(Ylist, dim)



lambda1vec <- 10^(seq(0, -6, length.out = 25))
ratiovec   <- 10^(seq(3, -3, length.out = 20))
Bcvprev <- array(0, c(q,d, 4*5))

res <- cv.modvar(
    Ylist, X = X, lambda1vec, ratiovec, multi = FALSE,
    cv.type = "bysubject", nfolds = 5)
opt_coeffs <- res$opt_coeffs


adaW <- adaweights(
    modvar_fit = res, 
    p = p,
    multi = FALSE, 
    alpha = 3, 
    inf = 1e10, 
    thr = 1e-4)

res.ada <- cv.modvar(
    Ylist, X = X, lambda1vec, c(1), weights = adaW,
    multi = FALSE, cv.type = "rolling")
opt_coeffs_ada <- res.ada$opt_coeffs

par(mar = c(5.1, 4.1, 4.1, 4.1), mfcol = c(3,1))
plot(abs(opt_coeffs) > 10^(-4))
plot(adaW)
plot(abs(opt_coeffs_ada) > 10^(-4))

par(mar = c(5.1, 4.1, 4.1, 4.1), mfcol = c(3,1))
plot(abs(opt_coeffs) > 10^(-4))
plot(adaW)
plot(abs(opt_coeffs_ada) > 10^(-4))

par(mar = c(5.1, 4.1, 4.1, 4.1), mfcol = c(4,1))
plot(log(res.ada$eval.mat, 10))
plot(t(adaW), breaks = 30)
plot(opt_coeffs_ada, breaks = 30)
plot(abs(opt_coeffs_ada) > 10^(-5))
res$lambda1_opt
res$lambda2_opt


par(mar = c(5.1, 4.1, 4.1, 4.1), mfrow = c(1,1))
plot(res.ada$joint_coeffs, breaks = seq(-0.6,0.6, length.out = 30))
res.ada$moderator_coeffs

par(mfrow = c(3,1)) 
for(i in 1:3) {
    plot(res.ada$bysubject_coeffs[[i]], breaks = seq(-1,1, length.out = 30))
} 

