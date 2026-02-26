library(Rcpp)
library(RcppArmadillo)
library(magrittr)
library(plot.matrix)
Rcpp::sourceCpp("my_implementation/matrix_fista.cpp")


set.seed(123)

N  <- 100
q  <- 7
d  <- 3

# Design matrix
X <- matrix(rnorm(N * q), N, q)

# True sparse coefficient matrix
Btrue <- matrix(0, q, d)
Btrue[1:d, ] <- 5 * diag(3)
Btrue[d + (1:d), ] <- -5 * diag(3)
Btrue

# Response
Y <- X %*% Btrue + matrix(rnorm(N * d, sd = 0.5), N, d)


lambda1vec <- seq(0.5, 0.1, length.out = 5)
ratiovec   <- seq(5, 0.2, length.out = 4)
Bcvprev <- array(0, c(7,d, 4*5))

Barray <- weighted_lasso_path(
  X = X,
  Y = Y,
  lambda1vec = lambda1vec,
  ratiovec = ratiovec,
  max_iter = 2000,
  Bcvprev = Bcvprev,
  tol = 1e-8
)

dim(Barray)
dim(Bcvprev)


dim(Barray) <- c(q, d, length(lambda1vec), length(ratiovec))
dim(Barray)

Barray[,,1,1]


objective <- function(B, X, Y, W, N) {
  loss <- sum((Y - X %*% B)^2) / (2 * N)
  pen  <- sum(W * abs(B))
  loss + pen
}


dim(Barray) <- c(q, d, length(lambda1vec), length(ratiovec))
dim(Barray)



i <- 4
j <- 4

lambda1 <- lambda1vec[i]
lambda2 <- lambda1 * ratiovec[j]

W <- matrix(lambda2, q, d)
W[1:d, ] <- lambda1
Btest <- Barray[, , i, j]
objective(Btest, X, Y, W, N)
image(abs(Barray[, , i, j]) > 1e-6)


i <- 1
j <- 1
lambda1 <- lambda1vec[i]
lambda2 <- lambda1 * ratiovec[j]
W <- matrix(lambda2, q, d)
W[1:d, ] <- lambda1
Btest <- Barray[, , i, j]
objective(Btest, X, Y, W, N)
image(abs(Barray[, , i, j]) > 1e-6)





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
Rcpp::sourceCpp("my_implementation/matrix_fista.cpp")
source("my_implementation/cv.modvar.R")



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

res <- cv_rolling.mod_var(
    Ylist, X = NULL, lambda1vec, ratiovec, multi = TRUE,
    cv.type = "rolling")

par(mar = c(5.1, 4.1, 4.1, 4.1), mfrow = c(2,1))
plot(log(res$CV.MSFE, 10))

library(plot.matrix)
class(res$Bopt)
dim(res$Bopt)
Bopt <- res$Bopt
dim(Bopt) <- c(25,5)

plot((t(as.matrix(Bopt))), breaks = 30)
res$lambda1_opt
res$lambda2_opt

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

Rcpp::sourceCpp("my_implementation/matrix_fista.cpp")
source("my_implementation/cv.modvar.R")

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

res <- cv_rolling.mod_var(
    Ylist, X = X, lambda1vec, ratiovec, multi = FALSE,
    cv.type = "rolling")

par(mar = c(5.1, 4.1, 4.1, 4.1), mfrow = c(3,1))
plot(res$CV.MSFE, 10)
plot(log(res$CV.MSFE, 10))

library(plot.matrix)
class(res$Bopt)
dim(res$Bopt)
Bopt <- res$Bopt
dim(Bopt) <- c(15,5)

plot((t(as.matrix(Bopt))), breaks = 30)
res$lambda1_opt
res$lambda2_opt






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
Rcpp::sourceCpp("my_implementation/matrix_fista.cpp")
source("my_implementation/cv.modvar.R")

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

res <- cv_rolling.mod_var(
    Ylist, X = X, lambda1vec, ratiovec, multi = TRUE,
    cv.type = "rolling")

par(mar = c(5.1, 4.1, 4.1, 4.1), mfrow = c(3,1))
plot(res$CV.MSFE, 10)
plot(log(res$CV.MSFE, 10))

library(plot.matrix)
class(res$Bopt)
dim(res$Bopt)
Bopt <- res$Bopt
dim(Bopt) <- dim(res$Bopt)[-3]

plot((t(as.matrix(Bopt))), breaks = 30, border = NA)
res$lambda1_opt
res$lambda2_opt
