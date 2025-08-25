
library(magrittr)
library(glmnet)
library(plot.matrix)
library(tidyverse)


#################################################################
#################################################################
#################################################################
#################################################################
## In the implementation of the R-VAR, I noticed that the estimators
## are row-wise sparse, but not overall sparse. I am working on this
## script to explore if there is any issue with my calculations, or
## if there are underlying assumptions in the "mgaussian" mode of
## glmnet.


set.seed(100)
n <- 100
p <- 20
q <- 3

bmin  <- 5
bmax  <- 10
prob <- 0.1


X <- matrix(rnorm(n * p), ncol = p)
signs <- c(-1, 1, 0)
Beta <- matrix(runif(p * q, bmin, bmax) *
                 sample(signs, 
                        size = p * q, 
                        replace = TRUE, 
                        prob = c(prob/2, prob/2, 1 - prob)), 
               ncol = q)
E <- matrix(rnorm(n * q, sd =  0.2), ncol = q)

Y <- X %*% Beta + E

plot(Beta, border = NA, breaks = 20)



glmnet_fit <- glmnet(x = X, y = Y, family = "mgaussian")

str(glmnet_fit$beta)
par(mfrow = c(1,4))
plot(as.matrix(Beta), border = NA, breaks = 100)
plot(as.matrix(glmnet_fit$beta$y1), border = NA, breaks = 100)
plot(as.matrix(glmnet_fit$beta$y2), border = NA, breaks = 100)
plot(as.matrix(glmnet_fit$beta$y3), border = NA, breaks = 100)


par(mfrow = c(1,4))
plot(as.matrix(Beta) != 0, border = NA)
plot(as.matrix(glmnet_fit$beta$y1) != 0, border = NA)
plot(as.matrix(glmnet_fit$beta$y2) != 0, border = NA)
plot(as.matrix(glmnet_fit$beta$y3) != 0, border = NA)


## The answer is yes: 
## When glmnet uses the "mgaussian family, it ensures that the solutions
## are "row-wise sparse", meaning that if an entry is non-zero for one
## of the regression problems, it will be non-zero for all of them. This
## will be a problem since we do not rely on this principle. 

## Does this mean that we need to code our own LASSO function?

#################################################################
#################################################################
#################################################################
#################################################################

## Should we try to solve the problem with the ADMM algorithm?
## With the LASSO coordinate descent? 
## With a different tool?

## Lets explore if we can do it with the RMTL package:

install.packages("RMTL")
library(RMTL)

#create simulated data
datar <- Create_simulated_data(Regularization="L21", type="Regression")
str(datar$X)
str(datar$Y)

#perform the cross validation
cvfitr <- cvMTL(datar$X, datar$Y, type="Regression", Regularization="L21", Lam1_seq=10^seq(1,-4, -1),  Lam2=0, opts=list(init=0,  tol=10^-6, maxIter=1500), nfolds=5, stratify=FALSE, parallel=FALSE)

# meta-information and results of CV 
#sequence of lam1
cvfitr$Lam1_seq
#> [1] 1e+01 1e+00 1e-01 1e-02 1e-03 1e-04

#value of lam2
cvfitr$Lam2
#> [1] 0

#the output lam1 value with minimum CV error
print (paste0("estimated lam1: ", cvfitr$Lam1.min))
#> [1] "estimated lam1: 0.1"

#plot CV errors across lam1 sequence in the log space
plot(cvfitr)


#################################################################
#################################################################
#################################################################
#################################################################

X1 <- matrix(rnorm(n * p), ncol = p)
X2 <- matrix(rnorm(n * p), ncol = p)
X3 <- matrix(rnorm(n * p), ncol = p)

X_list <- list(X, X, X)
Y_list <- list(matrix(Y[,1], ncol = 1),
               matrix(Y[,2], ncol = 1),
               matrix(Y[,3], ncol = 1))
str(X_list)
str(Y_list)
cvfitr <- cvMTL(X = X_list, Y = Y_list, type = "Regression", Regularization = "Lasso", 
                Lam1_seq=10^seq(1,-4, -1),  Lam2=0, opts=list(init=0,  tol=10^-6, maxIter=1500), 
                nfolds=10, parallel=FALSE)

# meta-information and results of CV 
#sequence of lam1
cvfitr$Lam1_seq
#> [1] 1e+01 1e+00 1e-01 1e-02 1e-03 1e-04

#value of lam2
cvfitr$Lam2
#> [1] 0

#the output lam1 value with minimum CV error
print (paste0("estimated lam1: ", cvfitr$Lam1.min))
#> [1] "estimated lam1: 0.1"

#plot CV errors across lam1 sequence in the log space
plot(cvfitr)



## After going through some exploration, I have realized that the 
## package RMTL will be really useful for solving our problem. 
##
## The main issue: it does not allow for weighted penalties, which
##                  are essential to our problem.
##
## Solution: extract the vital code we need from the package, and 
##            simply edit it to allow for a weighted version of the
##            problem. 
##
## Reflections: will this be easy? I DONT KNOW! WORTH THE TRY.


#################################################################
#################################################################
#################################################################
#################################################################
## 
