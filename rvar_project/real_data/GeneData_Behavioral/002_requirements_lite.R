
wd <- getwd()
req_lib_dir <- paste0(wd,"/req_lib")
packageVersion("rlang")


###################### Creating folders:
subfolder_new        <- paste0("req_lib/")

## Technical requirements:
.libPaths(req_lib_dir)
library(e1071)
library(Rcpp)
library(RcppArmadillo)
library(magrittr)
library(lubridate)
library(tidyverse)
library(lobstr)
library(remotes)
library(expm)
library(readxl)

## Model fitting packages:
library(glmnet)
library(mvtnorm)
# library(gimme)
library(multivar)
library(BigVAR)

## Plotting packages:
library(ggh4x)
library(gridExtra)
library(plot.matrix)







mq_na <- function(x, date, lag = 24, adj = 0, plot = TRUE){
  if (!is.matrix(x))
    x <- as.matrix(x)
  
  if (length(date) != nrow(x))
    stop("date must have one entry per row of x.")
  
  ## convert dates to integers
  if (inherits(date, "Date"))
    date <- as.integer(date)
  
  nc <- ncol(x)
  
  ## covariance at lag 0
  g0 <- var(x)
  ginv <- solve(g0)
  
  ## lookup table
  idx <- seq_along(date)
  names(idx) <- as.character(date)
  
  qm <- 0
  QM <- NULL
  df <- 0
  
  for(h in 1:lag){
    
    ## rows whose successor h days later exists
    id2 <- idx[match(as.character(date + h), names(idx))]
    keep <- !is.na(id2)
    
    if(sum(keep) <= 1){
      
      QM <- rbind(QM,
                  c(h, NA, df + nc^2 - adj, NA))
      next
    }
    
    x1 <- x[id2[keep], , drop=FALSE]
    x2 <- x[keep, , drop=FALSE]
    
    nh <- nrow(x1)
    
    ## lag-h covariance
    g <- cov(x1, x2)
    
    ## finite-sample correction
    g <- g * (nh-1)/nh
    
    hmat <- t(g) %*% ginv %*% g %*% ginv
    
    ## Ljung-Box contribution
    qm <- qm + nh^2 * sum(diag(hmat))/(nh-1)
    
    df <- df + nc^2
    dff <- df - adj
    
    pv <- 1
    if(dff > nc^2-1)
      pv <- 1 - pchisq(qm, dff)
    
    QM <- rbind(QM,
                c(h, qm, dff, pv))
  }
  
  colnames(QM) <- c("m","Q(m)","df","p-value")
  
  cat("Modified Ljung-Box Statistics (irregular dates)\n")
  printCoefmat(QM,digits=3)
  
  if(plot){
    plot(QM[,4],
         ylim=c(0,1),
         xlab="Lag",
         ylab="p-value",
         main="Modified Ljung-Box test")
    abline(h=.05,lty=2,col="blue")
  }
  
  invisible(QM)
}
