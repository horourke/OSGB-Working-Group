
wd <- getwd()
req_lib_dir <- paste0(wd,"/req_lib")
packageVersion("rlang")


###################### Creating folders:
subfolder_new        <- paste0("req_lib/")

if (!dir.exists(subfolder_new)) {
       dir.create(subfolder_new)
}


## Technical requirements:
if(!require(e1071, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("e1071", repos = "https://archive.linux.duke.edu/cran/")
  library(e1071)
}
if(!require(Rcpp, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("Rcpp", repos = "https://archive.linux.duke.edu/cran/")
  library(Rcpp)
}
if(!require(RcppArmadillo, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("RcppArmadillo", repos = "https://archive.linux.duke.edu/cran/")
  library(RcppArmadillo)
}
if(!require(magrittr, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("magrittr", repos = "https://archive.linux.duke.edu/cran/")
  library(magrittr)
}
if(!require(lubridate, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("lubridate", repos = "https://archive.linux.duke.edu/cran/")
  library(lubridate)
}
if(!require(tidyverse, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("tidyverse", repos = "https://archive.linux.duke.edu/cran/")
  library(tidyverse)
}
if(!require(lobstr, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("lobstr", repos = "https://archive.linux.duke.edu/cran/")
  library(lobstr)
}
if(!require(remotes, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("remotes", repos = "https://archive.linux.duke.edu/cran/")
  library(remotes)
}



## Model fitting packages:
if(!require(glmnet, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("glmnet", repos = "https://archive.linux.duke.edu/cran/")
  library(glmnet)
}
if(!require(mvtnorm, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("mvtnorm", repos = "https://archive.linux.duke.edu/cran/")
  library(mvtnorm)
}
#if(!require(gimme, lib = req_lib_dir)){
#  .libPaths(req_lib_dir)
#  install.packages("gimme", repos = "https://archive.linux.duke.edu/cran/")
#  library(gimme)
#}
if(!require(multivar, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  remotes::install_github("zackfisher/multivar", lib = req_lib_dir)
  library(multivar)
}
if(!require(BigVAR, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("BigVAR", repos = "https://archive.linux.duke.edu/cran/")
  library(BigVAR)
}

## Plotting packages:
if(!require(ggh4x, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("ggh4x", repos = "https://archive.linux.duke.edu/cran/")
  library(ggh4x)
}
if(!require(gridExtra, lib = req_lib_dir)){
  .libPaths(req_lib_dir)
  install.packages("gridExtra", repos = "https://archive.linux.duke.edu/cran/")
  library(gridExtra)
}
if (!require(plot.matrix, lib = req_lib_dir)) {
  .libPaths(req_lib_dir)
  install.packages("plot.matrix", repos = "https://archive.linux.duke.edu/cran/")
  library(plot.matrix)
}



## Other packages: necessary?
#if(!require(matrixcalc, lib = req_lib_dir)){
#  .libPaths(req_lib_dir)
#  install.packages("matrixcalc", repos = "https://archive.linux.duke.edu/cran/")
#  library(matrixcalc)
#}
#if(!require(pracma, lib = req_lib_dir)){
#  .libPaths(req_lib_dir)
#  install.packages("pracma", repos = "https://archive.linux.duke.edu/cran/")
#  library(pracma)
#}
#if(!require(MASS, lib = req_lib_dir)){
#  .libPaths(req_lib_dir)
#  install.packages("MASS", repos = "https://archive.linux.duke.edu/cran/")
#  library(MASS)
#}
#if(!require(Matrix, lib = req_lib_dir)){
#  .libPaths(req_lib_dir)
#  install.packages("Matrix", repos = "https://archive.linux.duke.edu/cran/")
#  library(Matrix)
#}
#if(!require(RSpectra, lib = req_lib_dir)){
#  .libPaths(req_lib_dir)
#  install.packages("RSpectra", repos = "https://archive.linux.duke.edu/cran/")
#  library(RSpectra)
#}
#if(!require(spam, lib = req_lib_dir)){
#  .libPaths(req_lib_dir)
#  install.packages("spam", repos = "https://archive.linux.duke.edu/cran/")
#  library(spam)
#}
#if(!require(readr, lib = req_lib_dir)){
#  .libPaths(req_lib_dir)
#  install.packages("readr", repos = "https://archive.linux.duke.edu/cran/")
#  library(readr)
#}

## bla bla bla