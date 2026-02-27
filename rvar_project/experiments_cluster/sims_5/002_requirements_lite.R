
wd <- getwd()
req_lib_dir <- paste0(wd,"/req_lib")
packageVersion("rlang")


###################### Creating folders:
subfolder_new        <- paste0("req_lib/")

## Technical requirements:
.libPaths(req_lib_dir)
library(e1071)
library(magrittr)
library(lubridate)
library(tidyverse)
library(lobstr)

## Model fitting packages:
library(glmnet)
library(mvtnorm)
library(gimme)
library(multivar)
library(BigVAR)

## Plotting packages:
library(ggh4x)
library(gridExtra)
library(plot.matrix)
