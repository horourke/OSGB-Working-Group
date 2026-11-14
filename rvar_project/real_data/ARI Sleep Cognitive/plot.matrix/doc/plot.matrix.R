## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.height=4, fig.width=4------------------------------------------------
library('plot.matrix')
# numeric matrix
x <- matrix(runif(35), ncol=5) # create a numeric matrix object
class(x)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(x)
# logical matrix
m <- matrix(runif(35)<0.5, ncol=7)
plot(m)
# text matrix
s <- matrix(sample(letters[1:10], 35, replace=TRUE), ncol=5)
plot(s)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
library('psych')
data <- na.omit(bfi[,1:25])
fa <- fa(data, 5, rotate="varimax")
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(loadings(fa), cex=0.5)

## ----echo=FALSE---------------------------------------------------------------
p   <-formals(plot.matrix:::plot.matrix)
p$na.col <- paste0('"', p$na.col, '"')
out <- matrix('', nrow=length(p), ncol=4)
out[1,1] <- 'plot.matrix('
out[,2] <- names(p)
cp  <- as.character(p)
ct  <- (cp!='')
out[ct,3] <- '=' 
out[ct,4] <- cp[ct] 
out[1:length(p),4] <- cp
out[1,2] <- paste0(out[1,2], ',')
out[2:(length(p)-1),4] <- paste0(out[2:(length(p)-1),4], ',')
out[length(p),2] <- paste0(out[length(p),2], ')')
#out <- cbind(out[,1], apply(out[,2:5], 1, function(e) paste0(e, collapse="")))
knitr::kable(out)

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
# omit all borders
plot(x, border=NA)

## ----fig.height=4, fig.width=4------------------------------------------------
x <- matrix(runif(35), ncol=5) # create a numeric matrix object
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins such that all labels are visible
plot(x, axis.col=list(side=1, las=2), axis.row = list(side=2, las=1))

## ----fig.height=4, fig.width=4------------------------------------------------
param_est <- matrix(runif(25), nrow=5)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins such that all labels are visible
res <- plot(param_est, digits=2, text.cell=list(pos=3, cex=0.75))
sderr_est <- matrix(runif(25)/10, nrow=5)
for (i in 1:nrow(param_est)) {
  for (j in 1:ncol(param_est)) {  
    args <- res$cell.text[[i,j]]
    args$labels <- paste0('(', fmt(sderr_est[i,j], 3), ')')
    args$cex    <- 0.5
    args$pos    <- 1
    do.call(text, args)
  }
}

## ----fig.height=4, fig.width=4------------------------------------------------
param_est <- matrix(runif(25), nrow=5)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins such that all labels are visible
res <- plot(param_est, digits=2, text.cell=list(cex=0.75))
sderr_est <- matrix(runif(25)/10, nrow=5)
for (i in 1:nrow(param_est)) {
  for (j in 1:ncol(param_est)) {  
    args <- res$cell.text[[i,j]]
    args$labels <- paste0('(', round(sderr_est[i,j], 3), ')')
    args$cex    <- 0.6
    args$y      <- args$y-0.3
    do.call(text, args)
  }
}

## ---- warning=FALSE-----------------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins such that all labels are visible
library("plot.matrix")
library("psych")
fa <- fa(iris[,-5], 2)
plot(fa$loadings) 
plot(fa$scores)

## ---- warning=FALSE-----------------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins such that all labels are visible
fa <- fa(iris[,-5], nfactor=2)
plot(fa$loadings) 
plot(as.data.frame(fa$scores))

## ---- warning=FALSE-----------------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins such that all labels are visible
fa <- fa(iris[,-5], nfactors=2)
plot(fa$loadings) 
plot.default(fa$scores)

## ---- warning=FALSE-----------------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins such that all labels are visible
library("plot.matrix")
library("psych")
fa <- fa(iris[,-5], 2)
plot(fa$loadings) 
devtools::unload('plot.matrix') # Package devtools must be installed!
plot(fa$scores)

## ---- include = FALSE---------------------------------------------------------
library('plot.matrix')

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
# we only want the range of x
plot(x, breaks=range(x))     
# we want seven colors
plot(x, breaks=7)  
# user defined breaks, out-of-range entries are colored white
plot(x, breaks=c(0.3,0.5,0.8))  

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
# use a color function
plot(x, col=topo.colors)
# user defined breaks and colors
plot(x, col=c('red', 'green'), breaks=c(0, 0.5, 1))   
# non-numeric matrix
# assign colors
plot(m, col=c('red', 'green'))
# assign colors and breaks directly 
plot(m, col=c('red', 'green'), breaks=c(TRUE, FALSE))

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
x[1,1] <- NA
plot(x, col=topo.colors)
plot(x, col=topo.colors, na.col="red")
plot(s, col=topo.colors, breaks=c('a', 'c', 'e', 'g', 'i'))

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
# delete key
plot(x, key=NULL)

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 5.1, 4.1))   # adapt margins
# move key to the top and make axis text smaller
plot(x, key=list(side=3, cex.axis=0.75), breaks=c(0,1))

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 5.1, 4.1))   # adapt margins
# move key to the top, make axis text smaller, and move key farther away from the plot
plot(x, key=list(side=3, cex.axis=0.75), breaks=c(0,1), spacing.key=2)

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 5.1, 4.1))   # adapt margins
xl  <- array(runif(1e4), c(1e2, 1e2))
brk <- 20
plot(xl, border=NA, breaks=brk, col=heat.colors(brk), 
     key=list(side=4,  font=2, cex.axis=0.75), fmt.key="%.2f", 
     polygon.key=NULL, axis.key=NULL, spacing.key=c(3,2,2))

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
# no plus sign
plot(x, fmt.key="%.3f", breaks=c(0,1))

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
plot(x, xlab="my x label", ylab="my y label")
# The HairEyeColor has its own names
tab <- apply(HairEyeColor,1:2, sum)
plot(tab)

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 5.1, 4.1))   # adapt margins
plot(x, axis.col=list(side=3, cex.axis=0.7), axis.row=list(cex.axis=0.7))
# or alternatively set cex.axis for all axes and use abbreviated positioning
plot(x, axis.col=3, cex.axis=0.7)

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
plot(x, axis.col=NULL, axis.row=NULL, xlab='', ylab='')

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
# change all text output
plot(x, digits=4, text.cell=list(cex=0.5))
# of alternatively use the global parameter cex
plot(x, digits=4, cex=0.5)
# change just matrix entries (no plus sign)
plot(x, fmt.cell='%.2f')

## ----fig.height=4, fig.width=4------------------------------------------------
x <- matrix(c(NA, 1, 2, 3), ncol=2)
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
plot(x, na.col='black', fmt.cell='%.0f')

## ----fig.height=4, fig.width=4------------------------------------------------
x <- matrix(c(NA, 1, 2, 3), ncol=2)
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
plot(x, na.col='black', fmt.cell='%.0f', na.print=FALSE)
plot(x, na.col='black', fmt.cell='%.0f', na.print='Missing')

## ----fig.height=4, fig.width=4------------------------------------------------
x <- matrix(c(NA, 1, 2, 3), ncol=2)
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
plot(x, fmt.cell='%.0f', na.cell=FALSE)

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
# Never replace text color by black or white
plot(x, digits=2, text.cell=list(col="yellow", cex=0.75), max.col=-1) 

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
plot(x, digits=2, text.cell=list(col="yellow", cex=0.75)) 

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
# In fewer cells the text color will be replaced by black
plot(x, digits=2, text.cell=list(col="yellow", cex=0.75), max.col=35) 

## ----fig.height=4, fig.width=4------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 4.1))   # adapt margins
# In more cells the text color will be replaced by black
plot(x, digits=2, text.cell=list(col="yellow", cex=0.75), max.col=140) 

## ----fig.height=4, fig.width=4------------------------------------------------
# numeric matrix
x <- matrix(runif(35), ncol=5) # create a numeric matrix object
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
res <- plot(x)

## -----------------------------------------------------------------------------
#####
names(res)
# parameters of polygon which was used to draw x[3,4]
res$cell.polygon[[3,4]]
#####
# parameters of text which was used for x[3,4]
res$cell.text[[3,4]]  # NULL since no text was drawn
#####
# parameters of polygon which was used to draw the second key element
res$key.polygon[[2]]  

## -----------------------------------------------------------------------------
x <- matrix(runif(35), ncol=5) # create a numeric matrix object
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
res <- plot(x)

for (i in 1:length(res$cell.polygon)) {
  args <- res$cell.polygon[[i]]
  args$col     <- NA       # use no color
  args$density <- 20*x[i]  # density depends on matrix entries
  args$angle   <- 45       
  do.call("polygon", args)
}

## -----------------------------------------------------------------------------
x <- matrix(runif(35), ncol=5) # create a numeric matrix object
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
res <- plot(x)
#
library("png")
# PNGs are created from wikimedia commons images, CC0 by C. Koltzenburg 
# https://commons.wikimedia.org/wiki/File:C.Koltzenburg_-_smiley-yes.xcf
# https://commons.wikimedia.org/wiki/File:C.Koltzenburg_-_smiley_no.xcf
happy <- readPNG(system.file('png', 'happy.png', package="plot.matrix"))
sad   <- readPNG(system.file('png', 'sad.png', package="plot.matrix"))

for (i in 1:length(res$cell.polygon)) {
  args <- res$cell.polygon[[i]]
  if (x[i]>0.5) {
    rasterImage(happy, args$x[1]+0.1, args$y[1]+0.1, args$x[3]-0.1, args$y[2]-0.1)
  } else {
    rasterImage(sad, args$x[1]+0.1, args$y[1]+0.1, args$x[3]-0.1, args$y[2]-0.1)
  }
}

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
library('psych')
data(bfi.2)
fa <- fa(bfi.2, 5, rotate="varimax")
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(loadings(fa), cex=0.5)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
library('psych')
data(bfi.2)
fa <- fa(bfi.2, 5, rotate="varimax")
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(loadings(fa), cex=0.5, gray=TRUE)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
library('psych')
data(bfi.2)
fa <- fa(bfi.2, 5, rotate="varimax")
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(loadings(fa), digits=NA)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
library('psych')
data(bfi.2)
fa <- fa(bfi.2, 5, rotate="varimax")
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(loadings(fa), reorder=FALSE, cex=0.5)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
data(Titanic.cramer)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.assoc(Titanic.cramer))

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
library('datasets')
c <- cor(airquality[,1:4], use="complete")
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.cor(c))

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
data(Titanic.cramer)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.assoc(Titanic.cramer), gray=TRUE)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
data(Titanic.cramer)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.assoc(Titanic.cramer), digits=NA)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
data(Titanic.cramer)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.assoc(Titanic.cramer), reorder=FALSE)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
data(air.pvalue)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.pvalue(air.pvalue))

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
data(air.pvalue)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.pvalue(air.pvalue), gray=TRUE)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
data(air.pvalue)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.pvalue(air.pvalue), digits=NA)

## ----fig.height=5, fig.width=5------------------------------------------------
library('plot.matrix')
data(air.pvalue)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(as.pvalue(air.pvalue), reorder=FALSE)

