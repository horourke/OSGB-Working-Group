install.packages("gimme")
library(gimme)

seed = 11287

# load & view sim data in gimme package
data("simData", package = "gimme") 


# there are some bugs in the source code
#i think the above is correctly coded, but error:
# Error in if (suppressWarnings(tseries::adf.test(time[c, ])$p.value) >  : 
#missing value where TRUE/FALSE needed
# use source code function and work thru:
# https://github.com/GatesLab/gimme/blob/master/R/simulateVAR.R
#see siepe paper for matrices input

# set up matrices for subgroups
a1 <- c(0,0,0,0,0,0.25,
        0,0.44,0,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0
        )

dim(a1) <- c(6,6)
a2 <- c(0,0,0,0,0,0.25,
        0,0.44,0,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0,
        0,0,0,0,0,0
)

dim(a2) <- c(6,6)

alist <- list(a1,a2)

phi1 <- c(0.2,0,0,0,0,0,
          0,0.2,0,0,0,0,
          0,0,0.2,0,0,0,
          0,0,0,0.2,0,0,
          0,0,0,0,0.2,0,
          0,0,0,0,0,0.2
)
dim(phi1) <- c(6,6)
phi2 <- c(0.2,0,0,0,0,0,
          0,0.2,0,0,0,0,
          0,0,0.2,0,0,0,
          0,0,0,0.2,0,0,
          0,0,0,0,0.2,0,
          0,0,0,0,0,0.2
)
dim(phi2) <- c(6,6)
philist <- list(phi1,phi2)

psi1 <- c(0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0
)
dim(psi1) <- c(6,6)
psi2 <- c(0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0,
          0,0,0,0,0,0
)
dim(psi2) <- c(6,6)
psilist <- list(psi1,psi2)

sub <- c(rep(0, 12),rep(1, 12))

# try to sim data with no subs
datsim <- simulateVAR(A = a1,
                      Phi = phi1,
                      Psi = psi1,
                      subAssign = NULL,
                      N = 24,
                      ASign = "random",
                      PhiSign = "random",
                      Obs = 100,
                      indA = 0.01,
                      indPhi = 0.01,
                      indPsi = 0.00)
