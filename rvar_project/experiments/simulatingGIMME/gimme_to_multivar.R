
example <- FALSE 


######################################################################
######################################################################
## Based on a USEM model of the form:
##
##  X(t) = A X(t) + SUM[1-q]{ phi[[q]] X(t-q) } + et
##
## Changes the parameters so the time series can be simulated as a
## VAR model,
##
##  X(t) = SUM[1-q] { Phi[q] X(t-q) } + Et
##
usem_to_var <- function(A, phi_list, q = 1, covD) {
  d         <- ncol(A)
  Inv       <- solve(diag(d) - A)
  
  Phi_list  <- lapply(phi_list, 
                      function(x, Inv) {Inv %*% x},
                      Inv = Inv)
  cov       <- Inv %*% covD %*% t(Inv)
  
  output    <- list(Phi_list = Phi_list, cov = cov)
  
  return(output)
}

usem_to_var1 <- function(A, phi, covD) {
  d         <- ncol(A)
  Inv       <- solve(diag(d) - A)
  
  Phi  <- Inv %*% phi
  cov       <- Inv %*% covD %*% t(Inv)
  
  output    <- list(Phi = Phi, cov = cov)
  
  return(output)
}


######################################################################
######################################################################
## Based on a GIMME model of the form:
##
## Xk(t) = (Ag + Ak) Xk(t) + SUM[1-q] { (phi_g + phi_k) Xk(t-q) } + et
##
## we change the parameters of the time series to be simulated as a
## VAR model
##
## Xk(t) = SUM[1-q] { Phi_k[q] Xk(t-q) } + Ek(t)
##
gimme_to_multivar <- function(Ag, Alist, phi_g, phi_list, q, covD) {
  
  K <- length(Alist)
  VARs <- list()
  for (k in 1:K) {
    Ak <- Ag + Alist[k]
    phik <- mapply(function(phi1, phi2) {return(phi1 + phi2)},
                   list(phi_g),
                   phi_list[[k]])
    
    VAR[[k]] <- usem_to_var(Ak, phik, q, covD)
  }
  
  return(VAR)
}
gimme_to_multivar1 <- function(Ag, Alist, phi_g, phi_list, covD) {
  
  K <- length(Alist)
  VARs <- list()
  for (k in 1:K) {
    Ak <- Ag + Alist[k]
    phik <- phi_g + phi_list[[k]]
    
    VARs[[k]] <- usem_to_var1(Ak, phik, covD)
  }
  
  Phis <- lapply(VARs, function(x) {return(x$Phi)})
  Covs <- lapply(VARs, function(x) {return(x$cov)})
  
  output <- list(K = K, Phis = Phis, Covs = Covs)
  return(output)
}


######################################################################
######################################################################
## Generates a single USEM model
##
generate_usem_pars <- function(d, q, prob_A, prob_phi, 
                               Amin, Amax, phimin, phimax,
                               vmin, vmax) {
  
  sort <- sample(1:d, d, replace = FALSE)
  
  n_entries <- d * (d-1) / 2
  A <- matrix(0, d, d)
  A[upper.tri(A)] <- rbinom(n_entries, 1, prob_A) * 
                        runif(n_entries, Amin, Amax) *
                        sample(c(-1,1), n_entries, TRUE)
  A <- A[sort, sort]
  
  phi_list <- list()
  for (ind in 1:q) {
    phi_temp <- matrix(0, d, d)
    phi_temp[upper.tri(phi_temp)] <- rbinom(n_entries, 1, prob_phi) * 
      runif(n_entries, phimin, phimax) *
      sample(c(-1,1), n_entries, TRUE)
    
    phi_list[[ind]] <- phi_temp[sort, sort]
  }
  
  covD <- diag(runif(d, vmin, vmax))
  
  output <- list(A = A, phi_list = phi_list, q = q, covD = covD)
  return(output)
}



if (example) {
  
  library(plot.matrix)
  par(mar = c(5.1, 4.1, 4.1, 3.5),
      mfrow = c(2,3))
  usem_pars <- generate_usem_pars(d = 10, q = 1, 
                             prob_A = 0.2, prob_phi = 0.2, 
                             Amin = 0.1, Amax = 0.5, 
                             phimin = 0.1, phimax = 0.5,
                             vmin = 0.3, vmax = 0.5)
  plot(usem_pars$A)
  plot(usem_pars$phi_list[[1]])
  plot(usem_pars$covD)
  

  
  var_pars <- usem_to_var(A = usem_pars$A, phi_list = usem_pars$phi_list, 
                           q = usem_pars$q, covD = usem_pars$covD)
  plot(var_pars$Phi_list[[1]])
  plot(var_pars$cov)
  plot(solve(diag(10) - usem_pars$A))

  
  
  usem_sims <- simulateVAR(A = usem_pars$A, ## Contemporaneous effect
                           Phi = usem_pars$phi_list[[1]], ## extemporaneous effects
                           Psi = usem_pars$covD, ## Variances.
                           N = 1,
                           Obs = 500)
  
  sum(usem_sims$A[[1]] != usem_pars$A)
  
  sum(usem_sims$Phi[[1]] != usem_pars$phi_list[[1]])
  
  sum(usem_sims$Psi[[1]] != usem_pars$covD)
  
}






if (example) {
  
  library(gimme)
  data("simData", package = "gimme") 
  outGIMME <- gimme(data = simData)
  plot(outGIMME)
  
  
  
  simulateVAR(A = NULL,
              Phi = NULL,
              Psi = NULL,
              subAssign = NULL,
              N = 1,
              ASign = "random",
              PhiSign = "random",
              Obs = NULL,
              indA = 0.01,
              indPhi = 0.01,
              indPsi = 0.00)
  
}

