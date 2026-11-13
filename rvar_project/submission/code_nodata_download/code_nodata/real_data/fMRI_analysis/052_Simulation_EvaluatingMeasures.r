##################################################
##################################################
##################################################
eval_forecast <- function(Y_forecast, B_est_list, range, horizon) {
  
  n <- length(Y_forecast)
  
  #msfe <- numeric(horizon)
  msfe <- matrix(0, ncol = n, nrow = horizon)
  for (h in 1:horizon) {
    
    bysubj_msfe <- numeric(n)
    for (k in 1:n) {
      Y <- Y_forecast[[k]]
      B <- B_est_list[[k]] %^% h
      
      Z <- Y[h + 1:(range - h),]
      X <- Y[1:(range - h),]
      
      err_mat <- Z - X %*% t(B)
      
      bysubj_msfe[k] <- sqrt(mean(err_mat^2))
      
    }
    msfe[h,] <- bysubj_msfe
    
  }
  
  rownames(msfe) <- paste0("msfe_step", 1:h)
  
  return(msfe)
  
}