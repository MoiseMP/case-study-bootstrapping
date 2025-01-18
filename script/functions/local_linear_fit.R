local_linear_fit <- function(y, x, h, oversmooth = FALSE, h_tilde = NULL) {
  n <- length(y)
  tau <- seq(1, n, 1)/n
  y_hat <- numeric(n)
  beta_hat <- numeric(n)
  
  # Use appropriate bandwidth
  h_effective <- if(oversmooth) h_tilde else h
  
  # For each time point
  for(i in 1:n) {
    # Calculate kernel weights
    u <- (tau - tau[i])/h_effective
    kernel_weights <- K(u)
    
    # Create weighted design matrix
    X <- cbind(x, x * (tau - tau[i]))
    W <- diag(kernel_weights)
    
    # Weighted least squares estimation
    theta <- lm.wfit(x = X, y = y, w = diag(W))$coefficients

    beta_hat[i] <- theta[1]
    y_hat[i] <- x[i] * beta_hat[i]
  }
  
  residuals <- y - y_hat
  
  return(list(
    y_hat = y_hat,
    residuals = residuals,
    time = tau,
    beta_hat = beta_hat
  ))
}