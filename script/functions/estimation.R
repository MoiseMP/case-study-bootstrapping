library(np)

# Epanechnikov kernel
K <- function(u) {
  ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
}

# Fitting TVP model using WLS
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
    
    # Weighted least squares estimation
    theta <- lm.wfit(x = X, y = y, w = kernel_weights)$coefficients
    
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

# Main bandwidth selection function
select_bandwidths <- function(y, x, n) {
  h_candidates <- seq(0.2, 0.3, by = 0.001)
  
  # MCV approach
  h <- mcv_bandwidth(y, x, h_candidates)
  
  # Calculate h_tilde using Bühlmann's approach
  h_tilde <- get_buhlmann_h_tilde(h, n)
  
  cat("\nBandwidth Selection Results:\n")
  cat("MCV bandwidth (h):", round(h, 4), "\n")
  
  return(list(h = h, h_tilde = h_tilde))
}

# Bühlmann (1998) oversmoothing bandwidth
get_buhlmann_h_tilde <- function(h, n, C = 1) {
  h_tilde <- C * h^(5/9)
  check_B1 <- max(h_tilde, n*h_tilde*h^4, h*log(n)/h_tilde)
  
  cat("Bühlmann oversmoothing bandwidth check:\n")
  cat("h_tilde:", round(h_tilde, 4), "\n")
  cat("Assumption B1 check value:", round(check_B1, 4), "\n")
  
  return(h_tilde)
}

# Modified Cross-Validation (MCV)
mcv_bandwidth <- function(y, x, h_candidates, l = 1) {
  n <- length(y)
  tau <- seq(1, n, 1)/n
  mcv_scores <- numeric(length(h_candidates))
  
  for(h_idx in seq_along(h_candidates)) {
    h <- h_candidates[h_idx]
    mcv_sum <- 0
    
    for(t in 1:n) {
      indices <- which(abs(seq_len(n) - t) > l)
      u <- (tau[indices] - tau[t])/h
      kernel_weights <- K(u)
      X <- cbind(x[indices], x[indices] * (tau[indices] - tau[t]))
      W <- diag(kernel_weights)
      
      theta <- lm.wfit(x = X, y = y[indices], w = diag(W))$coefficients

      y_hat <- x[t] * theta[1]
      mcv_sum <- mcv_sum + (y[t] - y_hat)^2
    }
    mcv_scores[h_idx] <- mcv_sum/n
  }
  
  return(h_candidates[which.min(mcv_scores)])
}




