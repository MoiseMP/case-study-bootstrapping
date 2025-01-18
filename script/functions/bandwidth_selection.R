library(np)

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




