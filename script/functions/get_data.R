# Generate data with structural break in beta1
get.data <- function(n, phi, psi, b0, delta) {
  # Generate time points
  t <- seq(0, 1, length.out = n)
  
  # Generate AR(1) process
  x <- numeric(n)
  std_ar <- sqrt(1 / (1 - 0.3^2))  # Standard deviation of stationary AR(1)
  x[1] <- rnorm(1, mean = 0, sd = std_ar)
  
  # Generate AR(1) process with coefficient 0.3
  for (i in 2:n) {
    x[i] <- 0.3 * x[i - 1] + rnorm(1)
  }
  
  # Unconditional variance AR(1) process
  std_ar <- sqrt(1 / (1 - 0.3 ^ 2))
  
  # Generate ARMA(1,1) errors
  sigma_eps <- sqrt((1 - phi ^ 2) / (2 * (1 + psi ^ 2 + 2 * phi * psi)))
  eps <- rnorm(n, 0, sigma_eps)
  u <- numeric(n)
  u[1] <- eps[1]
  
  for (i in 2:n) {
    u[i] <- phi * u[i - 1] + eps[i] + psi * eps[i - 1]
  }
  
  # Introduce discontinuity
  b1 <- b0 + delta * std_ar

  # Define coefficient function β₁(·) with structural break at t = 0.5 size 0.1
  beta1 <- function(t) {
      ifelse(t <= 0.5, b0, b1)  #break
    }
  
  # Coefficient values DGP
  beta1_vals <- beta1(t)
  
  # Generate dependent variable y
  y <- numeric(n)
  for (i in 1:n) {
    y[i] <- beta1_vals[i] * x[i] + u[i]
  }
  
  return(list(
    y = y,
    time = t,
    beta1_vals = beta1_vals,
    x = x,
    u = u,
    beta_dgp = c(b0, b1)
  ))
}
