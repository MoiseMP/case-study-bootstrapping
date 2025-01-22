### Practical Case Data Science ###
# Authors:
# Emma Arussi, Julia Kroon, Moise Mpongo, Tess Scholtus
# Group 12

#### Clean ####
rm(list = ls())

### Bootstrapping with structural breaks ###

## Load packages ##
library(MASS)
library(foreach)
library(doParallel)
#library(np)
library(profvis)

# Set seed for reproducibility
set.seed(123)

# Starting time
start_time <- Sys.time()

# Generate data
n <- 300
phi <- 0.3
psi <- 0.3
b0 <- -0.1
b1 <- -0.9

# Number of Monte Carlo replications
M <- 1
B <- 10
alpha <- 0.05

# Generate data with structural break in beta1
get.data <- function(n, phi, psi, b0, b1) {
  # Generate time points
  t <- seq(0, 1, length.out = n)
  
  # Define coefficient function β₁(·) with structural break at t = 0.5 size 0.1
  beta1 <- function(t) {
    ifelse(t <= 0.5, b0, b1)  #break
  }
  
  # Generate AR(1) process
  x <- numeric(n)
  x[1] <- rnorm(1, 10, 1)
  
  # Generate AR(1) process with coefficient 0.3
  for(i in 2:n) {
    x[i] <- 0.3 * x[i-1] + rnorm(1)
  }
  
  # Generate ARMA(1,1) errors
  sigma_eps <- sqrt((1 - phi^2)/(2*(1 + psi^2 + 2*phi*psi)))
  eps <- rnorm(n, 0, sigma_eps)
  u <- numeric(n)
  u[1] <- eps[1]
  
  for(i in 2:n) {
    u[i] <- phi * u[i-1] + eps[i] + psi * eps[i-1]
  }
  
  # Coeff values DGP
  beta1_vals <- beta1(t)
  
  # Generate dependent variable y
  y <- numeric(n)
  for(i in 1:n) {
    y[i] <- beta1_vals[i] * x[i] + u[i]
  }
  
  return(list(
    y = y,
    time = t,
    beta1_vals = beta1_vals,
    x = x,
    u = u
  ))
}

sim_data <- get.data(n, phi, psi, b0, b1)

###########################
### Bandwidth Selection ###
###########################

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

# Calculate standard deviation of y for break sizes
sd_y <- sd(sim_data$y)
cat("Standard deviation of y:", sd_y, "\n")

# Define break sizes as multiples of sd(y)
break_sizes <- c(0.1, 0.5, 1) * sd_y

# Epanechnikov kernel
K <- function(u) {
  ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
}

# Generate data with structural break in beta1
get.data <- function(n, phi, psi, b0, b1) {
  # Generate time points
  t <- seq(0, 1, length.out = n)
  
  # Define coefficient function β₁(·) with structural break at t = 0.5 size 0.1
  beta1 <- function(t) {
    ifelse(t <= 0.5, b0, b1)  #break
  }
  
  # Generate AR(1) process
  x <- numeric(n)
  x[1] <- rnorm(1, 10, 1)
  
  # Generate AR(1) process with coefficient 0.3
  for(i in 2:n) {
    x[i] <- 0.3 * x[i-1] + rnorm(1)
  }
  
  # Generate ARMA(1,1) errors
  sigma_eps <- sqrt((1 - phi^2)/(2*(1 + psi^2 + 2*phi*psi)))
  eps <- rnorm(n, 0, sigma_eps)
  u <- numeric(n)
  u[1] <- eps[1]
  
  for(i in 2:n) {
    u[i] <- phi * u[i-1] + eps[i] + psi * eps[i-1]
  }
  
  # Coeff values DGP
  beta1_vals <- beta1(t)
  
  # Generate dependent variable y
  y <- numeric(n)
  for(i in 1:n) {
    y[i] <- beta1_vals[i] * x[i] + u[i]
  }
  
  return(list(
    y = y,
    time = t,
    beta1_vals = beta1_vals,
    x = x,
    u = u
  ))
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

# Calculate optimal bandwidths
bandwidths <- select_bandwidths(sim_data$y, sim_data$x, n)
h <- bandwidths$h
h_tilde <- bandwidths$h_tilde

########################
### Local Linear Fitting ###
########################

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

sim_data <- get.data(1000, phi, psi, b0, b1)
profvis(local_linear_fit(sim_data$y, sim_data$x, h = h))

########################
### Monte Carlo Simulation ###
########################

# Setup multiprocessing
cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)


# Monte Carlo loop
mc_results <- foreach(m = 1:M, .combine = rbind) %dopar% {
  
 
  
  # Generate new data for each Monte Carlo replication
  sim_data <- get.data(n, phi, psi, b0, b1)
  
  # Get bandwidths for this replication
  #bandwidths <- select_bandwidths(sim_data$y, sim_data$x, n)
  h <- 0.1
  h_tilde <- h^(5/9)
  
  # Bootstrap procedure
  initial_fit <- local_linear_fit(sim_data$y, sim_data$x, h = h)
  beta_hat <- initial_fit$beta_hat
  
  oversmooth_fit <- local_linear_fit(sim_data$y, sim_data$x, h = h,
                                     oversmooth = TRUE, h_tilde = h_tilde)
  beta_tilde <- oversmooth_fit$beta_hat
  
  # Calculate residuals
  z_hat <- numeric(n)
  for(t in 1:n) {
    z_hat[t] <- sim_data$y[t] - sim_data$x[t] * beta_tilde[t]
  }
  
  # Fit AR(p) model
  ar_fit <- ar(z_hat, aic = TRUE)
  p <- ar_fit$order
  phi_hat <- ar_fit$ar
  
  # Calculate centered residuals
  e_hat_p <- numeric(n-p)
  for(t in (p+1):n) {
    ar_part <- sum(phi_hat * z_hat[(t-1):(t-p)])
    e_hat_p[t-p] <- z_hat[t] - ar_part
  }
  e_tilde_p <- e_hat_p - mean(e_hat_p)
  
  # Bootstrap procedure
  beta_star <- matrix(0, n, B)
  
  for(b in 1:B) {
    # Bootstrap steps
    e_star <- sample(e_tilde_p, n, replace = TRUE)
    
    z_star <- numeric(n)
    z_star[1:p] <- z_hat[1:p]
    
    for(t in (p+1):n) {
      z_star[t] <- sum(phi_hat * z_star[(t-1):(t-p)]) + e_star[t]
    }
    
    y_star <- sim_data$x * beta_tilde + z_star
    
    boot_fit <- local_linear_fit(y_star, sim_data$x, h = h)
    
    beta_star[,b] <- boot_fit$beta_hat
  }
  
  # Calculate confidence bands
  beta_centered <- beta_star - beta_tilde
  #ci_beta <- beta_tilde + t(apply(beta_centered, 1, quantile, probs = c(alpha/2, 1 - alpha/2)))
  
  ci_beta <- matrix(0, n, 2)
  for(t in 1:n) {
    ci_beta[t,] <- beta_tilde[t] + quantile(beta_centered[t,], c(alpha/2, 1-alpha/2))
  }
  
  # Write results to data frame
  data.frame(
    coverage = mean(sim_data$beta1_vals >= ci_beta[,1] & 
                         sim_data$beta1_vals <= ci_beta[,2]),
    mse = mean((sim_data$beta1_vals - beta_tilde)^2),
    h = h,
    h_tilde = h_tilde,
    p = p
  )
  
  
}

# Exit Cluster
stopCluster(cl)

# Calculate Monte Carlo summary statistics

# Write output to a file
filename = paste0('monte_carlo_result_',format(start_time, '%d_%m_%Y_%H-%M-%S'), '.txt')
file_name_plots <- paste0("monte_carlo_plots_", format(start_time, "%d_%m_%Y_%H-%M-%S"), ".png")
sink(filename)
cat("Results Monte Carlo Anaylsis on date:", format(start_time, '%d %m %Y %H:%M:%S'))

cat('\n=== DGP Setup ====')
cat('\nn=',n)
cat('\nphi=', phi)
cat('\npsi=', psi)

cat("\n\n=== Monte Carlo and Bootstrapping Setup ===")
cat('\nM =', M)
cat('\nB =', B)
cat('\nalpha =', alpha)
cat('\nbeta0 =', b0)
cat('\nbeta1 =', b1)

cat("\n\n=== Monte Carlo Results ===")
cat("\nAverage Coverage Probability:", round(mean(mc_results$coverage), 3))
cat("\nStd Dev of Coverage:", round(sd(mc_results$coverage)), 3)
cat("\nAverage MSE:", round(mean(mc_results$mse), 3))
cat("\nStd Dev of MSE:", round(sd(mc_results$mse), 3))
cat("\nAverage Selected h:", round(mean(mc_results$h), 4))
cat("\nAverage Selected h_tilde:", round(mean(mc_results$h_tilde), 4))
cat("\nAverage Selected AR order:", round(mean(mc_results$p), 2))

cat("\nCorresponding plots can be found at: ", file_name_plots)

sink()

cat('\nResults written to file: ', filename)

# Save plots to file
png(file_name_plots, width=3000, height=2000, res = 300) 

# Optional: Plot Monte Carlo results
par(mfrow=c(2,2))

# Coverage probability distribution
hist(mc_results$coverage, main="Distribution of Coverage Probabilities",
     xlab="Coverage Probability", breaks=20)
abline(v=0.95, col="red", lty=2)

# MSE distribution
hist(mc_results$mse, main="Distribution of MSE",
     xlab="Mean Squared Error", breaks=20)

# Bandwidth distributions
hist(mc_results$h, main="Distribution of h",
     xlab="Selected Bandwidth h", breaks=20)
hist(mc_results$h_tilde, main="Distribution of h_tilde",
     xlab="Selected Bandwidth h_tilde", breaks=20)

par(mfrow=c(1,1))

dev.off()
cat("Plots saved to", file_name_plots, "\n")
cat("Done, total time taken: " , Sys.time() - start_time)
