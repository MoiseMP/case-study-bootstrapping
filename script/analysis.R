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
library(np)

# Load custom functions
source('script/functions/get_data.R')
source('script/functions/estimation.R')

# Set seed for reproducibility
set.seed(123)

# Starting time
start_time <- Sys.time()

# Generate data
n <- 300
phi <- 0.3
psi <- 0.3
b0 <- -0.1
break_sizes <- c(0.1, 0.5, 1)
sim_data <- get.data(n, phi, psi, b0, break_sizes[1])

# Number of Monte Carlo replications
M <- 5
B <- 100
alpha <- 0.05

###########################
### Bandwidth Selection ###
###########################

# Calculate optimal bandwidths
bandwidths <- select_bandwidths(sim_data$y, sim_data$x, n)
h <- bandwidths$h
h_tilde <- bandwidths$h_tilde

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
  bandwidths <- select_bandwidths(sim_data$y, sim_data$x, n)
  h <- bandwidths$h
  h_tilde <- bandwidths$h_tilde
  
  # Bootstrap procedure
  initial_fit <- local_linear_fit(sim_data$y, sim_data$x, h = h)
  beta_hat <- initial_fit$beta_hat
  
  oversmooth_fit <- local_linear_fit(sim_data$y, sim_data$x, h = h,
                                     oversmooth = TRUE, h_tilde = h_tilde)
  beta_tilde <- oversmooth_fit$beta_hat
  
  # Calculate residuals
  z_hat <- sim_data$y - sim_data$x * beta_tilde
  
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
filename = paste0('output/data/monte_carlo_result_',format(start_time, '%d_%m_%Y_%H-%M-%S'), '.txt')
file_name_plots <- paste0("output/plots/monte_carlo_plots_", format(start_time, "%d_%m_%Y_%H-%M-%S"), ".png")
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

cat("\nCorresponding plots can be found at:", file_name_plots)
cat("\nTotal time taken:" , Sys.time() - start_time)

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

