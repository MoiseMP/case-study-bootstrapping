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
library(tidyr)
library(here)

# Load custom functions
source(here("Desktop", "rstudio", "get_data.R"))
source(here("Desktop", "rstudio", "estimation.R"))

# Set seed for reproducibility
set.seed(123)

# Starting time
start_time <- Sys.time()

# Generate data
n <- 200
phi <- 0.3
psi <- 0.3
b0 <- -0.1

# Break sizes are determined as constant times long run variance of independent variable x
break_sizes <- c(0.1, 0.5, 1)
delta <- break_sizes[1]

sim_data <- get.data(n, phi, psi, b0, delta)
b1 <- sim_data$beta_dgp[2]
plot(sim_data$time, sim_data$beta1_vals)

# Number of Monte Carlo replications
M <- 500
B <- 500
alpha <- 0.05

###########################
### Bandwidth Selection ###
###########################

# Calculate optimal bandwidths
h <- 0.6
h_tilde <- h ^ (5/9)

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
  sim_data <- get.data(n, phi, psi, b0, delta)
  
  # Get bandwidths for this replication
  h_tilde <- h^(5/9)
  
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
  tidyr::tibble(
    coverage = mean(sim_data$beta1_vals >= ci_beta[,1] & 
                      sim_data$beta1_vals <= ci_beta[,2]),
    mse = mean((sim_data$beta1_vals - beta_tilde)^2),
    h = h,
    h_tilde = h_tilde,
    p = p,
    beta_tilde = list(beta_tilde),
    ci_beta = list(tidyr::as_tibble(ci_beta))
  )
}

# Exit Cluster
stopCluster(cl)

# Calculate Monte Carlo summary statistics
# Write output to a file
filename <- here("Desktop", "rstudio", paste0("monte_carlo_result_", format(start_time, "%d_%m_%Y_%H-%M-%S"), ".txt"))
file_name_plots <- here("Desktop", "rstudio", paste0("monte_carlo_plots_", format(start_time, "%d_%m_%Y_%H-%M-%S"), ".png"))

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

# Interactive plots in R
par(mfrow=c(2,2))

hist(mc_results$coverage, main="Distribution of Coverage Probabilities",
     xlab="Coverage Probability", breaks=20)
abline(v=0.95, col="red", lty=2)

hist(mc_results$mse, main="Distribution of MSE",
     xlab="Mean Squared Error", breaks=20)

hist(mc_results$h, main="Distribution of h",
     xlab="Selected Bandwidth h", breaks=20)
hist(mc_results$h_tilde, main="Distribution of h_tilde",
     xlab="Selected Bandwidth h_tilde", breaks=20)

par(mfrow=c(1,1))

# Save to file
png(file_name_plots, width=3000, height=2000, res=300)

par(mfrow=c(2,2))

hist(mc_results$coverage, main="Distribution of Coverage Probabilities",
     xlab="Coverage Probability", breaks=20)
abline(v=0.95, col="red", lty=2)

hist(mc_results$mse, main="Distribution of MSE",
     xlab="Mean Squared Error", breaks=20)

hist(mc_results$h, main="Distribution of h",
     xlab="Selected Bandwidth h", breaks=20)
hist(mc_results$h_tilde, main="Distribution of h_tilde",
     xlab="Selected Bandwidth h_tilde", breaks=20)

par(mfrow=c(1,1))

dev.off()
cat("Plots saved to", file_name_plots, "\n")

# Extract the first Monte Carlo replication's results for visualization
first_rep_beta <- unlist(mc_results$beta_tilde[1])  # Estimated beta
first_rep_ci <- mc_results$ci_beta[[1]]            # Confidence intervals

# Ensure confidence intervals are numeric vectors
lower_ci <- as.numeric(first_rep_ci[,1]$V1)        # Lower confidence interval
upper_ci <- as.numeric(first_rep_ci[,2]$V2)        # Upper confidence interval

# --- Interactive Plot in R ---
# Add n, h, M, and B to the title
plot(sim_data$time,                                 
     sim_data$beta1_vals,                          
     type = "l",                                   
     col = "black",                                
     lwd = 2,                                      
     ylim = range(c(lower_ci,                     
                    upper_ci,
                    sim_data$beta1_vals,
                    first_rep_beta)),
     main = paste("True Beta, Estimated Beta, and Confidence Intervals\n",
                  "n =", n, ", h =", round(h, 3), ", M =", M, ", B =", B),
     xlab = "Time",
     ylab = "Beta")

# Add estimated beta
lines(sim_data$time, 
      first_rep_beta, 
      col = "blue", 
      lwd = 2,
      lty = 2)  # dashed line

# Add confidence intervals
lines(sim_data$time, lower_ci, col = "red", lty = 3)    # lower CI
lines(sim_data$time, upper_ci, col = "red", lty = 3)    # upper CI

# Add legend
legend("bottomright", 
       legend = c("True Beta", "Estimated Beta", "95% CI"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 3),
       lwd = c(2, 2, 1))

# --- Save Plot to File ---
file_name_beta_plot <- here("Desktop", "rstudio", paste0("beta_plot_", format(start_time, "%d_%m_%Y_%H-%M-%S"), ".png"))
png(file_name_beta_plot, width = 3000, height = 2000, res = 300)

# Same plot saved to file
# Add n, h, M, and B to the title
plot(sim_data$time,                                 
     sim_data$beta1_vals,                          
     type = "l",                                   
     col = "black",                                
     lwd = 2,                                      
     ylim = range(c(lower_ci,                     
                    upper_ci,
                    sim_data$beta1_vals,
                    first_rep_beta)),
     main = paste("True Beta, Estimated Beta, and Confidence Intervals\n",
                  "n =", n, ", h =", round(h, 3), ", M =", M, ", B =", B),
     xlab = "Time",
     ylab = "Beta")

lines(sim_data$time, first_rep_beta, col = "blue", lwd = 2, lty = 2)
lines(sim_data$time, lower_ci, col = "red", lty = 3)    # lower CI
lines(sim_data$time, upper_ci, col = "red", lty = 3)    # upper CI

legend("bottomright", 
       legend = c("True Beta", "Estimated Beta", "95% CI"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 3),
       lwd = c(2, 2, 1))

dev.off()

cat("Beta plot saved to", file_name_beta_plot, "\n")
