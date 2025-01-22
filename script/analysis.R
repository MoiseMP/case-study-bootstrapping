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
source(here("script", "functions", "get_data.R"))
source(here("script", "functions", "estimation.R"))

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
delta <- break_sizes[3]

sim_data <- get.data(n, phi, psi, b0, delta)
b1 <- sim_data$beta_dgp[2]

# Number of Monte Carlo replications
M <- 5
B <- 5
alpha <- 0.05

# Bandwidth Selection
h_vector <- c(0.06, 0.15, 0.24)
h_tilde_vector <- h_vector ^ (5/9)
n_loops <- length(h_vector)


########################
### Monte Carlo Simulation ###
########################
print("Starting Monte Carlo simulation")

for (i in 1:n_loops) {
  
  cat("H", i, "out of", n_loops, '\n')
  
  h <- h_vector[i]
  h_tilde <- h_tilde_vector[i]
  
  # Setup multiprocessing
  cores <- detectCores()
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  
  # Monte Carlo loop
  mc_results <- foreach(m = 1:M, .combine = rbind) %dopar% {
    
    # Generate new data for each Monte Carlo replication
    sim_data <- get.data(n, phi, psi, b0, delta)
    
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
  
  ##############
  # plots of coverage probabilities
  ##############
  
  file_name_plots <- paste0("n", n, 
                            "_h", round(h, 3), 
                            "_M", M, 
                            "_B", B, 
                            "_delta", round(delta, 1), 
                            "_", format(start_time, "%d_%m_%Y_%H-%M-%S"), 
                            ".png")
  
  file_name_coverage <- here("output", "plots", paste0("coverage_failures_", file_name_plots))
  file_name_plot_distributions_coverage <- here("output", "plots", paste0("monte_carlo_plots_", file_name_plots))
  file_name_beta_plot <- here("output", "plots", paste0("beta_plot_", file_name_plots))
  filename_results <- here("output", "data", paste0("monte_carlo_result_", format(start_time, "%d_%m_%Y_%H-%M-%S"), ".txt"))
  
  # Interactive plots in R
  par(mfrow=c(2,2))
  
  # Plot 1: Distribution of Coverage probabilities, MSE, h and h_tilde
  
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
  
  
  # Save plot to file
  png(file_name_plot_distributions_coverage, width=3000, height=2000, res=300)
  
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
  cat("Plots saved to", file_name_plot_distributions_coverage, "\n")
  
  ####################
  # plots of estimated beta + CI
  ###################
  
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
  png(file_name_beta_plot, width = 3000, height = 2000, res = 300)
  
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
  
  ##############
  # Failure rate 
  ##############
  
  # Create a matrix to store coverage failures at each time point
  coverage_failures <- matrix(0, nrow = M, ncol = n)
  
  # For each Monte Carlo replication
  for(m in 1:M) {
    beta_tilde <- unlist(mc_results$beta_tilde[m])       # Estimated beta from replication
    ci_beta <- mc_results$ci_beta[[m]]                  # Confidence intervals
    lower_ci <- as.numeric(ci_beta[,1]$V1)             # Lower CI
    upper_ci <- as.numeric(ci_beta[,2]$V2)             # Upper CI
    
    # Check coverage at each time point
    coverage_failures[m,] <- !(sim_data$beta1_vals >= lower_ci & 
                                 sim_data$beta1_vals <= upper_ci)
  }
  
  # Calculate failure rate at each time point
  failure_rates <- colMeans(coverage_failures)          # Mean failure rates over replications
  
  # --- Interactive Plot: Failure Rates ---
  plot(sim_data$time, failure_rates,
       type = "l",
       col = "blue",
       lwd = 2,
       main = sprintf("Coverage Failure Rates Over Time\nn=%d, h=%.3f, M=%d, B=%d, delta=%.1f",
                      n, h, M, B, delta),
       xlab = "Time",
       ylab = "Failure Rate",
       ylim = c(0, 1))                                  # Failure rate between 0 and 1
  
  # Add expected failure rate (alpha) and structural break point
  abline(h = alpha, col = "red", lty = 2)               # Expected failure rate (alpha)
  abline(v = 0.5, col = "green", lty = 2)               # Structural break 
  
  # Add legend
  legend("topright",
         legend = c("Failure Rate", "Expected Rate (alpha)", "Structural Break"),
         col = c("blue", "red", "green"),
         lty = c(1, 2, 2),
         lwd = c(2, 1, 1))
  
  # --- Save Plot to File ---
  png(file_name_coverage, width = 3000, height = 2000, res = 300)
  
  # Re-plot for saving
  plot(sim_data$time, failure_rates,
       type = "l",
       col = "blue",
       lwd = 2,
       main = sprintf("Coverage Failure Rates Over Time\nn=%d, h=%.3f, M=%d, B=%d, delta=%.1f",
                      n, h, M, B, delta),
       xlab = "Time",
       ylab = "Failure Rate",
       ylim = c(0, 1))
  abline(h = alpha, col = "red", lty = 2)
  abline(v = n/2, col = "green", lty = 2)
  legend("topright",
         legend = c("Failure Rate", "Expected Rate (alpha)", "Structural Break"),
         col = c("blue", "red", "green"),
         lty = c(1, 2, 2),
         lwd = c(2, 1, 1))
  
  # --- Save Monte Carlo and Coverage Failure Results to File ---
  sink(filename_results, split = TRUE)
  
  ###################################
  # --- Monte Carlo Results ---
  cat("Results Monte Carlo Analysis on date:", format(start_time, '%d %m %Y %H:%M:%S'), "\n")
  
  cat('\n=== DGP Setup ====\n')
  cat('n =', n, "\n")
  cat('phi =', phi, "\n")
  cat('psi =', psi, "\n")
  cat('beta0 =', b0, "\n")
  cat('beta1 =', b1, "\n")
  
  cat("\n=== Monte Carlo and Bootstrapping Setup ===\n")
  cat('M (Monte Carlo Replications) =', M, "\n")
  cat('B (Bootstrap Samples) =', B, "\n")
  cat('alpha (Significance Level) =', alpha, "\n")
  
  cat("\n=== Monte Carlo Results ===\n")
  cat("Average Coverage Probability:", round(mean(mc_results$coverage), 3), "\n")
  cat("Std Dev of Coverage:", round(sd(mc_results$coverage), 3), "\n")
  cat("Average MSE:", round(mean(mc_results$mse), 3), "\n")
  cat("Std Dev of MSE:", round(sd(mc_results$mse), 3), "\n")
  cat("Average Selected h:", round(mean(mc_results$h), 4), "\n")
  cat("Average Selected h_tilde:", round(mean(mc_results$h_tilde), 4), "\n")
  cat("Average Selected AR order:", round(mean(mc_results$p), 2), "\n")
  
  # --- Coverage Failure Characteristics ---
  cat("\n=== Coverage Failure Analysis ===\n")
  cat("Average failure rate:", mean(failure_rates), "\n")           # Mean failure rate
  cat("Maximum failure rate:", max(failure_rates), "\n")            # Maximum failure rate
  cat("Time point of maximum failure rate:", which.max(failure_rates), "\n") # Time of max failure
  cat("Failure rate at break point (n/2):", failure_rates[round(n/2)], "\n") # Failure rate at break
  
  cat("\nCorresponding plots can be found at:\n")
  cat("Beta plot:", file_name_beta_plot, "\n")
  cat("Coverage failure plot:", file_name_coverage, "\n")
  
  cat("Total time taken:", Sys.time() - start_time, "\n")
  
  sink()
  
  # Final message in console
  cat('\nResults written to file:', filename_results, "\n")
}