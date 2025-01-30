# One off script to generate the the DGP plots

source("script/functions/get_data.R")

n <- 200
phi <- 0.3
psi <- 0.3
b0 <- -0.1

delta <- c(0.1, 0.5, 1)

f1 <- get.data(n, phi, psi , b0, delta[1])
f2 <- get.data(n, phi, psi , b0, delta[2])
f3 <- get.data(n, phi, psi , b0, delta[3])


n_vec <- seq(0, 1, length.out = n)

png("output/plots/dgp_plots.png", width = 3000, height = 2000, res = 300)

plot(n_vec, f3$beta1_vals, 
     type = 'l', 
     col = 'red',
     lwd = 2,
     xlab = "Time",
     ylab = expression(beta[t]),)

lines(n_vec, f2$beta1_vals, col = 'green', lwd = 2)
lines(n_vec, f1$beta1_vals, col ='blue', lwd = 2)

legend("topleft",
       legend = c(expression(lambda == 0.1), expression(lambda == 0.5), expression(lambda == 1)),
       col = c("blue", "red", "green"),
       lty = c(1, 2, 2),
       lwd = c(2, 1, 1))

dev.off()

