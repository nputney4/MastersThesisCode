SIR_discrete <- odin.dust::odin_dust({
  ## Individual probabilities of transition:
  p_SI <- 1 - exp(-beta * I / N) # S to I
  p_IR <- 1 - exp(-gamma) # I to R
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SI <- rbinom(S, p_SI)
  n_IR <- rbinom(I, p_IR)
  
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  
  beta <- user()
  gamma <- user()
  
  initial(S) <- 1000
  initial(I) <- 1
  initial(R) <- 0
  
  N = S + I + R
})

## Running the model
mod_SIR_discrete <- SIR_discrete$new(pars = list(beta = 0.4, gamma = 0.1),
                                     time = 1,
                                     n_particles = 10)
n_times <- 200
x <- mod_SIR_discrete$simulate(1:n_times)

# Plotting particles for each compartment
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
time <- seq(1, n_times)
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x))
matlines(time, t(x[2, , ]), col = cols[["I"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["R"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")