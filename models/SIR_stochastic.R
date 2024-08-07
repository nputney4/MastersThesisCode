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