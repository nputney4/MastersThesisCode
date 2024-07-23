SIR_continuous <- odin::odin({
  ## This code is not R code!
  deriv(S) <- - (beta * I/N) * S 
  deriv(I) <- (beta * I/N) * S - gamma * I
  deriv(R) <- gamma * I 
  
  beta <- user()
  gamma <- user()
  
  initial(S) <- 1000
  initial(I) <- 1
  initial(R) <- 0
  
  N = S + I + R
})

## Running the model
mod_SIR_continuous <- SIR_continuous$new(beta = 0.3, gamma = 0.2)

mod_SIR_continuous$run(1:20)