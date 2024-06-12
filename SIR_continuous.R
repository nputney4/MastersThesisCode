########################################################################################
######## Estimating parameters of deterministic SIR using pMCMC with Odin-Dust #########
########################################################################################
# R version 4.10
library(dust)
library(odin.dust)
library(mcstate)
library(ggplot2)
library(reshape2) # for melt function
library(coda) # for assessing markov chains

## First step is defining the model using the syntax of odin 
## which will then be parsed and compiled into C++ 
gen_odin_det <- odin::odin({
  
  ## Defining first order differential equations
  deriv(S) <- - (beta * I/N) * S 
  deriv(I) <- (beta * I/N) * 
    S - gamma * I
  deriv(R) <- gamma * I 
  
  ## Model parameters - user defined when creating object from generator
  beta <- user()
  gamma <- user()
  
  ## Initial conditions - user defined, defaults in parenthesis
  initial(S) <- 1000
  initial(I) <- 1
  initial(R) <- 0
  
  ## Total population size
  N = S + I + R
  
})

mod1 <- gen_odin_det$new(beta=0.3, gamma=0.2) # ensure beta / gamma > 1 so we exit equilibrium

mod1$contents() # inspect parameters and initial states

## Determining the state of the system after 50 days with time steps of 0.25 days
times <- seq(from = 0, to = 150, by = 1)
SIR_results <- as.data.frame(mod1$run(times))
SIR_results_melted <- melt(SIR_results, 
                           id = "t", 
                           variable.name = c("Compartment"))

## Plotting S, I, R compartments over time
ggplot(data=SIR_results_melted, aes(x = t, y = value, group=Compartment)) +
  geom_line(aes(color=Compartment), size=1.5) + 
  xlab("Days") + ylab("Number of People") + 
  ggtitle("Simulation of Epidemic Using Simple SIR Model") + 
  theme(plot.title = element_text(hjust=0.5, size=16)) +
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.title.y = element_text(size=14))

## Number of new infected cases (number of people leaving the susceptible compartment) 
## Creating simulated observed data - could add noise later
set.seed(10) # for reproducibility
#cases_day <- c(0, -diff(SIR_results$S)) + rgamma(max(times)+1, 0.5, 1)
cases_day <- c(0, -diff(SIR_results$S))
observed_data <- as.data.frame(cbind(0:(length(cases_day)-1), round(cases_day, 0)))

colnames(observed_data) <- c("t", "n_cases")

## Plotting "observed" data
ggplot(data = observed_data, aes(x = t, y = n_cases)) +
  geom_line() + geom_point(color = "blue", size = 2, shape = 1) +  
  xlab("Days") + ylab("New Cases") + 
  ggtitle("Simulated Cases Per Day Using SIR Model")+ 
  theme(plot.title = element_text(hjust=0.5, size=16)) +
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.title.y = element_text(size=14))

## Define new model with unknown beta and gamma and try to estimate 
## parameters using mcstate package - let's find out how (idk)!

## Note - "mcstate is primarily designed to work with discrete-time stochastic models written
## in odin.dist ... it is often possible to interact with continuous deterministic models written in
## odin." - https://mrc-ide.github.io/mcstate/articles/continuous.html

## Must use odin.dust to compile model in order to estimate parameters with mcstate
## Note that incidence or other model-derived metrics that will be compared 
## to observations for inference with mcstate must be defined within the model

gen_dust_det <- odin.dust::odin_dust({
  
  ## Defining first order differential equations
  deriv(S) <- - (beta * I/N) * S
  deriv(I) <- (beta * I/N) * S - gamma * I
  deriv(R) <- gamma * I
  output(n_cases) <- beta * I/N * S
  
  ## Model parameters - user defined when creating object from generator
  beta <- user()
  gamma <- user()
  b.A <- 0.5
  ## Initial conditions - user defined, defaults in parenthesis
  initial(S) <- 1000
  initial(I) <- 1
  initial(R) <- 0
  
  ## Total population size
  N = S + I + R
})

dev.off()

n_particles = 1
mod2 <- gen_dust_det$new(pars=list(beta = 0.4, gamma = 0.3),
                         time=1,
                         n_particles=n_particles)

n_times <- max(times)
x <- array(NA, dim = c(mod2$info()$len, n_particles, n_times))

for (t in seq_len(n_times)) {
  x[ , , t] <- mod2$run(t)
}
S <- x[1,,]
cases <- x[4,,]


plot(0:(n_times), observed_data$n_cases, col="green", lwd=2,
     main="Daily Incidence of SIR Model",
     xlab="time", ylab="New Cases", type = "l")
lines(0:(n_times-1), cases, col = "red", lwd=2)
legend("topright", legend = c("observed", "model"), 
       col = c("green", "red"), lty = 1, lwd = 2)


# COMPARISON FUNCTION:
# Must take arguments state, observed and pars as arguments 
# (though the arguments may have different names). state is the simulated model state 
# (a matrix with as many rows as there are state variables and as many columns as 
# there are particles, data is a list of observed data corresponding to the current 
# time's row in the data object provided here in the constructor. pars is any additional 
# parameters passed through to the comparison function (via the pars argument to $run). 


# note that 'prev_state' is no longer supported even though still found in their tutorials
# this means incidence must be calculated within the model itself

incidence_compare <- function(state, observed, pars = NULL){
  exp_noise <- 1e6
  incidence_observed <- observed$n_cases
  incidence_modelled <- state[4,,drop=TRUE] 
  lambda <- incidence_modelled + rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}

## Put observed data into format usable for mcstate
filt_data <- mcstate::particle_filter_data(data = observed_data,
                              time = "t",
                              initial_time = 0, rate = NULL)
head(filt_data)

## Currently not too sure of the purpose of this - apparently optional
index <- function(info) {
  idx <- unlist(info$index)
  list(run = idx, state = idx)
}

## Alternatively, one could define the comparison function inside the model
## in which case, compare would be set to 'NULL'
filter <- mcstate::particle_deterministic$new(data = filt_data, model = gen_dust_det,
                                              index = index, compare = incidence_compare)

#filter$run(save_history = TRUE, pars=list(beta = 0.4, gamma = 0.3)) 

beta <- mcstate::pmcmc_parameter("beta", initial = 0.5, min = 0, max = 2)
gamma <- mcstate::pmcmc_parameter("gamma", initial = 0.5, min = 0, max = 2)

## This will from a multivariate normal with the specified cov matrix - 
## adjust cov based on acceptance rate (or use adaptive_proposal which is giving error sometimes)
proposal_matrix <- diag(0.0001, 2)

mcmc_pars <- mcstate::pmcmc_parameters$new(list(beta = beta, gamma = gamma), 
                                           proposal_matrix)
#mcmc_pars$fix(c("gamma" = 0.2))

control <- mcstate::pmcmc_control(n_steps = 5e3, n_burnin = 500, progress = TRUE,
                                  n_chains = 4) # this parallelizes over particles by default (~ 3.5 mins)

control <- mcstate::pmcmc_control(n_steps = 10e3, n_burnin = 5000, progress = TRUE,
                                  n_chains = 4, n_workers = 4, n_threads_total = 8)
                                  #adaptive_proposal = TRUE) # this parallelizes over chains (~ 36 seconds)


mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)
mcmc_sample <- mcstate::pmcmc_sample(mcmc_run, n_sample = 500)

# Calculating acceptance rate and summary of posterior
coda_pars <- as.mcmc(cbind(mcmc_run$probabilities, mcmc_run$pars))
acc_rate <- 1 - rejectionRate(coda_pars)
sum_post <- summary(coda_pars)
med_beta <- sum_post$quantiles[4,3]
med_gamma <- sum_post$quantiles[5,3] 
c1_beta <- sum_post$quantiles[4,1]
c2_beta <- sum_post$quantiles[4,5]
c1_gamma <- sum_post$quantiles[5,1]
c2_gamma <- sum_post$quantiles[5,5]

# Assessing convergence of parameter estimates
par(mfrow=c(1,2))
#traceplot(coda_pars)
plot(mcmc_run$pars[, "beta"], type = "l", xlab = "Iteration",
     ylab = expression(beta))
plot(mcmc_run$pars[, "gamma"], type = "l", xlab = "Iteration",
     ylab = expression(gamma))

# Plotting posterior
par(mfrow = c(1,2))
hist(mcmc_sample$pars[, "beta"], main = "", xlab = expression(beta),
     freq = FALSE, breaks = "fd")
abline(v = mod1$contents()$beta, lty = 1, col = "darkred", lwd = 3)
abline(v = c(med_beta, c1_beta, c2_beta), lty = 2, 
       col = c("black", "blue", "blue"), lwd = 3)
hist(mcmc_sample$pars[, "gamma"], main = "", xlab = expression(gamma),
     freq = FALSE, breaks = "fd")
abline(v = mod1$contents()$gamma, lty = 1, col = "darkred", lwd = 3)
abline(v = c(med_gamma, c1_gamma, c2_gamma), lty = 2, 
       col = c("black", "blue", "blue"), lwd = 3)
legend("topright", inset=c(0,0.15), 
       legend = c("True value", "Median", "2.5/97.5 Quants"), 
       col = c("darkred", "black", "blue"), 
       lty = c(1,2,2), bty = "o", cex=0.5)

##-------------------- MLE OF BETA AND GAMMA ---------------------------##
cases_obs <- observed_data$n_cases

likelihood <- function(params, cases_obs){
  times <- 80
  particles <- 1
  model <- gen_dust_det$new(pars=list(beta = params[1], gamma = params[2]),
                            time=0,
                            n_particles = particles)
  
  states <- array(NA, dim = c(model$info()$len, particles, times))
  for (t in seq_len(times)) {
    states[ , , t] <- model$run(t)
  }
  cases_pred <- states[4,,]
  neg_ll <- -sum(dpois(x = cases_obs, lambda = cases_pred, log = TRUE))
  return(neg_ll)
}

opt_params <- optim(par = c(0.5, 0.5), fn = likelihood, cases_obs = cases_obs)

## True 'beta' is 0.3, true 'gamma' is 0.2
print(paste("MLE beta:", round(opt_params$par[1], 3), ";", 
            "MLE gamma:", round(opt_params$par[2], 3)))

##---------- REDEFINING MODEL USING WEEKLY INCIDENCE (in progress) ------------##
# Is a change even needed? A test!

## Taking every 7th observation as simulated weeky data
observed_weekly <- observed_data[observed_data$t %% 7 == 0,]
ggplot(data = observed_weekly, aes(x = t, y = n_cases)) +
  geom_line() + geom_point(color = "blue", size = 3, shape = 1) +  
  xlab("Days") + ylab("New Cases") + 
  ggtitle("Simulated Cases Per Week Using SIR Model")+ 
  theme(plot.title = element_text(hjust=0.5, size=16)) +
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.title.y = element_text(size=14))

## Inference with weekly data
filt_data2 <- mcstate::particle_filter_data(data = observed_weekly,
                                           time = "t",
                                           initial_time = 0, rate = NULL)
head(filt_data2)
filter2 <- mcstate::particle_deterministic$new(data = filt_data2, model = gen_dust_det,
                                              index = index, compare = incidence_compare)

control2 <- mcstate::pmcmc_control(n_steps = 10e3, n_burnin = 5000, progress = TRUE,
                                  n_chains = 4, n_workers = 4, n_threads_total = 8)


mcmc_run2 <- mcstate::pmcmc(mcmc_pars, filter2, control = control2)
mcmc_sample2 <- mcstate::pmcmc_sample(mcmc_run2, n_sample = 500)

# Looking at summary of marginal posteriors
coda_pars2 <- as.mcmc(cbind(mcmc_run2$probabilities, mcmc_run2$pars))
summary(coda_pars2)

## ---------------- INCORPORATING FORCING INTO THE MODEL ---------------- ##
## Same as before but now we say that the infectivity, beta, is 
## function of two parameters, one of which is time varying.
## Maybe could think like 'rate' is rainfall and 'm' is 
## mosquito biting rate, both of which could be correlated with
## infectivity?
gen_dust_det3 <- odin.dust::odin_dust({
  
  ## Defining first order differential equations
  deriv(S) <- - (beta * I/N) * S
  deriv(I) <- (beta * I/N) * S - gamma * I
  deriv(R) <- gamma * I
  output(n_cases) <- beta * I/N * S
  #output(rate_adj) <- rate_adj
  
  ## Model parameters - user defined when creating object from generator
  gamma <- user()
  m <- user()
  time[] <- user()
  rate[] <- user()
  #rate_adj[] <- rate[i] * m # rate is array, m is scalar
  rate_t <- interpolate(time, rate, "constant")
  beta <- rate_t * m
  
  ## Dimension of arrays
  dim(rate_adj) <- length(rate)
  dim(time) <- user()
  dim(rate) <- user()
  ## Initial conditions - user defined, defaults in parenthesis
  initial(S) <- 1000
  initial(I) <- 1
  initial(R) <- 0
  
  ## Total population size
  N = S + I + R
})

## Simulation happens over 50 days
time <- seq(1, 50, 1)

## Simulating some data, three different regimes (maybe rainy season vs not?)
## this is overly complicated/not smart but okay
rate <- sin(1:50) + 30 
rate[1:8] <- rate[1:8] + rexp(n=8, rate = 1)
rate[9:25] <- rate[9:25] + rexp(n=17, rate = 0.1) + 20
rate[26:50] <- rate[26:50] + rexp(n=25, rate = 1) + 10
rate <- rate / 100
m = 0.8
pars <- list(gamma = 0.2, m = m, time = time, rate = rate)

## Defining model - starting at time = 1 and only using 1 particle as it is 
## deterministic
mod4 <- gen_dust_det3$new(pars = pars,
                          time = 1,
                          n_particles = 1)
xForc <- array(data = NA, dim = c(mod4$info()$len, 1, length(time)))
for (t in 1:50) {xForc[ , , t] <- mod4$run(t)}

## 'mod2' defined earlier is simply the same model without forcing
mod2 <- gen_dust_det$new(pars= list(beta = 0.5, gamma = 0.2), # mod 2 is without forcing
                         time=1,
                         n_particles=n_particles)
xNoForc <- array(data = NA, dim = c(mod2$info()$len, 1, length(time)))
for (t in 1:50) {xNoForc[ , , t] <- mod2$run(t)}

## Plotting forcing vs no forcing
par(mfrow = c(1,2))
plot(xNoForc[4,,], type = "l", col = "blue", 
     xlab = "time", ylab = "Incidence")
lines(xForc[4,,], type = "l")
legend("topright", legend = c("No Forcing", "Forcing"), 
       col = c("black", "blue"), lty = 1, cex = 0.5)
plot(approx(time, rate*m, method = "linear"), type = "l",
     xlab = "time", ylab = "Beta = rate * m")