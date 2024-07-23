################################################################################
### --------------------- RUNNING MOISSALA INFERENCE ----------------------- ###
################################################################################
rm(list=ls())
library(dust)
library(odin.dust)
library(mcstate)
library(coda)
library(ggplot2)
library(reshape2) # for melt function
library(wesanderson)
library(patchwork)
library(coda)
library(GGally)
library(dplyr)
library(readxl)
library(forcats)

################################################################################
### -------------------- LOADING AUXILIARY FUNCTIONS ----------------------- ###
################################################################################
# process for saving the loaded data can be found in the script 
# 'savingMoissalaModelInputs.R' in this repository 
# https://github.com/SwissTPH/moissala-smc/tree/main/misc
data_dir <- "C:/Users/putnni/switchdrive/Chad/"
code_dir <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/"
dir_vcv <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/identifiability/vcv_saved/"
source(paste(code_dir, "aux_functions.R", sep = ""))
source(paste(code_dir, "misc/inf_run_stoch.R", sep = ""))

################################################################################
### --------------------- LOADING STOCHASTIC MODEL ------------------------- ###
################################################################################
path_full_model_SMC <- paste(code_dir, 
                             "models/full_model_SMC_stochastic.R", 
                             sep = "")
model_SMC <- odin.dust::odin_dust(path_full_model_SMC)

################################################################################
### ----------------- LOADING PROCESSED MODEL INPUTS ----------------------- ###
################################################################################
# process for saving the loaded data can be found in the script 
# 'savingMoissalaModelInputs.R' in this repository 
# https://github.com/SwissTPH/moissala-smc/tree/main/misc
load(paste(data_dir, "Data/model-inputs/moiss_inputs.Rdata", sep = ""))
source(paste(code_dir, "inputs-to-load/MCMCParameters.R", sep = ""))
source(paste(code_dir, "inputs-to-load/start_and_proposal_values.R", sep = ""))

################################################################################
### --------------------- DEFINING FIXED PARAMETERS ------------------------ ###
################################################################################
inputs_moiss <- read_excel(paste(data_dir, 
                                 "estimated-and-fixed-values/moiss_fixed.xlsx", 
                                 sep = ""))
inputs_moiss <- setNames(as.list(as.numeric(inputs_moiss$Value)), inputs_moiss$Name)
inputs_moiss[c("decay", "SMC", "cov_SMC", "c_R_D", "temp")] <- 
  c(list(decay_moiss), list(SMC_moiss), list(cov_moiss), list(anom_360_moiss), 
    list(temp_360_moiss))

################################################################################
### ----------------- DEFINING PARAMETERS TO ESTIMATE ---------------------- ###
################################################################################
params_to_estimate_moiss <- c(a_R = "a_R", b_R = "b_R", s = "s", 
                              qR = "qR", eff_SMC = "eff_SMC", 
                              size = "size",
                              z = "z", phi = "phi")


################################################################################
### ---------------- DEFINING PROPOSAL AND STARTING VALUES ----------------- ###
################################################################################
# defining starting values for chains
# these have to be the same order as what is defined in paramList in the inference functio
control_params_1$n_chains <- 2
start_values_moiss <- create_start_values(params_to_estimate_moiss, 
                                          control_params_1, paramList_name, 
                                          random = TRUE)

# choose proposal matrix
proposal_matrix_moiss <- create_proposal_matrix(proposal_variance, 
                                                params_to_estimate_moiss, 
                                                paramList_name)

################################################################################
### ---------------------- INFERENCE STAGE 1 ------------------------------- ###
################################################################################
# Dates of simulation, these are run a bit before observed data to 
# lower influence of starting conditions
dates <- c("2006-01-01", "2022-12-31") 

control_params_1$n_chains <- 3
control_params_1$n_steps <- 5000
start_values <- readRDS(paste(dir_vcv, "start_values_1.rds", sep = ""))
start_values <- t(apply(start_values, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
proposal_matrix <- readRDS(paste(dir_vcv, "vcv_1.rds", sep = ""))

# Running
results_1 <- inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                         control_params = control_params_1, 
                         params_to_estimate = params_to_estimate_moiss, 
                         proposal_matrix = proposal_matrix, 
                         start_values = start_values, noise = FALSE, seed = 24, 
                         month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                         synthetic = FALSE, incidence_df = cases_moiss, 
                         n_particles = 10, save_trajectories = FALSE, 
                         save_state = FALSE, rerun_n = 1000, rerun_random = TRUE, 
                         dir = code_dir)

# Saving new proposal and starting values for next stage
file_proposal_1 <- paste(dir_vcv, "vcv_1_stoch.rds", sep = "")
file_start_1 <- paste(dir_vcv, "start_values_1_stoch.rds", sep = "")
vcv_results_1 <- extract_vcv(results_2, S_prev = 20000,
                             file_proposal = file_proposal_1, file_start = file_start_1)

saveRDS(results_1, file = paste(data_dir, "MCMCRuns/results_moiss_stoch_1", sep = ""))

################################################################################
### ---------------------- INFERENCE STAGE 2 ------------------------------- ###
################################################################################
start_values_2 <- readRDS(paste(dir_vcv, "start_values_1_stoch.rds", sep = ""))
start_values_2 <- t(apply(start_values_2, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
proposal_matrix_2 <- readRDS(paste(dir_vcv, "vcv_1_stoch.rds", sep = ""))

# Running
control_params_1$n_steps <- 50000
control_params_1$n_burnin <- 10000
results_2 <- inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                           control_params = control_params_1, 
                           params_to_estimate = params_to_estimate_moiss, 
                           proposal_matrix = proposal_matrix_2, 
                           start_values = start_values_2, noise = FALSE, seed = 24, 
                           month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                           synthetic = FALSE, incidence_df = cases_moiss, 
                           n_particles = 10, save_trajectories = FALSE, 
                           save_state = FALSE, rerun_n = 1000, rerun_random = TRUE, 
                           dir = code_dir)

saveRDS(results_2, file = paste(data_dir, "MCMCRuns/results_moiss_stoch_2", sep = ""))
################################################################################
### ---------------------- SAVING ESTIMATED VALUES ------------------------- ###
################################################################################
post_quants_moiss <- as.data.frame(summary(results_2$coda_pars[,-c(1,2,3)])[2])
best_params <- results_2$coda_pars[which.max(results_2$coda_pars[,2]),-c(1,2,3)]
post_quants_moiss$best_params <- best_params
colnames(post_quants_moiss) <- c("2.5%", "25%", "50%", "75%", "97.5%", "best_params")
saveRDS(post_quants_moiss, "C:/Users/putnni/switchdrive/Chad/estimated-and-fixed-values/moiss_estimated_stoch.rds")

################################################################################
### ------------------ PARTICLE EFFECT ON PERFORMANCE ---------------------- ###
################################################################################
start.time <- Sys.time()
particle_1 <- inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                            control_params = control_params_1, 
                            params_to_estimate = params_to_estimate_moiss, 
                            proposal_matrix = proposal_matrix_2, 
                            start_values = start_values_2, noise = FALSE, seed = 24, 
                            month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                            synthetic = FALSE, incidence_df = cases_moiss, 
                            n_particles = 1, save_trajectories = TRUE, 
                            save_state = TRUE, rerun_n = 1000, rerun_random = TRUE, 
                            dir = code_dir)
end.time <- Sys.time()


particle_5 <- inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                               control_params = control_params_1, 
                               params_to_estimate = params_to_estimate_moiss, 
                               proposal_matrix = proposal_matrix_2, 
                               start_values = start_values_2, noise = FALSE, seed = 24, 
                               month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                               synthetic = FALSE, incidence_df = cases_moiss, 
                               n_particles = 5, save_trajectories = TRUE, 
                               save_state = TRUE, rerun_n = 1000, rerun_random = TRUE, 
                               dir = code_dir)

particle_50 <- inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                                control_params = control_params_1, 
                                params_to_estimate = params_to_estimate_moiss, 
                                proposal_matrix = proposal_matrix_2, 
                                start_values = start_values_2, noise = FALSE, seed = 24, 
                                month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                                synthetic = FALSE, incidence_df = cases_moiss, 
                                n_particles = 50, save_trajectories = TRUE, 
                                save_state = TRUE, rerun_n = 1000, rerun_random = TRUE, 
                                dir = code_dir)

# Define a function to wrap the call and measure time
run_and_time <- function(call) {
  start_time <- Sys.time()
  result <- call()
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  list(result = result, time = elapsed_time)
}

# Define the function calls
run_particle_1 <- function() {
  inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                control_params = control_params_1, 
                params_to_estimate = params_to_estimate_moiss, 
                proposal_matrix = proposal_matrix_2, 
                start_values = start_values_2, noise = FALSE, seed = 24, 
                month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                synthetic = FALSE, incidence_df = cases_moiss, 
                n_particles = 1, save_trajectories = TRUE, 
                save_state = TRUE, rerun_n = 1000, rerun_random = TRUE, 
                dir = code_dir)
}

run_particle_5 <- function() {
  inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                control_params = control_params_1, 
                params_to_estimate = params_to_estimate_moiss, 
                proposal_matrix = proposal_matrix_2, 
                start_values = start_values_2, noise = FALSE, seed = 24, 
                month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                synthetic = FALSE, incidence_df = cases_moiss, 
                n_particles = 5, save_trajectories = TRUE, 
                save_state = TRUE, rerun_n = 1000, rerun_random = TRUE, 
                dir = code_dir)
}

run_particle_10 <- function() {
  inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                control_params = control_params_1, 
                params_to_estimate = params_to_estimate_moiss, 
                proposal_matrix = proposal_matrix_2, 
                start_values = start_values_2, noise = FALSE, seed = 24, 
                month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                synthetic = FALSE, incidence_df = cases_moiss, 
                n_particles = 10, save_trajectories = TRUE, 
                save_state = TRUE, rerun_n = 1000, rerun_random = TRUE, 
                dir = code_dir)
}

run_particle_50 <- function() {
  inf_run_stoch(model_SMC = model_SMC, param_inputs = inputs_moiss, 
                control_params = control_params_1, 
                params_to_estimate = params_to_estimate_moiss, 
                proposal_matrix = proposal_matrix_2, 
                start_values = start_values_2, noise = FALSE, seed = 24, 
                month = TRUE, dates = dates, age_for_inf = "sep_ages", 
                synthetic = FALSE, incidence_df = cases_moiss, 
                n_particles = 50, save_trajectories = TRUE, 
                save_state = TRUE, rerun_n = 1000, rerun_random = TRUE, 
                dir = code_dir)
}

# Run the function calls and measure times
result_1 <- run_and_time(run_particle_1)
result_5 <- run_and_time(run_particle_5)
result_10 <- run_and_time(run_particle_10)
result_50 <- run_and_time(run_particle_50)

# Store the results in a list
results <- list(
  particle_1 = result_1$result,
  particle_1_time = result_1$time,
  particle_5 = result_5$result,
  particle_5_time = result_5$time,
  particle_10 = result_10$result,
  particle_10_time = result_10$time,
  particle_50 = result_50$result,
  particle_50_time = result_50$time
)

# Alternatively, store times in a data frame
times <- data.frame(
  Simulation = c("particle_1", "particle_5", "particle_10","particle_50"),
  Time = c(result_1$time, result_5$time, 
           result_10$time,result_50$time)
)

MCMC_diag(result_1$result)
MCMC_diag(result_5$result)
MCMC_diag(result_10$result)
MCMC_diag(result_50$result)
################################################################################
### ----------------------------- TESTING ---------------------------------- ###
################################################################################
start_values <- readRDS(paste(dir_vcv, "start_values_1.rds", sep = ""))
start_values <- t(apply(start_values, 1, function(x) runif(n = 2,  min = x * 0.8, max = x)))
proposal_matrix <- readRDS(paste(dir_vcv, "vcv_1.rds", sep = ""))

model_SMC <- model_SMC
param_inputs <- inputs_moiss
control_params <- control_params_1
params_to_estimate <- params_to_estimate_moiss
adaptive_param <- adaptive_param_1
#start_values <- start_values_1
month <- TRUE
dates <- dates
age_for_inf <- "sep_ages"
synthetic <- FALSE
incidence_df <- cases_moiss
clim_model = FALSE
save_trajectories = FALSE
rerun_n <- 2000
rerun_random = TRUE
dir <- code_dir
n_particles <- 5
