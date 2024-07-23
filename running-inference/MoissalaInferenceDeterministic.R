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
r_D_df <- readRDS(paste(data_dir, "Data/model-inputs/daily_growth_rates.rds", sep = ""))
r_C <- 1 + r_D_df$u5
r_A <- 1 + r_D_df$o5

################################################################################
### -------------------- LOADING DETERMINISTIC MODEL ----------------------- ###
################################################################################
path_full_model_SMC <- paste(code_dir, 
                             "models/full_model_SMC_deterministic.R", 
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

#inputs_moiss$r_C = r_C
#inputs_moiss$r_A = r_A
################################################################################
### ----------------- DEFINING PARAMETERS TO ESTIMATE ---------------------- ###
################################################################################
params_to_estimate_moiss <- c(a_R = "a_R", b_R = "b_R", s = "s", 
                              qR = "qR", z = "z", eff_SMC = "eff_SMC", "phi",
                              size = "size", p_surv = "p_surv")

params_to_estimate_moiss <- c(a_R = "a_R", b_R = "b_R", s = "s", 
                              qR = "qR", eff_SMC = "eff_SMC", 
                              size = "size",
                              z = "z", phi = "phi")

################################################################################
### ---------------- DEFINING PROPOSAL AND STARTING VALUES ----------------- ###
################################################################################
# defining starting values for chains
# these have to be the same order as what is defined in paramList in the inference function
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

# defining starting values for chains
# these have to be the same order as what is defined in paramList in the inference function
start_values_1 <- create_start_values(params_to_estimate_moiss, 
                                          control_params_1, paramList_name, 
                                          random = TRUE)

# choose proposal matrix
proposal_matrix_1 <- create_proposal_matrix(proposal_variance, 
                                                params_to_estimate_moiss, 
                                                paramList_name)

#control_params_1$n_steps <- 5000
results_moiss_1 <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss,
                           control_params = control_params_1, 
                           params_to_estimate = params_to_estimate_moiss, 
                           proposal_matrix = proposal_matrix_1, 
                           adaptive_param = adaptive_param_1, 
                           start_values = start_values_1, month = TRUE, 
                           dates = dates, age_for_inf = 'sep_ages',
                           synthetic = FALSE, incidence_df = cases_moiss, 
                           save_trajectories = FALSE, 
                           rerun_n = 1000, rerun_random = TRUE, dir = code_dir)

saveRDS(results_moiss_1, file = paste(data_dir, "MCMCRuns/results_moiss_det_1", sep = ""))
#MCMC_diag(results_moiss_1)
# plot_corr(results_moiss_1)
results_moiss_1$coda_pars[,2] %>% max() # -1712

 # Saving new proposal and starting values for next stage
file_proposal_1 <- paste(dir_vcv, "vcv_1.rds", sep = "")
file_start_1 <- paste(dir_vcv, "start_values_1.rds", sep = "")
vcv_results_1 <- extract_vcv(results_moiss_1, S_prev = 3000,
                                    file_proposal = file_proposal_1, file_start = file_start_1)

################################################################################
### ---------------------- INFERENCE STAGE 2 ------------------------------- ###
################################################################################
start_values_2 <- readRDS(paste(dir_vcv, "start_values_1.rds", sep = ""))
start_values_2 <- t(apply(start_values_2, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
#start_values_2[4,] <- rep(-1.54, 4)
proposal_matrix_2 <- readRDS(paste(dir_vcv, "vcv_1.rds", sep = ""))

#control_params_2$n_steps <- 15000
#control_params_2$n_burnin <- 5000

results_moiss_2 <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss,
                           control_params = control_params_2, 
                           params_to_estimate = params_to_estimate_moiss, 
                           proposal_matrix = proposal_matrix_2, 
                           adaptive_param = adaptive_param_2, 
                           start_values = start_values_2, month = TRUE, 
                           dates = dates, age_for_inf = 'sep_ages',
                           synthetic = FALSE, incidence_df = cases_moiss, 
                           save_trajectories = FALSE, 
                           rerun_n = 1000, rerun_random = TRUE, dir = code_dir)

#MCMC_diag(results_moiss_2)
#plot_corr(results_moiss_2)
results_moiss_2$coda_pars[,2] %>% max() # -1712
saveRDS(results_moiss_2, file = paste(data_dir, "MCMCRuns/results_moiss_det_2", sep = ""))

# Saving new proposal and starting values for next stage
file_proposal_2 <- paste(dir_vcv, "vcv_2.rds", sep = "")
file_start_2 <- paste(dir_vcv, "start_values_2.rds", sep = "")
vcv_results_2 <- extract_vcv(results_moiss_2, S_prev = 3000,
                             file_proposal = file_proposal_2, file_start = file_start_2)

# ################################################################################
# ### ---------------------- INFERENCE STAGE 3 ------------------------------- ###
# ################################################################################
# # Dates of simulation, these are run a bit before observed data to 
# # lower influence of starting conditions
# dates <- c("2006-01-01", "2022-12-31") 
# 
# # Loading in start values and proposal variance based on previous runs
# # This will start the MCMC run closer to the posterior mode and with
# # a covariance matrix closer to that of the true posterior
# start_values_3 <- readRDS(paste(dir_vcv, "start_values_2.rds", sep = ""))
# proposal_matrix_3 <- readRDS(paste(dir_vcv, "vcv_2.rds", sep = ""))
# 
# results_moiss_3 <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss,
#                             control_params = control_params_3, 
#                             params_to_estimate = params_to_estimate_moiss, 
#                             proposal_matrix = proposal_matrix_3, 
#                             adaptive_param = adaptive_param_3, 
#                             start_values = start_values_3, month = TRUE, 
#                             dates = dates, age_for_inf = 'sep_ages',
#                             synthetic = FALSE, incidence_df = cases_moiss, 
#                             save_trajectories = FALSE, 
#                             rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
# 
# saveRDS(results_moiss_3, file = paste(data_dir, "MCMCRuns/results_moiss_070224", sep = ""))
# MCMC_diag(results_moiss_3)
# results_moiss_3$coda_pars[,2] %>% max() # -1712
################################################################################
### ---------------------- SAVING ESTIMATED VALUES ------------------------- ###
################################################################################
post_quants_moiss <- as.data.frame(summary(results_moiss_2$coda_pars[,-c(1,2,3)])[2])
best_params <- results_moiss_2$coda_pars[which.max(results_moiss_2$coda_pars[,2]),-c(1,2,3)]
post_quants_moiss$best_params <- best_params
colnames(post_quants_moiss) <- c("2.5%", "25%", "50%", "75%", "97.5%", "best_params")
saveRDS(post_quants_moiss, "C:/Users/putnni/switchdrive/Chad/estimated-and-fixed-values/moiss_estimated_det.rds")
