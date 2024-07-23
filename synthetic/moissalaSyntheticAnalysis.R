################################################################################
### ------------------------ SYNTHETIC ANALYSIS ---------------------------- ###
################################################################################

################################################################################
### ------------------- LOADING FUNCTIONS AND INPUTS ----------------------- ###
################################################################################
dir <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/"
dir_data <- "C:/Users/putnni/switchdrive/Chad/Data/cases-data/"
dir_clim <-  "C:/Users/putnni/switchdrive/Chad/Data/climate-data/"
data_dir <- "C:/Users/putnni/switchdrive/Chad/"
source(paste(dir, "aux_functions.R", sep = ""))
load("C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/synthetic/inputs_synthetic")

################################################################################
### ----------------- INFERRING ONE PARAMETER AT A TIME -------------------- ###
################################################################################
params_to_estimate_1 <- c(eff_SMC = "eff_SMC")

# defining MCMC control parameters
control_params <- list(n_steps = 20000, n_burnin = 0, n_chains = 3,
                       n_workers = 3, n_threads_total = 6)

# defining starting values for chains
# these have to be the same order as what is defined in paramList in the inference function
start_values_moiss <- create_start_values(params_to_estimate_moiss, 
                                          control_params_1, paramList_name, 
                                          random = TRUE)

# choose proposal matrix
proposal_matrix_moiss <- create_proposal_matrix(proposal_variance, 
                                                params_to_estimate_moiss, 
                                                paramList_name)

results_moiss_2 <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss_real,
                           control_params = control_params, 
                           params_to_estimate = params_to_estimate_1, 
                           proposal_matrix = proposal_matrix_2, 
                           adaptive_param = NULL, 
                           start_values = start_values_2, month = TRUE, 
                           dates = dates, age_for_inf = 'sep_ages',
                           synthetic = FALSE, incidence_df = cases_moiss, 
                           save_trajectories = FALSE, 
                           rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
