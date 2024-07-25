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
data_dir <- "C:/Users/putnni/switchdrive/Chad/"
code_dir <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/"
dir_vcv <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/identifiability/vcv_saved/"
source(paste(dir, "aux_functions.R", sep = ""))
load("C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/synthetic/inputs_synthetic")
source(paste(code_dir, "inputs-to-load/MCMCParameters.R", sep = ""))
source(paste(code_dir, "inputs-to-load/start_and_proposal_values.R", sep = ""))
source("C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/createResultsPlots.R")
dates <- c("2014-01-01", "2022-12-31")
source("C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/misc/posterior_plots_helper.R")

################################################################################
### -------------------- LOADING DETERMINISTIC MODEL ----------------------- ###
################################################################################
path_full_model_SMC <- paste(code_dir, 
                             "models/full_model_SMC_deterministic.R", 
                             sep = "")
model_SMC <- odin.dust::odin_dust(path_full_model_SMC)

################################################################################
### ----------------- INFERRING ONE PARAMETER AT A TIME -------------------- ###
################################################################################
params_to_estimate_full <- c(eff_SMC = "eff_SMC", a_R = "a_R", b_R = "b_R",
                             s = "s", qR = "qR", size = "size", phi = "phi", z = "z",
                             p_surv = "p_surv")

S_prev = 3000
# Defining MCMC control parameters
control_params_init <- list(n_steps = 3000, n_burnin = 0, n_chains = 3,
                            n_workers = 3, n_threads_total = 6)
control_params <- list(n_steps = 1000, n_burnin = 0, n_chains = 3,
                       n_workers = 3, n_threads_total = 6)

for (i in 1:length(params_to_estimate_full)) {
  params_to_estimate <- params_to_estimate_full[1:i]
  
  for (noise in c(FALSE, TRUE)) {
    noise_str <- ifelse(noise, "with_noise", "no_noise")
    file_suffix <- paste0(i, "_params_", noise_str)
    
    # Defining starting values for chains
    start_values <- create_start_values(params_to_estimate, 
                                        control_params, paramList_name, 
                                        random = TRUE)
    
    # Choose proposal matrix
    proposal_matrix_init <- create_proposal_matrix(proposal_variance, 
                                                   params_to_estimate, 
                                                   paramList_name)
    
    # Initial run
    results_init <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss_real_mo,
                            control_params = control_params_init, 
                            params_to_estimate = params_to_estimate, 
                            proposal_matrix = proposal_matrix_init, 
                            adaptive_param = NULL, noise = noise,
                            start_values = start_values, month = TRUE, 
                            dates = dates, age_for_inf = 'sep_ages',
                            synthetic = TRUE, incidence_df = NULL, 
                            save_trajectories = FALSE, 
                            rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
    
    proposal_matrix <- extract_vcv(results_init, S_prev = S_prev, save = FALSE)[[2]]
    
    # Main run
    results <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss_real_mo,
                       control_params = control_params, 
                       params_to_estimate = params_to_estimate, 
                       proposal_matrix = proposal_matrix, 
                       adaptive_param = NULL, noise = noise,
                       start_values = start_values, month = TRUE, 
                       dates = dates, age_for_inf = 'sep_ages',
                       synthetic = TRUE, incidence_df = NULL, 
                       save_trajectories = FALSE, 
                       rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
    
    # Analyze results
    analyze_results_synthetic(results, model_SMC, file_suffix = file_suffix, 
                              params_to_estimate, true_values = true_values)
  }
}

################################################################################
### -------------------- EXPLORING ADAPTIVE PARAMETERS --------------------- ###
################################################################################
# we use adaptive starting from inferring parameters where we have troubles (from qR)
# Proposed adaptive parameters to try
adaptive_param_init_1 <- adaptive_proposal_control()
adaptive_param_1 <- adaptive_proposal_control()

adaptive_param_init_2 <- adaptive_proposal_control(initial_vcv_weight = 1)
adaptive_param_2 <- adaptive_proposal_control(initial_vcv_weight = 100)

adaptive_param_init_3 <- adaptive_proposal_control(initial_vcv_weight = 1, 
                                                   initial_scaling = 2, 
                                                   initial_scaling_weight = NULL,
                                                   min_scaling = 0.1, 
                                                   scaling_increment = NULL, 
                                                   log_scaling_update = TRUE,
                                                   acceptance_target = 0.234, 
                                                   forget_rate = 0.6, 
                                                   forget_end = Inf,
                                                   adapt_end = Inf, 
                                                   pre_diminish = 10000)

adaptive_param_3 <- adaptive_proposal_control(
  initial_vcv_weight = 100, initial_scaling = 2, initial_scaling_weight = NULL,
  min_scaling = 0.1, scaling_increment = NULL, log_scaling_update = TRUE,
  acceptance_target = 0.234, forget_rate = 0.4, forget_end = Inf,
  adapt_end = Inf, pre_diminish = 5000)

# List of adaptive parameters for both initial and main runs
adaptive_params_init <- list(adaptive_param_init_1, adaptive_param_init_2, adaptive_param_init_3)
adaptive_params <- list(adaptive_param_1, adaptive_param_2, adaptive_param_3)

params_to_estimate <- c(eff_SMC = "eff_SMC", a_R = "a_R", b_R = "b_R",
                        s = "s", qR = "qR", size = "size", phi = "phi", z = "z")

noise <- TRUE  # Only run with noise

# Defining MCMC control parameters
control_params_init <- list(n_steps = 5000, n_burnin = 0, n_chains = 3,
                            n_workers = 3, n_threads_total = 6)
control_params <- list(n_steps = 20000, n_burnin = 0, n_chains = 3,
                       n_workers = 3, n_threads_total = 6)

for (i in seq_along(adaptive_params)) {
  adaptive_param_init <- adaptive_params_init[[i]]
  adaptive_param <- adaptive_params[[i]]
  file_suffix <- paste0("full_params_with_noise_adaptive_", i)
  
  # Defining starting values for chains
  start_values <- create_start_values(params_to_estimate, 
                                      control_params, paramList_name, 
                                      random = TRUE)
  
  # Choose proposal matrix
  proposal_matrix_init <- create_proposal_matrix(proposal_variance, 
                                                 params_to_estimate, 
                                                 paramList_name)
  
  # Initial run
  results_init <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss_real_mo,
                          control_params = control_params_init, 
                          params_to_estimate = params_to_estimate, 
                          proposal_matrix = proposal_matrix_init, 
                          adaptive_param = adaptive_param_init, noise = noise,
                          start_values = start_values, month = TRUE, 
                          dates = dates, age_for_inf = 'sep_ages',
                          synthetic = TRUE, incidence_df = NULL, 
                          save_trajectories = FALSE, 
                          rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
  
  proposal_matrix <- extract_vcv(results_init, S_prev = S_prev, save = FALSE)[[2]]
  start_values <- extract_vcv(results_init, S_prev = S_prev, save = FALSE)[[1]]
  start_values <- t(apply(start_values, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
  
  # Adjust start_values for adaptive_param_3
  # if (i == 3) {
  #   start_values <- extract_vcv(results_init, S_prev = S_prev, save = FALSE)[[1]]
  #   start_values <- t(apply(start_values, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
  # }
  # 
  # Main run
  results <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss_real_mo,
                     control_params = control_params, 
                     params_to_estimate = params_to_estimate, 
                     proposal_matrix = proposal_matrix, 
                     adaptive_param = adaptive_param, noise = noise,
                     start_values = start_values, month = TRUE, 
                     dates = dates, age_for_inf = 'sep_ages',
                     synthetic = TRUE, incidence_df = NULL, 
                     save_trajectories = FALSE, 
                     rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
  
  # Analyze results
  analyze_results_synthetic(results, model_SMC, file_suffix = file_suffix, 
                            params_to_estimate, true_values = true_values)
}

################################################################################
### ------------------ FIT TO SEP AGE GROUPS VS TOTAL  --------------------- ###
################################################################################
# List of age_for_inf values to try
age_for_inf_values <- c('sep_ages', 'total')

# Loop through each age_for_inf value
for (age_for_inf in age_for_inf_values) {
  file_suffix <- paste0("full_params_with_noise_", age_for_inf, "_adaptive_3")
  
  # Defining starting values for chains
  start_values <- create_start_values(params_to_estimate, 
                                      control_params, paramList_name, 
                                      random = TRUE)
  
  # Choose proposal matrix
  proposal_matrix_init <- create_proposal_matrix(proposal_variance, 
                                                 params_to_estimate, 
                                                 paramList_name)
  
  # Initial run
  results_init <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss_real_mo,
                          control_params = control_params_init, 
                          params_to_estimate = params_to_estimate, 
                          proposal_matrix = proposal_matrix_init, 
                          adaptive_param = adaptive_param_init_3, noise = noise,
                          start_values = start_values, month = TRUE, 
                          dates = dates, age_for_inf = age_for_inf,
                          synthetic = TRUE, incidence_df = NULL, 
                          save_trajectories = FALSE, 
                          rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
  
  #proposal_matrix <- extract_vcv(results_init, S_prev = S_prev, save = FALSE)[[2]]
  proposal_matrix <- extract_vcv(results_init, S_prev, save = FALSE)[[2]]
  # Adjust start values based on results_init
  start_values <- extract_vcv(results_init, S_prev, save = FALSE)[[1]]
  #start_values <- extract_vcv(results_init, S_prev = S_prev, save = FALSE)[[1]]
  start_values <- t(apply(start_values, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
  
  # Main run
  results <- inf_run(model_SMC = model_SMC, param_inputs = inputs_moiss_real_mo,
                     control_params = control_params, 
                     params_to_estimate = params_to_estimate, 
                     proposal_matrix = proposal_matrix, 
                     adaptive_param = adaptive_param_3, noise = noise,
                     start_values = start_values, month = TRUE, 
                     dates = dates, age_for_inf = age_for_inf,
                     synthetic = TRUE, incidence_df = NULL, 
                     save_trajectories = FALSE, 
                     rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
  
  # Analyze results
  analyze_results_synthetic(results, model_SMC, file_suffix = file_suffix, 
                            params_to_estimate, true_values = true_values)
}

################################################################################
### --------------------------- WEEKLY VS MONTHLY -------------------------- ###
################################################################################
# List of month values to try
month_values <- c(TRUE, FALSE)

# Loop through each month value
for (month in month_values) {
  month_str <- ifelse(month, "with_month", "no_month")
  
  # Assign the correct param_inputs based on the month value
  if (month) {
    param_inputs <- inputs_moiss_real_mo
  } else {
    param_inputs <- inputs_moiss_real_wk
  }
  
  file_suffix <- paste0("full_params_with_noise_", month_str, "_adaptive_3")
  
  # Defining starting values for chains
  start_values <- create_start_values(params_to_estimate, 
                                      control_params, paramList_name, 
                                      random = TRUE)
  
  # Choose proposal matrix
  proposal_matrix_init <- create_proposal_matrix(proposal_variance, 
                                                 params_to_estimate, 
                                                 paramList_name)
  
  # Initial run
  results_init <- inf_run(model_SMC = model_SMC, param_inputs = param_inputs,
                          control_params = control_params_init, 
                          params_to_estimate = params_to_estimate, 
                          proposal_matrix = proposal_matrix_init, 
                          adaptive_param = adaptive_param_init_3, noise = noise,
                          start_values = start_values, month = month, 
                          dates = dates, age_for_inf = 'sep_ages',
                          synthetic = TRUE, incidence_df = NULL, 
                          save_trajectories = FALSE, 
                          rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
  
  proposal_matrix <- extract_vcv(results_init, S_prev, save = FALSE)[[2]]
  
  # Adjust start values based on results_init
  start_values <- extract_vcv(results_init, S_prev, save = FALSE)[[1]]
  start_values <- t(apply(start_values, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
  
  # Main run
  results <- inf_run(model_SMC = model_SMC, param_inputs = param_inputs,
                     control_params = control_params, 
                     params_to_estimate = params_to_estimate, 
                     proposal_matrix = proposal_matrix, 
                     adaptive_param = adaptive_param_3, noise = noise,
                     start_values = start_values, month = month, 
                     dates = dates, age_for_inf = 'sep_ages',
                     synthetic = TRUE, incidence_df = NULL, 
                     save_trajectories = FALSE, 
                     rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
  
  # Analyze results
  analyze_results_synthetic(results, model_SMC, file_suffix = file_suffix, 
                            params_to_estimate, true_values = true_values)
}

################################################################################
### ------------------- BENEFIT OF NATURAL EXPERIMENT  --------------------- ###
################################################################################
# List of input sets and corresponding month values
inputs_list <- list(inputs_moiss_alt_mo = inputs_moiss_alt_mo, inputs_moiss_alt_wk = inputs_moiss_alt_wk, 
                    inputs_moiss_real_mo = inputs_moiss_real_mo, inputs_moiss_real_wk = inputs_moiss_real_wk)
month_values <- list(inputs_moiss_alt_mo = TRUE, inputs_moiss_alt_wk = FALSE, 
                     inputs_moiss_real_mo = TRUE, inputs_moiss_real_wk = FALSE)


# Loop through each input set, month value, and age_for_inf value
for (input_name in names(inputs_list)) {
  param_inputs <- inputs_list[[input_name]]
  month <- month_values[[input_name]]
  
  for (age_for_inf in c("sep_ages", "total")) {
    month_str <- ifelse(month, "with_month", "no_month")
    file_suffix <- paste0("full_params_with_noise_", month_str, "_", age_for_inf, "_adaptive_3_", input_name)
    
    # Defining starting values for chains
    start_values <- create_start_values(params_to_estimate, 
                                        control_params, paramList_name, 
                                        random = TRUE)
    
    # Choose proposal matrix
    proposal_matrix_init <- create_proposal_matrix(proposal_variance, 
                                                   params_to_estimate, 
                                                   paramList_name)
    
    # Initial run
    results_init <- inf_run(model_SMC = model_SMC, param_inputs = param_inputs,
                            control_params = control_params_init, 
                            params_to_estimate = params_to_estimate, 
                            proposal_matrix = proposal_matrix_init, 
                            adaptive_param = adaptive_param_init_3, noise = noise,
                            start_values = start_values, month = month, 
                            dates = dates, age_for_inf = age_for_inf,
                            synthetic = TRUE, incidence_df = NULL, 
                            save_trajectories = FALSE, 
                            rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
    
    proposal_matrix <- extract_vcv(results_init, S_prev, save = FALSE)[[2]]
    
    # Adjust start values based on results_init
    start_values <- extract_vcv(results_init, S_prev, save = FALSE)[[1]]
    start_values <- t(apply(start_values, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
    
    # Main run
    results <- inf_run(model_SMC = model_SMC, param_inputs = param_inputs,
                       control_params = control_params, 
                       params_to_estimate = params_to_estimate, 
                       proposal_matrix = proposal_matrix, 
                       adaptive_param = adaptive_param_3, noise = noise,
                       start_values = start_values, month = month, 
                       dates = dates, age_for_inf = age_for_inf,
                       synthetic = TRUE, incidence_df = NULL, 
                       save_trajectories = FALSE, 
                       rerun_n = 1000, rerun_random = TRUE, dir = code_dir)
    
    # Analyze results
    analyze_results_synthetic(results, model_SMC, file_suffix = file_suffix, 
                              params_to_estimate, true_values = true_values)
  }
}
################################################################################
start_values <- create_start_values(params_to_estimate, control_params, paramList_name, random = TRUE)

dir_save <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/synthetic/test_plots/"
# Save results to a specified directory
save_results <- function(results, file_suffix, dir_save) {
  save_path <- file.path(dir_save, paste0("results_", file_suffix, ".rds"), fsep = "")
  saveRDS(results, save_path)
}

# Run simulations and save results
run_simulation <- function(model_SMC, param_inputs, control_params_init, control_params, params_to_estimate, proposal_variance, adaptive_param_init, adaptive_param, noise, start_values, month, dates, age_for_inf, dir, file_suffix, dir_save, S_prev) {
  save_path <- file.path(dir_save, paste0("results_", file_suffix, ".rds"), fsep = "")
  if (file.exists(save_path)) {
    cat("Results file already exists for", file_suffix, ". Skipping...\n")
    return()
  }
  
  # Initial run
  results_init <- inf_run(
    model_SMC = model_SMC, 
    param_inputs = param_inputs,
    control_params = control_params_init, 
    params_to_estimate = params_to_estimate, 
    proposal_matrix = create_proposal_matrix(proposal_variance, params_to_estimate, paramList_name), 
    adaptive_param = adaptive_param_init, 
    noise = noise,
    start_values = start_values, 
    month = month, 
    dates = dates, 
    age_for_inf = age_for_inf,
    synthetic = TRUE, 
    incidence_df = NULL, 
    save_trajectories = FALSE, 
    rerun_n = 1000, 
    rerun_random = TRUE, 
    dir = dir
  )
  
  # Extract proposal matrix and start values for main run
  proposal_matrix <- extract_vcv(results_init, S_prev = S_prev, save = FALSE)[[2]]
  start_values <- extract_vcv(results_init, S_prev = S_prev, save = FALSE)[[1]]
  start_values[start_values <= 0] <- 0.01
  start_values <- t(apply(start_values, 1, function(x) runif(n = 3,  min = x * 0.8, max = x)))
  
  # Main run
  results <- inf_run(
    model_SMC = model_SMC, 
    param_inputs = param_inputs,
    control_params = control_params, 
    params_to_estimate = params_to_estimate, 
    proposal_matrix = proposal_matrix, 
    adaptive_param = adaptive_param, 
    noise = noise,
    start_values = start_values, 
    month = month, 
    dates = dates, 
    age_for_inf = age_for_inf,
    synthetic = TRUE, 
    incidence_df = NULL, 
    save_trajectories = FALSE, 
    rerun_n = 1000, 
    rerun_random = TRUE, 
    dir = dir
  )
  
  # Save the results
  save_results(results, file_suffix, dir_save)
}
# Load and analyze results
load_and_analyze_results <- function(file_suffix, model_SMC, params_to_estimate, true_values, dir_save) {
  load_path <- file.path(dir_save, paste0("results_", file_suffix, ".rds"), fsep = "")
  #load_path <- file.path(dir_save, "saved_results", paste0("results_", file_suffix, ".rds"), sep = "")
  results <- readRDS(load_path)
  analyze_results_synthetic(results, model_SMC, file_suffix, params_to_estimate, true_values)
}

# Define parameters for the simulation
params_to_estimate_full <- c(eff_SMC = "eff_SMC", a_R = "a_R", b_R = "b_R", s = "s", qR = "qR", size = "size", phi = "phi", z = "z", p_surv = "p_surv")
params_to_estimate <- c(eff_SMC = "eff_SMC", a_R = "a_R", b_R = "b_R",
                        s = "s", qR = "qR", size = "size", phi = "phi", z = "z")
S_prev <- 100
control_params_init <- list(n_steps = 5000, n_burnin = 0, n_chains = 3, n_workers = 3, n_threads_total = 6)
control_params <- list(n_steps = 40000, n_burnin = 0, n_chains = 3, n_workers = 3, n_threads_total = 6)

### INFERRING ONE PARAMETER AT A TIME ###
one_param_file_suffixes <- c()
for (i in 1:length(params_to_estimate_full)) {
  params_to_estimate <- params_to_estimate_full[1:i]
  for (noise in c(FALSE, TRUE)) {
    noise_str <- ifelse(noise, "with_noise", "no_noise")
    file_suffix <- paste0(i, "_params_", noise_str)
    one_param_file_suffixes <- c(one_param_file_suffixes, file_suffix)
    
    # Run simulation and save results
    run_simulation(model_SMC, inputs_moiss_real_mo, control_params_init, control_params, params_to_estimate, proposal_variance, NULL, NULL, noise, start_values, TRUE, dates, 'sep_ages', dir, file_suffix, dir_save = dir_save, S_prev = S_prev)
  }
}

### EXPLORING ADAPTIVE PARAMETERS ###
adaptive_param_init_1 <- adaptive_proposal_control()
adaptive_param_1 <- adaptive_proposal_control()

adaptive_param_init_2 <- adaptive_proposal_control(initial_vcv_weight = 1)
adaptive_param_2 <- adaptive_proposal_control(initial_vcv_weight = 100)

adaptive_param_init_3 <- adaptive_proposal_control(initial_vcv_weight = 1, initial_scaling = 2, initial_scaling_weight = NULL, min_scaling = 0.1, scaling_increment = NULL, log_scaling_update = TRUE, acceptance_target = 0.234, forget_rate = 0.6, forget_end = Inf, adapt_end = Inf, pre_diminish = 10000)
adaptive_param_3 <- adaptive_proposal_control(initial_vcv_weight = 100, initial_scaling = 2, initial_scaling_weight = NULL, min_scaling = 0.1, scaling_increment = NULL, log_scaling_update = TRUE, acceptance_target = 0.234, forget_rate = 0.4, forget_end = Inf, adapt_end = Inf, pre_diminish = 5000)

adaptive_params_init <- list(adaptive_param_init_1, adaptive_param_init_2, adaptive_param_init_3)
adaptive_params <- list(adaptive_param_1, adaptive_param_2, adaptive_param_3)

adaptive_params_file_suffixes <- c()
for (i in seq_along(adaptive_params)) {
  adaptive_param_init <- adaptive_params_init[[i]]
  adaptive_param <- adaptive_params[[i]]
  file_suffix <- paste0("full_params_with_noise_adaptive_", i)
  adaptive_params_file_suffixes <- c(adaptive_params_file_suffixes, file_suffix)
  
  start_values <- create_start_values(params_to_estimate, control_params, paramList_name, random = TRUE)
  
  run_simulation(model_SMC, inputs_moiss_real_mo, control_params_init, control_params, params_to_estimate, proposal_variance, adaptive_param_init, adaptive_param, TRUE, start_values, TRUE, dates, 'sep_ages', dir, file_suffix, dir_save = dir_save, S_prev = S_prev)
}

### FIT TO SEP AGE GROUPS VS TOTAL ###
age_for_inf_values <- c('sep_ages', 'total')
age_groups_file_suffixes <- c()
for (age_for_inf in age_for_inf_values) {
  file_suffix <- paste0("full_params_with_noise_", age_for_inf)
  age_groups_file_suffixes <- c(age_groups_file_suffixes, file_suffix)
  
  start_values <- create_start_values(params_to_estimate, control_params, paramList_name, random = TRUE)
  
  run_simulation(model_SMC, inputs_moiss_real_mo, control_params_init, control_params, params_to_estimate, proposal_variance, adaptive_param_init_3, adaptive_param_3, TRUE, start_values, TRUE, dates, age_for_inf, dir, file_suffix, dir_save = dir_save, S_prev = S_prev)
}

### WEEKLY VS MONTHLY ###
month_values <- c(TRUE, FALSE)
monthly_vs_weekly_file_suffixes <- c()
for (month in month_values) {
  month_str <- ifelse(month, "with_month", "no_month")
  month_str <- ifelse(month, "with_month", "no_month")
  # Assign the correct param_inputs based on the month value
  if (month) {
    param_inputs <- inputs_moiss_real_mo
  } else {
    param_inputs <- inputs_moiss_real_wk
  }
  file_suffix <- paste0("full_params_with_noise_", month_str)
  monthly_vs_weekly_file_suffixes <- c(monthly_vs_weekly_file_suffixes, file_suffix)
  
  start_values <- create_start_values(params_to_estimate, control_params, paramList_name, random = TRUE)
  
  run_simulation(model_SMC, param_inputs, control_params_init, control_params, params_to_estimate, proposal_variance, adaptive_param_init_3, adaptive_param_3, TRUE, start_values, month, dates, 'sep_ages', dir, file_suffix, dir_save = dir_save, S_prev = S_prev)
}

### BENEFIT OF NATURAL EXPERIMENT ###
inputs_list <- list(inputs_moiss_alt_mo = inputs_moiss_alt_mo, inputs_moiss_alt_wk = inputs_moiss_alt_wk, inputs_moiss_real_mo = inputs_moiss_real_mo, inputs_moiss_real_wk = inputs_moiss_real_wk)
month_values <- list(inputs_moiss_alt_mo = TRUE, inputs_moiss_alt_wk = FALSE, inputs_moiss_real_mo = TRUE, inputs_moiss_real_wk = FALSE)

natural_experiment_file_suffixes <- c()
for (input_name in names(inputs_list)) {
  param_inputs <- inputs_list[[input_name]]
  month <- month_values[[input_name]]
  
  for (age_for_inf in c("sep_ages", "total")) {
    month_str <- ifelse(month, "with_month", "no_month")
    file_suffix <- paste0("full_params_with_noise_", month_str, "_", age_for_inf, input_name)
    natural_experiment_file_suffixes <- c(natural_experiment_file_suffixes, file_suffix)
    
    start_values <- create_start_values(params_to_estimate, control_params, paramList_name, random = TRUE)
    
    run_simulation(model_SMC, param_inputs, control_params_init, control_params, params_to_estimate, proposal_variance, adaptive_param_init_3, adaptive_param_3, TRUE, start_values, month, dates, age_for_inf, dir, file_suffix, dir_save = dir_save, S_prev = S_prev)
  }
}

file_suffixes <- list(
  one_param_file_suffixes,
  adaptive_params_file_suffixes,
  age_groups_file_suffixes,
  monthly_vs_weekly_file_suffixes,
  natural_experiment_file_suffixes
)

K <- ceiling(seq_along(file_suffixes[[1]])/2) # repeating same params for noise + no noise
for (i in seq_along(file_suffixes[[1]])) {
  params_to_estimate <- params_to_estimate_full[1:K[i]]
  load_and_analyze_results(file_suffixes[[1]][i], model_SMC, params_to_estimate, true_values, dir_save)
}

for (section_suffixes in file_suffixes[-1]) {
  for (file_suffix in section_suffixes) {
    load_and_analyze_results(file_suffix, model_SMC, params_to_estimate, true_values, dir_save)
  }
}

################################################################################
model_SMC = model_SMC
param_inputs = inputs_moiss_real
control_params = control_params
params_to_estimate = params_to_estimate_1
proposal_matrix = proposal_matrix_1
adaptive_param = NULL
start_values = start_values_1
month = TRUE
dates = dates
age_for_inf = 'sep_ages'
synthetic = TRUE
incidence_df = NULL
save_trajectories = FALSE
rerun_n = Inf
rerun_random = FALSE
dir = code_dir
noise = TRUE
seed = 24
