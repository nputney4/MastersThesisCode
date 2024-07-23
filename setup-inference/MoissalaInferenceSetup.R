################################################################################
### --------------------- MOISSALA INFERENCE SETUP ------------------------- ###
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
load(paste(data_dir, "Data/model-inputs/moiss_inputs_10.Rdata", sep = ""))
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
                              size = "size", p_surv = "p_surv",
                              z = "z", phi = "phi")

################################################################################
### ---------------- DEFINING PROPOSAL AND STARTING VALUES ----------------- ###
################################################################################
# Dates of simulation, these are run a bit before observed data to 
# lower influence of starting conditions
dates <- c("2006-01-01", "2022-12-31") 

# defining starting values for chains
# these have to be the same order as what is defined in paramList in the inference function
start_values_moiss <- create_start_values(params_to_estimate_moiss, 
                                          control_params_1, paramList_name, 
                                          random = TRUE)

# choose proposal matrix
proposal_matrix_moiss <- create_proposal_matrix(proposal_variance, 
                                                params_to_estimate_moiss, 
                                                paramList_name)