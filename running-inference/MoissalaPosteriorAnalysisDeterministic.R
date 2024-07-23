################################################################################
### ------------------------ POSTERIOR ANALYSIS ---------------------------- ###
################################################################################
library(flextable)
library(dplyr)
data_dir <- "C:/Users/putnni/switchdrive/Chad/"
# results_moiss <- readRDS(file = paste(data_dir, "MCMCRuns/results_moiss_070224", sep = ""))
# results_moiss <- readRDS(file = paste(data_dir, "MCMCRuns/results_moiss_070724", sep = ""))
results_moiss <- readRDS(file = paste(data_dir, "MCMCRuns/results_moiss_det_2", sep = ""))
params_to_estimate_moiss <- c(a_R = "a_R", b_R = "b_R", s = "s", 
                              qR = "qR", z = "z", phi = "phi", eff_SMC = "eff_SMC")
params_to_estimate_moiss <- c(a_R = "a_R", b_R = "b_R", s = "s", 
                              qR = "qR", eff_SMC = "eff_SMC", 
                              size = "size", phi = "phi", z = "z")

theme_plots <- theme(plot.title = element_text(hjust=0.5, size=20),
                     axis.text.x = element_text(size = 14),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size = 16),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size = 14)) 

################################################################################
### -------------------- LOADING DETERMINISTIC MODEL ----------------------- ###
################################################################################
path_full_model_SMC <- paste(code_dir, 
                             "models/full_model_SMC_deterministic.R", 
                             sep = "")
model_SMC <- odin.dust::odin_dust(path_full_model_SMC)

################################################################################
### --------------------- FUNCTION SAVING ALL PLOTS ------------------------ ###
################################################################################
source("C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/createResultsPlots.R")
analyze_results(results_moiss, model_SMC, stochastic = FALSE)

################################################################################
### -------------------- LOADING AUXILIARY FUNCTIONS ----------------------- ###
################################################################################
# process for saving the loaded data can be found in the script 
code_dir <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/"
plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"
source(paste(code_dir, "aux_functions.R", sep = ""))

################################################################################
### ----------------- IMPROVING READABILITY OF PLOTS ----------------------- ###
################################################################################
theme_plots <- theme(plot.title = element_text(hjust=0.5, size=20),
                     axis.text.x = element_text(size = 14),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size = 16),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size = 14)) 

################################################################################
### --------------------- INITIAL CHECK OF MCMC ---------------------------- ###
################################################################################
MCMC_diag(results_moiss) # trace plots, effective sample size
plot_corr(results_moiss) # examining correlation of draws

################################################################################
### -------------------- MARGINAL POSTERIOR DISTRIBUTIONS ------------------ ###
################################################################################
source(paste(code_dir, "misc/posterior_plots_helper.R", sep = ""))
dim_plot <- c(5,2)
p1 <- post_plot(results_moiss, params_to_estimate_moiss, dim_plot = dim_plot,
          show_true = FALSE)

ggsave(p1, filename = paste(plot_loc, "results/moiss_post_distributions.pdf", sep = ""), device = "pdf",
       height = 20, width = 6, units = "in")

################################################################################
### -------------------- POSTERIOR PREDICTIVE CHECK ------------------------ ###
################################################################################
dates_sim <- c("2006-01-01", "2022-12-31") # dates to be used for simulated data
dates_obs <- c("2014-01-01", "2022-12-31") # dates corresponding ot observed data
sim_df <- create_sim_df(results_moiss, 50, dates_sim, dates_obs, model_SMC)
clim_df <- create_clim_df(moiss_met_2014_2023)
p2 <- post_pred_plot(results_moiss, sim_df, dates_obs, dates_obs, show_clim = TRUE, clim_df = clim_df,
               title = "Moissala", theme = theme_plots, ages = c("o5", "u5", "total"))

ggsave(p2, filename = paste(plot_loc, "results/moiss_post_predictive_check.pdf", sep = ""), device = "pdf",
       height = 6, width = 9, units = "in")

################################################################################
### -------------------- CALCULATING CASES AVERTED ------------------------- ###
################################################################################
# Simulate model with SMC Active
estimated_values <- readRDS("C:/Users/putnni/switchdrive/Chad/estimated-and-fixed-values/moiss_estimated.rds")
inputs_moiss[rownames(estimated_values)] <- estimated_values$`50%` # fixing values to those from estimated for inference procedure
sim_SMC <- data_sim(model_SMC, inputs_moiss,
                    start_date = dates_sim[1], end_date = dates_sim[2],
                    month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)
sim_SMC <- sim_SMC[which(sim_SMC$date_ymd %in% seq(as.Date(dates_obs[1]), as.Date(dates_obs[2]), by = "day")), ]

plot(cases_moiss$month, cases_moiss$inc_C)
lines(sim_SMC$date_ymd, sim_SMC$inc_C, type = "l")

# Simulate model without SMC
inputs_moiss_2 <- inputs_moiss
inputs_moiss_2$eff_SMC <- 0
sim_noSMC <- data_sim(model_SMC, inputs_moiss_2,
                      start_date = dates_sim[1], end_date = dates_sim[2],
                      month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)
sim_noSMC <- sim_noSMC[which(sim_noSMC$date_ymd %in% seq(as.Date(dates_obs[1]), as.Date(dates_obs[2]), by = "day")), ]

# Visualize cases averted
p3 <- cases_averted_plot(sim_SMC, sim_noSMC, show_observed = FALSE)
ggsave(p3, filename = paste(plot_loc, "results/case_averted_plot.pdf", sep = ""), device = "pdf",
       height = 6, width = 11, units = "in")

# Calculate cases averted
total_noSMC <- sim_noSMC$inc %>% sum()
u5_noSMC <- sim_noSMC$inc_C %>% sum()
total_SMC <- sim_SMC$inc %>% sum()
u5_SMC <- sim_SMC$inc_C %>% sum()
total_obs <- cases_moiss$inc
u5_obs <- cases_moiss$inc_C

# Create table from case averted calculations
total_CA <- data.frame("Age Group" = "total",
                       "Number of Cases" = round(total_SMC) , 
                       "Cases Averted" = round(total_noSMC - total_SMC), 
                       "Percent Reduction" = 100 * round((total_noSMC - total_SMC) / total_SMC, 4))


u5_CA <- data.frame("Age Group" = "under 5",
                    "Number of Cases" = round(u5_SMC), 
                    "Cases Averted" = round(u5_noSMC - u5_SMC), 
                    "Percent Reduction" = 100 * round((u5_noSMC - u5_SMC) / u5_SMC, 4))

comb_CA <- rbind(u5_CA, total_CA)

table_flextable <- flextable(comb_CA) %>%
  colformat_num(j = c("Number.of.Cases", "Cases.Averted"), big.mark = ",", digits = 0) %>%
  colformat_num(j = "Percent.Reduction", digits = 2) %>%
  set_header_labels(
    `Number of Cases` = "Number of Cases",
    `Cases Averted` = "Cases Averted",
    `Percent Reduction` = "% Reduction"
  ) %>%
  theme_box() %>%
  autofit() %>%
  set_table_properties(layout = "autofit") %>%
  add_header_row(values = c("", "Cases Averted in Moissala"), colwidths = c(1, 3))

# Save table as .svg
save_as_image(table_flextable, 
              path = paste(plot_loc, "results/cases_averted_table.svg", sep = ""))


################################################################################
### ------------------- CHANGING DEPLOYMENT STRATEGY ----------------------- ###
################################################################################
start_date <- ymd("2014-01-01")
end_date <- ymd("2022-12-31")
years <- c(2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)
times_per_year <- c(4, 4, 4, 4, 4, 4, 4, 4, 4)
months_active <- matrix(c(0,0,0,0,0,0,1,1,1,1,0,0,
                          0,0,0,0,0,0,1,1,1,1,0,0,
                          0,0,0,0,0,0,1,1,1,1,0,0,
                          0,0,0,0,0,0,1,1,1,1,0,0,
                          0,0,0,0,0,0,1,1,1,1,0,0,
                          0,0,0,0,0,0,1,1,1,1,0,0,
                          0,0,0,0,0,0,1,1,1,1,0,0,
                          0,0,0,0,0,0,1,1,1,1,0,0,
                          0,0,0,0,0,0,1,1,1,1,0,0), nrow = length(years), byrow = TRUE)

# Generate SMC schedule dataframe
schedule_df_1 <- generate_smc_schedule_df_synthetic(start_date, end_date, years, times_per_year, months_active,
                                          months_30_days = FALSE)

inputs_moiss_3 <- inputs_moiss
sim_SMC2 <- data_sim(model_SMC, inputs_moiss_3,
                     start_date = dates_sim[1], 
                     end_date = dates_sim[2],
                                  month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)