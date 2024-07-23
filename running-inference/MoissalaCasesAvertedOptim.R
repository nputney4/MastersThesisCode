################################################################################
### ------------------------------ SETUP ----------------------------------- ###
################################################################################
library(flextable)
library(lubridate)
library(xtable)
source("C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/setup-inference/MoissalaInferenceSetup.R")
data_dir <- "C:/Users/putnni/switchdrive/Chad/"
dates_sim <- c("2006-01-01", "2022-12-31") # dates to be used for simulated data
dates_obs <- c("2014-01-01", "2022-12-31") # dates corresponding ot observed data
plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"
theme_plots <- theme(plot.title = element_blank(),
                     axis.text.x = element_text(size = 14),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size = 16),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size = 14))
# Simulate model with SMC Active
estimated_values <- readRDS("C:/Users/putnni/switchdrive/Chad/estimated-and-fixed-values/moiss_estimated.rds")

################################################################################
### -------------------- CALCULATING CASES AVERTED ------------------------- ###
################################################################################
inputs_moiss[rownames(estimated_values)] <- estimated_values$best_params # fixing values to those from estimated for inference procedure
sim_SMC <- data_sim(model_SMC, inputs_moiss,
                    start_date = dates_sim[1], end_date = dates_sim[2],
                    month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)
sim_SMC <- sim_SMC[which(sim_SMC$date_ymd %in% seq(as.Date(dates_obs[1]), as.Date(dates_obs[2]), by = "day")), ]

# Simulate model without SMC
inputs_moiss_2 <- inputs_moiss
inputs_moiss_2$eff_SMC <- 0
sim_noSMC <- data_sim(model_SMC, inputs_moiss_2,
                      start_date = dates_sim[1], end_date = dates_sim[2],
                      month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)
sim_noSMC <- sim_noSMC[which(sim_noSMC$date_ymd %in% seq(as.Date(dates_obs[1]), as.Date(dates_obs[2]), by = "day")), ]

# Visualize cases averted
p3 <- cases_averted_plot(sim_SMC, sim_noSMC, show_observed = FALSE, theme = theme_plots)
ggsave(p3, filename = paste(plot_loc, "results/deterministic/optim/case_averted_plot_det.pdf", sep = ""), device = "pdf",
       height = 6, width = 11, units = "in")

# Calculate cases averted
months_2019 <- ((2019 - 2014)*12) : ((2020 - 2014) * 12)

total_noSMC <- sim_noSMC$inc[-months_2019] %>% sum()
u5_noSMC <- sim_noSMC$inc_C[-months_2019] %>% sum()
total_SMC <- sim_SMC$inc[-months_2019] %>% sum()
u5_SMC <- sim_SMC$inc_C[-months_2019] %>% sum()
total_obs <- cases_moiss$inc[-months_2019]
u5_obs <- cases_moiss$inc_C[-months_2019]

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
comb_CA_tex <- xtable(comb_CA)
print(comb_CA_tex, file=paste(plot_loc, "results/deterministic/optim/cases_averted_table_det.txt", sep = ""))

# table_flextable <- flextable(comb_CA) %>%
#   colformat_num(j = c("Number.of.Cases", "Cases.Averted"), big.mark = ",", digits = 0) %>%
#   colformat_num(j = "Percent.Reduction", digits = 2) %>%
#   set_header_labels(
#     `Number of Cases` = "Number of Cases",
#     `Cases Averted` = "Cases Averted",
#     `Percent Reduction` = "% Reduction"
#   ) %>%
#   theme_box() %>%
#   autofit() %>%
#   set_table_properties(layout = "autofit") %>%
#   add_header_row(values = c("", "Cases Averted in Moissala"), colwidths = c(1, 3))
# 
# # Save table as .svg
# save_as_image(table_flextable, 
#               path = paste(plot_loc, "results/cases_averted_table.svg", sep = ""))

################################################################################
### ------------------- CHANGING DEPLOYMENT STRATEGY ----------------------- ###
################################################################################
# When adding additional rounds, assumptions on coverage have to be made
# If at the beginning, repeat coverage of first round twice
# If at the end, repeat coverage of last round twice
start_date <- ymd("2014-01-01")
end_date <- ymd("2022-12-31")
dates_sim <- ymd("2012-01-01")
dates_obs <- ymd("2022-12-31")

years <- c(2014:2022)

################################################################################
### ------------------- AVERAGING CLIMATE PARAMETERS ----------------------- ###
################################################################################
moiss_met_2014_2023 <- readRDS(paste(data_dir, "Data/climate-data/moiss_met_2014_2023_20.rds", sep = ""))
moiss_met_360 <- climate_to_30_day_months(moiss_met_2014_2023, start_year = 2014, end_year = 2022)
moiss_met_360$day <- rep(1:360, 9)
moiss_met_360_1yr <- moiss_met_360 %>% group_by(day) %>% summarize(rollmean = mean(rollmean), temp = mean(temp))
moiss_met_360_1yr$anom <- (moiss_met_360_1yr$rollmean - mean(moiss_met_360_1yr$rollmean)) / (sd(moiss_met_360_1yr$rollmean))
moiss_met_360_1yr$temp <- mean(moiss_met_360_1yr$temp)
moiss_met_360_1yr
moiss_met_360_rep <- rbind(moiss_met_360_1yr, do.call(rbind, replicate(8, moiss_met_360_1yr, simplify = FALSE)))
moiss_met_360_rep$day <- 1:nrow(moiss_met_360_rep)

################################################################################
### ---------------------- AVERAGING COVERAGE DATA ------------------------- ###
################################################################################
smc_cov <- read_excel("C:/Users/putnni/switchdrive/Chad/Data/data-raw/CPS/CPS_coverage_imput_2018.xlsx")
smc_cov_df <- data.frame(date_start = smc_cov$date_start, 
                         cov_total = smc_cov$smc_couv_tot)
smc_cov_df$month <- month(smc_cov_df$date_start)
smc_cov_avg <- smc_cov_df %>% group_by(month) %>% summarize(smc_couv_tot = mean(cov_total))
smc_cov_avg[6,1] <- 6
smc_cov_avg[6,2] <- 0.890
smc_cov_avg <- smc_cov_avg[order(smc_cov_avg$month),]
smc_cov_avg$date_start <- as.Date(c("2014-06-14", "2014-07-14", "2014-08-11", 
                            "2014-09-15", "2014-10-13", "2014-11-13"))

# Repeat the dataframe for each year and adjust the dates
smc_cov_avg_expanded <- bind_rows(
  lapply(years, function(y) {
    smc_cov_avg %>%
      mutate(date_start = as.Date(format(date_start, paste0(y, "-%m-%d"))))
  })
)

# Add the original year 2014 data
smc_cov_avg_final <- bind_rows(smc_cov_avg, smc_cov_avg_expanded)
smc_cov_avg <- smc_cov_avg_final[-1]

inputs_moiss$c_R_D <- moiss_met_360_rep$anom
inputs_moiss$temp <- moiss_met_360_rep$temp

################################################################################
### --------------------- 4 ROUNDS START IN JULY --------------------------- ###
################################################################################
smc_cov_1 <- subset(smc_cov_avg, !(month(smc_cov_avg$date_start) %in% c(6, 11)))
smc_results_1 <- generate_smc_schedule_df(start_date, end_date, years, months_30_days = TRUE, smc_cov = smc_cov_1, const = -0.1806)
smc_1 <- smc_results_1$SMC
cov_1 <- smc_results_1$cov
decay_1 <- smc_results_1$decay
inputs_moiss[c("SMC", "cov_SMC", "decay")] <- c(list(smc_1), list(cov_1), list(decay_1))
sim_1 <- data_sim(model_SMC, inputs_moiss,
                    start_date = start_date, end_date = end_date,
                    month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)

################################################################################
### -------------------- 5 ROUNDS START IN JULY ---------------------------- ###
################################################################################
smc_cov_2 <- subset(smc_cov_avg, !(month(smc_cov_avg$date_start) %in% c(6)))
smc_results_2 <- generate_smc_schedule_df(start_date, end_date, years, months_30_days = TRUE, smc_cov = smc_cov_2, const = -0.1806)
smc_2 <- smc_results_2$SMC
cov_2 <- smc_results_2$cov
decay_2 <- smc_results_2$decay
inputs_moiss[c("SMC", "cov_SMC", "decay")] <- c(list(smc_2), list(cov_2), list(decay_2))
sim_2 <- data_sim(model_SMC, inputs_moiss,
                  start_date = start_date, end_date = end_date,
                  month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)

################################################################################
### -------------------- 5 ROUNDS START IN JUNE ---------------------------- ###
################################################################################
smc_cov_3 <- subset(smc_cov_avg, !(month(smc_cov_avg$date_start) %in% c(11)))
smc_results_3 <- generate_smc_schedule_df(start_date, end_date, years, months_30_days = TRUE, smc_cov = smc_cov_3, const = -0.1806)
smc_3 <- smc_results_3$SMC
cov_3 <- smc_results_3$cov
decay_3 <- smc_results_3$decay
inputs_moiss[c("SMC", "cov_SMC", "decay")] <- c(list(smc_3), list(cov_3), list(decay_3))
sim_3 <- data_sim(model_SMC, inputs_moiss,
                  start_date = start_date, end_date = end_date,
                  month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)

################################################################################
### ---------------------- SUMMARIZING RESULTS ----------------------------- ###
################################################################################
# Plotting trajectories for each simulation
colnames(sim_1) <- c("date_ymd", "month_no", "o5", "u5", "total")
colnames(sim_2) <- c("date_ymd", "month_no", "o5", "u5", "total")
colnames(sim_3) <- c("date_ymd", "month_no", "o5", "u5", "total")

# Add simulation identifier
sim_1$simulation <- "4RoundsJuly"
sim_2$simulation <- "5RoundsJuly"
sim_3$simulation <- "5RoundsJune"

# Combine into one dataframe
combined_sim <- bind_rows(sim_1[-c(1:5),], sim_2[-(1:5),], sim_3[-(1:5),])

# Pivot longer to reshape data
combined_sim_long <- combined_sim %>%
  pivot_longer(cols = c(u5, total), names_to = "type", values_to = "value")

# Create the plot
plot <- ggplot(combined_sim_long, aes(x = date_ymd, y = value, color = type)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ simulation, scales = "free_y", ncol = 1) +
  labs(title = "",
       x = "Date",
       y = "Monthly Cases",
       color = "Type of Incidence") +
  theme_bw() + theme_plots

ggsave(paste(plot_loc, "results/deterministic/optim/incidence_simulations_plot_det.pdf", sep = ""), plot, width = 10, height = 6)

# Calculating total cases for each simulation
calc_total_cases <- function(sim) {
  #months_2019 <- ((2019 - 2014) * 12) : ((2020 - 2014) * 12)
  total <- sum(sim$total)
  u5 <- sum(sim$u5)
  return(list(total = total, u5 = u5))
}

results <- tibble(
  Simulation = c("4RoundsJuly", "5RoundsJuly", "5RoundsJune"),
  Total_Cases = c(calc_total_cases(sim_1)$total, calc_total_cases(sim_2)$total, calc_total_cases(sim_3)$total),
  U5_Cases = c(calc_total_cases(sim_1)$u5, calc_total_cases(sim_2)$u5, calc_total_cases(sim_3)$u5)
)

# Convert results to long format for easier plotting
results_long <- results %>%
  pivot_longer(cols = c(Total_Cases, U5_Cases), names_to = "Metric", values_to = "Cases")

# Create the plot
plot <- ggplot(results_long, aes(x = Simulation, y = Cases, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Simulation",
       y = "Monthly Cases",
       fill = "Metric") +
  theme_minimal() + theme_plots

ggsave(paste(plot_loc, "results/deterministic/simulations_summary_plot_det.pdf", sep = ""), plot, width = 8, height = 6)

results_tex <- xtable(results)
print(results_tex, file=paste(plot_loc, "results/deterministic/optim/cases_optim_table_det.txt", sep = ""))

# # Making summary table
# table_flextable <- flextable(results) %>%
#   colformat_num(j = c("Total_Cases", "U5_Cases"), big.mark = ",", digits = 0) %>%
#   set_header_labels(
#     Simulation = "Simulation",
#     Total_Cases = "Total Cases",
#     U5_Cases = "Under-5 Cases"
#   ) %>%
#   theme_box() %>%
#   autofit() %>%
#   set_table_properties(layout = "autofit") %>%
#   add_header_row(values = c("", "Summary of Cases by Simulation"), colwidths = c(1, 2))
# 
# # Save table as .svg
# save_as_image(table_flextable, 
#               path = paste(plot_loc, "results/SMC_simulations_summary_table.svg", sep = ""))
