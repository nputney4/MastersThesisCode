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
# When adding additional rounds, assumptions on coverage have to be made
# If at the beginning, repeat coverage of first round twice
# If at the end, repeat coverage of last round twice

start_date <- ymd("2014-01-01")
end_date <- ymd("2022-12-31")
years <- c(2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)

################################################################################
### --------------------- TRUE DEPLOYMENT STRATEGY ------------------------- ###
################################################################################


################################################################################
### ---------  -- IF ALL CYCLES WERE 5 ROUNDS EXTRA AT START --------------- ###
################################################################################

################################################################################
### -------------- IF ALL CYCLES WERE 5 ROUNDS EXTRA AT END ---------------- ###
################################################################################

################################################################################
### --------------------- IF ALL CYCLES WERE 4 ROUNDS ---------------------- ###
################################################################################



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