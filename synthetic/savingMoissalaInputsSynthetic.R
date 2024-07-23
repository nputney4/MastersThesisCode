################################################################################
### ----------- LOADING MODEL AND AUXILIARY FUNCTIONS ---------------------- ###
################################################################################
dir <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/"
dir_data <- "C:/Users/putnni/switchdrive/Chad/Data/cases-data/"
dir_clim <-  "C:/Users/putnni/switchdrive/Chad/Data/climate-data/"
data_dir <- "C:/Users/putnni/switchdrive/Chad/"
source(paste(dir, "aux_functions.R", sep = ""))

################################################################################
### -------------- GENERATING SYNTHETIC SMC SCHEDULES ---------------------- ###
################################################################################
dates <- c("2006-01-01", "2022-12-31")
start_date <- ymd("2006-01-01")
end_date <- ymd("2022-12-31")
smc_cov <- read_excel("C:/Users/putnni/switchdrive/Chad/Data/data-raw/CPS/CPS_coverage_imput_2018.xlsx")
years <- c(2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)

################################################################################
### ------------------ OUTLINING TWO SMC SCHEDULES ------------------------- ###
################################################################################
months_active_alt <- matrix(data = 0, nrow = length(years), ncol = 12) # SMC every year
months_active_alt[, c(7, 8, 9, 10)] <- 1
months_active_real <- months_active_alt # SMC every year but 2019
months_active_real[which(years == 2019),] <- 0  

################################################################################
### ---------------------- CREATING SCHEDULES ------------------------------ ###
################################################################################
SMC_schedule_real <- generate_smc_schedule_df_synthetic(start_date, end_date, years, months_active = months_active_real, months_30_days = TRUE, coverage = 0.90)
SMC_schedule_alt <- generate_smc_schedule_df_synthetic(start_date, end_date, years, months_active = months_active_alt, months_30_days = TRUE, coverage = 0.90)

################################################################################
### ------------------- READING IN CLIMATE INPUTS -------------------------- ###
################################################################################
moiss_met_2006_2023 <- readRDS(file = paste(dir_clim, "moiss_met_2006_2023_20.rds", sep = ""))
moiss_met_2014_2023 <- readRDS(file = paste(dir_clim, "moiss_met_2014_2023_20.rds", sep = ""))
moiss_met_360 <- climate_to_30_day_months(moiss_met_2014_2023, start_year = 2014, end_year = 2022)

# 365 day years
temp_moiss <- moiss_met_2006_2023$temp
anom_moiss <- moiss_met_2006_2023$anom

# 360 day years
temp_360_moiss <- moiss_met_360
anom_360_moiss <- moiss_met_360$anom

# repeating mean temperature across every day
mean_temp_moiss <- mean(temp_moiss)
temp_360_moiss <- rep(mean_temp_moiss, nrow(temp_360_moiss))

################################################################################
### ------------------- FIXED AND ESTIMATED VALUES ------------------------- ###
################################################################################
inputs_moiss_real <- read_excel(paste(data_dir, 
                                 "estimated-and-fixed-values/moiss_fixed.xlsx", 
                                 sep = ""))
inputs_moiss_real <- setNames(as.list(as.numeric(inputs_moiss_real$Value)), inputs_moiss_real$Name)
estimated_values <- readRDS("C:/Users/putnni/switchdrive/Chad/estimated-and-fixed-values/moiss_estimated.rds")
inputs_moiss_real[rownames(estimated_values)] <- estimated_values$best_params # fixing values to those from estimated for inference procedure

inputs_moiss_alt <- inputs_moiss_real

inputs_moiss_real[c("decay", "SMC", "cov_SMC", "c_R_D", "temp")] <- 
  c(list(SMC_schedule_real$decay), list(SMC_schedule_real$SMC), list(SMC_schedule_real$cov), list(anom_360_moiss), 
    list(temp_360_moiss))

inputs_moiss_alt[c("decay", "SMC", "cov_SMC", "c_R_D", "temp")] <- 
  c(list(SMC_schedule_alt$decay), list(SMC_schedule_alt$SMC), list(SMC_schedule_alt$cov), list(anom_360_moiss), 
    list(temp_360_moiss))

################################################################################
### ---------------------- SAVING MOISSALA INPUTS -------------------------- ###
################################################################################
data_dir <- "C:/Users/putnni/switchdrive/Chad/"
save(inputs_moiss_real, inputs_moiss_alt, moiss_met_2014_2023, 
     anom_360_moiss, temp_360_moiss, file = "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/synthetic/inputs_synthetic")
