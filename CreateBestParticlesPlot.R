################################################################################
### ----------- Simulating multiple particles using best params ------------ ###
################################################################################
dir_clim <-  "C:/Users/putnni/switchdrive/Chad/Data/climate-data/"
data_dir <- "C:/Users/putnni/switchdrive/Chad/"
code_dir <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/"
plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"
source(paste(dir, "aux_functions.R", sep = ""))
path_full_model_SMC <- paste(code_dir, "models/full_model_SMC_stochastic.R", sep = "")
model_SMC_stoch <- odin.dust::odin_dust(path_full_model_SMC)

dates <- c("2006-01-01", "2022-12-31")
start_date <- ymd("2006-01-01")
end_date <- ymd("2022-12-31")
years <- c(2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 
           2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)

################################################################################
### ------------------------- SMC Inputs ----------------------------------- ###
################################################################################
load(paste(data_dir, "Data/model-inputs/moiss_inputs.Rdata", sep = ""))
inputs_moiss <- read_excel(paste(data_dir, 
                                 "estimated-and-fixed-values/moiss_fixed.xlsx", 
                                 sep = ""))
inputs_moiss <- setNames(as.list(as.numeric(inputs_moiss$Value)), inputs_moiss$Name)
inputs_moiss[c("decay", "SMC", "cov_SMC", "c_R_D", "temp")] <- 
  c(list(decay_moiss), list(SMC_moiss), list(cov_moiss), list(anom_360_moiss), 
    list(temp_360_moiss))
estimated_values <- readRDS("C:/Users/putnni/switchdrive/Chad/estimated-and-fixed-values/moiss_estimated_stoch.rds")
inputs_moiss[rownames(estimated_values)] <- estimated_values$best_params # fixing values to those from estimated for inference procedure

################################################################################
### ------------------------- Running Model -------------------------------- ###
################################################################################
plot_particles <- function(inputs_moiss, title){
  months <- seq(ymd("2014-01-01"), ymd("2022-12-31"), by = "month")
  obs_start <- 6120 - 9 * 360
  mod <- model_SMC_stoch$new(pars = inputs_moiss, 
                             time = obs_start, 
                             n_particles = 10)
  
  n_steps <- length(anom_360_moiss)
  results <- mod$simulate((obs_start + 1):n_steps) # Simulate for 1 year
  ind_month <- seq(30, n_steps - obs_start, by = 30)
  ind_total <- mod$info()$index$month_inc_total
  ind_C <- mod$info()$index$month_inc_C
  ind_A <- mod$info()$index$month_inc_A
  
  month_inc <- results[c(ind_C, ind_A, ind_total),,][,,ind_month]
  
  # Plot the monthly incidence
  matplot(months, t(month_inc[3,,]), type = "l",
          xlab = "", ylab = "Monthly Cases", main = title,
          lty = 1, ylim = range(month_inc))
}
# Combine the plots using par(mfrow)
pdf(paste(plot_loc, "best_particle_plots.pdf", sep = "")) # Save to PDF
par(mfrow = c(2, 2)) # 2x2 layout

# Generate and plot for different values of N
inputs_moiss$N <- 800000.00
plot_particles(inputs_moiss, "N = 800,000")

inputs_moiss$N <- 100000
plot_particles(inputs_moiss, "N = 100,000")

inputs_moiss$N <- 10000
plot_particles(inputs_moiss, "N = 10,000")

inputs_moiss$N <- 1000
plot_particles(inputs_moiss, "N = 1,000")

dev.off() # Close the PDF device
