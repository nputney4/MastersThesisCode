library(GGally)
library(ggplot2)
library(dplyr)
library(flextable)
library(xtable)
library(coda)
library(scales)

analyze_results <- function(results, model_SMC, stochastic, corr = TRUE, ppc = TRUE) {
  data_dir <- "C:/Users/putnni/switchdrive/Chad/"
  code_dir <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/"
  plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"
  
  params_to_estimate_moiss <- c(a_R = "a_R", b_R = "b_R", s = "s", size = "size",
                                qR = "qR", z = "z", phi = "phi", eff_SMC = "eff_SMC")
  
  if (stochastic) {
    result_path <- "results/stochastic/"
    file_suffix <- "stoch"
  } else {
    result_path <- "results/deterministic/"
    file_suffix <- "det"
  }
  
  # Posterior plots
  source(paste(code_dir, "misc/posterior_plots_helper.R", sep = ""))
  dim_plot <- c(3, 3)
  p1 <- post_plot(results, params_to_estimate_moiss, dim_plot = dim_plot, show_true = FALSE)
  ggsave(p1, filename = paste0(plot_loc, result_path, "post_distributions_", file_suffix, ".pdf"), device = "pdf",
         height = 4, width = 6, units = "in")
  
  # Trace plots
  pdf(paste0(plot_loc, result_path, "trace_plots_", file_suffix, ".pdf"), width = 6, height = 4)
  # Adjust plot layout parameters
  par(mfrow = c(3, 3),   # Example: 2 rows and 2 columns layout (adjust as needed)
      mar = c(3, 3, 3, 3), # Margin size: bottom, left, top, right (adjust as needed)
      oma = c(0, 0, 0, 0)) # Outer margins (space around the whole plot area)
  
  # Generate plots from the mcmc.list object
  chains <- plot_chains(results[[1]])
  
  # Plot each chain
  plot(chains, density = FALSE, auto.layout = TRUE)
  # Close the PDF device
  dev.off()
  
  # Gelman Rubin diagnostic
  extracted_chains <- lapply(chains, function(chain) {
    as.mcmc(chain[, 1, drop = FALSE])
  })
  extracted_chains <- mcmc.list(extracted_chains)
  gelman_diag <- gelman.diag(chains)
  gel_latex <- xtable(gelman_diag[[1]])
  
  if(stochastic){caption <- "Stochastic model"}else{caption <- "Deterministic model"}
  print(gel_latex, file = paste0(plot_loc, result_path, "gelman_diag_", file_suffix, ".txt"))
  
  
  pdf(paste0(plot_loc, result_path, "gelman_plot_", file_suffix, ".pdf"), width = 6, height = 4)
  par(cex.main = 2,  # Size of the main title
      cex.lab = 1.6,   # Size of the axis labels
      cex.axis = 1.6)  # Size of the axis tick labels
  gelman.plot(extracted_chains)
  dev.off()
  
  # Effective sample size
  eff_SS <- data.frame(effectiveSize(results$coda_pars[, -c(1, 2, 3)]))
  colnames(eff_SS) <- "Effective SS"
  eff_SS_latex <- xtable(t(eff_SS))
  
  print(eff_SS_latex, file = paste0(plot_loc, result_path, "eff_SS_", file_suffix, ".txt"))
  
  # Autocorrelation plots
  pdf(paste0(plot_loc, result_path, "acf_plot_", file_suffix, ".pdf"), width = 6, height = 4)
  par(mfrow = c(3, 3),   # Example: 2 rows and 2 columns layout (adjust as needed)
      mar = c(3, 3, 3, 3), # Margin size: bottom, left, top, right (adjust as needed)
      oma = c(0, 0, 0, 0)) # Outer margins (space around the whole plot area)
  autocorr.plot(results$coda_pars[,-c(1,2,3)])
  dev.off()
  # Correlation plots
  p2 <- plot_corr(results)
  ggsave(p2, filename = paste0(plot_loc, result_path, "corr_plot_", file_suffix, ".pdf"), device = "pdf",
         height = 4.4, width = 6, units = "in")
  
  # Posterior predictive check
  dates_sim <- c("2006-01-01", "2022-12-31")
  dates_obs <- c("2014-01-01", "2022-12-31")
  sim_df <- create_sim_df(results, 300, dates_sim, dates_obs, model_SMC)
  clim_df <- create_clim_df(moiss_met_2014_2023)
  p3 <- post_pred_plot(results, sim_df, dates_obs, dates_obs, show_clim = TRUE, clim_df = clim_df,
                       title = "", theme = theme_plots, ages = c("o5", "u5", "total"))
  ggsave(p3, filename = paste0(plot_loc, result_path, "post_predictive_check_", file_suffix, ".pdf"), device = "pdf",
         height = 4, width = 8.3, units = "in")
}



analyze_results_synthetic <- function(results, model_SMC, file_suffix, params_to_estimate, true_values) {
  data_dir <- "C:/Users/putnni/switchdrive/Chad/"
  code_dir <- "C:/Users/putnni/Documents/MastersThesisTPH/ProjectGithub/moissala-smc/"
  plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"
  result_path <- "synthetic/test_plots/"
  # Posterior plots
  source(paste(code_dir, "misc/posterior_plots_helper.R", sep = ""))
  dim_plot <- c(3, 3)
  p1 <- post_plot(results, params_to_estimate, dim_plot = dim_plot, show_true = TRUE, true_value = true_values)
  ggsave(p1, filename = paste0(plot_loc, result_path, "post_distributions_", file_suffix, ".pdf"), device = "pdf",
         height = 4, width = 6, units = "in")
  
  # Trace plots
  pdf(paste0(plot_loc, result_path, "trace_plots_", file_suffix, ".pdf"), width = 6, height = 4)
  # Adjust plot layout parameters
  par(mfrow = c(3, 3),   # Example: 2 rows and 2 columns layout (adjust as needed)
      mar = c(3, 3, 3, 3), # Margin size: bottom, left, top, right (adjust as needed)
      oma = c(0, 0, 0, 0)) # Outer margins (space around the whole plot area)
  
  # Generate plots from the mcmc.list object
  chains <- plot_chains(results[[1]])
  
  # Plot each chain
  plot(chains, density = FALSE, auto.layout = TRUE)
  # Close the PDF device
  dev.off()
  
  # Gelman Rubin diagnostic
  extracted_chains <- lapply(chains, function(chain) {
    as.mcmc(chain[, 1, drop = FALSE])
  })
  extracted_chains <- mcmc.list(extracted_chains)

  # Attempt to calculate the Gelman-Rubin diagnostic
  gelman_diag <- tryCatch({
    gelman.diag(chains)
  }, error = function(e) {
    # Handle the error (print a message or take any other action)
    message("Warning in gelman.diag: ", e$message)
    NULL  # Return NULL or any other default value
  })
  
  # You can then check if gelman_diag is NULL before using it
  if (!is.null(gelman_diag)) {
    gel_latex <- xtable(gelman_diag[[1]])
    print(gel_latex, file = paste0(plot_loc, result_path, "gelman_diag_", file_suffix, ".txt"))
  }
  
  #gelman_diag <- gelman.diag(chains)
  #gel_latex <- xtable(gelman_diag[[1]])

  #print(gel_latex, file = paste0(plot_loc, result_path, "gelman_diag_", file_suffix, ".txt"))
  
  
  pdf(paste0(plot_loc, result_path, "gelman_plot_", file_suffix, ".pdf"), width = 6, height = 4)
  par(cex.main = 2,  # Size of the main title
      cex.lab = 1.6,   # Size of the axis labels
      cex.axis = 1.6)  # Size of the axis tick labels
  gelman.plot(extracted_chains)
  dev.off()
  
  # Effective sample size
  eff_SS <- data.frame(effectiveSize(results$coda_pars[, -c(1, 2, 3)]))
  colnames(eff_SS) <- "Effective SS"
  eff_SS_latex <- xtable(t(eff_SS))
  
  print(eff_SS_latex, file = paste0(plot_loc, result_path, "eff_SS_", file_suffix, ".txt"))
  
  # Autocorrelation plots
  pdf(paste0(plot_loc, result_path, "acf_plot_", file_suffix, ".pdf"), width = 6, height = 4)
  par(mfrow = c(3, 3),   # Example: 2 rows and 2 columns layout (adjust as needed)
      mar = c(3, 3, 3, 3), # Margin size: bottom, left, top, right (adjust as needed)
      oma = c(0, 0, 0, 0)) # Outer margins (space around the whole plot area)
  autocorr.plot(results$coda_pars[,-c(1,2,3)])
  dev.off()
  
  # Correlation plots
  # p2 <- plot_corr(results)
  # ggsave(p2, filename = paste0(plot_loc, result_path, "corr_plot_", file_suffix, ".pdf"), device = "pdf",
  #        height = 4.4, width = 6, units = "in")
  
}


# # Set up PDF output with width for two plots side by side and a compact height
# pdf(paste(dir, "trace_plots.pdf", sep = ""), width = 4, height = 3.5)  # Width adjusted for 2 columns, height compact
# 
# par(mfrow = c(3, 3),   # Example: 2 rows and 2 columns layout (adjust as needed)
#     mar = c(2.5, 2.5, 2.5, 2.5), # Margin size: bottom, left, top, right (adjust as needed)
#     oma = c(0, 0, 0, 0)) # Outer margins (space around the whole plot area)
# 
# plot(chains_list, 
#      density = FALSE,              # Disable density plots
#      cex.main = 0.5,                # Main title font size
#      cex.lab = 0.7,                 # Axis labels font size
#      cex.axis = 1.4)                # Axis values font size
# 
# # Close PDF device
# dev.off()
# 
# 
# dir <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/synthetic/test_plots/"
# pdf(paste(dir, "trace_plots.pdf", sep = ""), width = 10, height = 10)  # Adjust dimensions as needed
# plot(chains_list, density = FALSE)
# dev.off()