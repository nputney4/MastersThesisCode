# Load necessary libraries
library(ggplot2)
library(gridExtra)

theme_pop <- theme(plot.title = element_text(hjust=0.5, size=16),
                   axis.text.x = element_text(size = 10),
                   axis.title.y = element_text(size = 16),
                   legend.title = element_text(size=16),
                   legend.text = element_text(size = 14))

# Function to plot beta distribution
plot_beta <- function(alpha, beta, title) {
  x <- seq(0, 1, length.out = 1000)
  y <- dbeta(x, alpha, beta)
  data <- data.frame(x = x, y = y)
  ggplot(data, aes(x, y)) +
    geom_line() +
    ggtitle(title) +
    xlab('') +
    ylab('Density') + theme_bw() + theme_pop
}

# Function to plot normal distribution
plot_normal <- function(mu, sigma, title, xlim_vals) {
  x <- seq(xlim_vals[1], xlim_vals[2], length.out = 1000)
  y <- dnorm(x, mu, sigma)
  data <- data.frame(x = x, y = y)
  ggplot(data, aes(x, y)) +
    geom_line() +
    ggtitle(title) +
    xlab('') +
    ylab('Density') + theme_bw() + theme_pop
}

# Function to plot gamma distribution
plot_gamma <- function(alpha, beta, title) {
  x <- seq(0, 1, length.out = 1000)
  y <- dgamma(x, alpha, scale = 1/beta)
  data <- data.frame(x = x, y = y)
  ggplot(data, aes(x, y)) +
    geom_line() +
    ggtitle(title) +
    xlab('') +
    ylab('Density') + theme_bw() + theme_pop
}

# Function to plot uniform distribution
plot_uniform <- function(min_val, max_val, title) {
  x <- seq(min_val, max_val, length.out = 1000)
  y <- dunif(x, min_val, max_val)
  data <- data.frame(x = x, y = y)
  ggplot(data, aes(x, y)) +
    geom_line() +
    ggtitle(title) +
    xlab('') +
    ylab('Density') + theme_bw() + theme_pop
}

# Generate plots for each prior
p1 <- plot_beta(2, 2, 'SMC Effectiveness')
p2 <- plot_beta(6, 12, 'Mosquito Egg-adult Sensitivity to Moisture')
p3 <- plot_normal(2, 2, 'Rainfall for 50% Mosquito Survival', c(-5, 5))
p4 <- plot_gamma(20, 20, 'Population Scaling Factor')
p5 <- plot_beta(2, 18, 'Relative Infectivity of Sub-Patent Infection')
p6 <- plot_beta(4, 2, 'Relative Proportion of Adults Receiving Treatment')
p7 <- plot_beta(0.25, 2.065, 'Relative Force of Infection of Adults')
p8 <- plot_normal(0.91, 0.01, 'Daily Probability of Mosquito Survival', c(0.88, 0.95))
p9 <- plot_uniform(4, 10, 'Negative Binomial Dispersion Parameter')

# Combine all plots into one
p_all <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 2)

# Save plots
plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"
ggsave(p_all, filename = paste(plot_loc, "results/prior_plots.pdf", sep = ""), device = "pdf",
       heigh = 12, width = 11, units = "in")

