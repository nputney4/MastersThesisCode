################################################################################
### ----------- EXAMINING LAG BETWEEN RAINFALL AND INCIDENCE --------------- ###
################################################################################
library(dplyr)
library(ggplot2)
library(patchwork)

dir_data <- "C:/Users/putnni/switchdrive/Chad/Data/cases-data/"
dir_clim <-  "C:/Users/putnni/switchdrive/Chad/Data/climate-data/"
plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"

cases_goundi <- readRDS(paste(dir_data, "cases_GOUNDI.rds", sep = ""))
colnames(cases_goundi) <- c("date", colnames(cases_goundi[2:5]))
goundi_met_30 <- readRDS(file = paste(dir_clim, "goundi_met_2014_2023.rds", sep = ""))

# Min-max scale function
min_max_scale <- function(x, min_val, max_val, plus = 0) {
  max_val = max_val + plus
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max_val - min_val) + min_val
}

# Create a sequence of dates from the minimum to the maximum date in your data
all_dates <- seq.Date(from = min(goundi_met_30$date),
                      to = max(goundi_met_30$date),
                      by = "day")

# Expand cases_goundi to include all dates (with NA for missing cases)
exp_cases <- expand.grid(date = all_dates) %>%
  left_join(cases_goundi, by = c("date"))

# Calculate min and max for cases_total and incidence
min_cases <- min(cases_goundi$inc, na.rm = TRUE)
max_cases <- max(cases_goundi$inc, na.rm = TRUE)

# Scaling rainfall data to match cases
exp_cases$rollmean <- min_max_scale(goundi_met_30$rollmean, min_cases, max_cases, plus = 3500)
exp_cases$rainfall <- min_max_scale(goundi_met_30$rainfall, min_cases, max_cases, plus = 3500)

# Plot with daily rainfall
p1 <- exp_cases %>% ggplot(aes(x = date, y = inc)) +
  geom_bar(stat = "identity", position = "dodge", width = 24) +
  geom_line(aes(group = 1, y = rainfall), color = "blue", linewidth = 0.8) + 
  theme_bw() + xlab("") + ylab("monthly cases") + ggtitle("Daily Rainfall vs 30 Day Moving Average") +
  theme(plot.title = element_text(hjust=0.5, size = 16),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12))


# Plot with 30 day rolling average
p2 <- exp_cases %>% ggplot(aes(x = date, y = inc)) +
  geom_bar(stat = "identity", position = "dodge", width = 24) +
  geom_line(aes(group = 1, y = rollmean), color = "blue", linewidth = 0.8) + 
  theme_bw() + xlab("") + ylab("monthly cases") +
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12))

clim_plot <- p1 / p2 + plot_layout(axis_titles = "collect")

ggsave(clim_plot, filename = paste(plot_loc, "rainfall_lag.pdf", sep = ""), device = "pdf",
       heigh = 4, width = 7, units = "in")
