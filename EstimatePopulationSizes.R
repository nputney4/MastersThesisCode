################################################################################
### ------------- ESTIMATING DISTRICT LEVEL POPULATION SIZES --------------- ###
################################################################################

library(readxl)
library(reshape2)
# assuming uniform population growth across age groups across Chad between 2014 and 2022
# assuming growth between 2014 and 2022 is approximately linear
# estimate previous year total population and 5 and under population in districts of interest
dir_demo <- "C:/Users/putnni/switchdrive/Chad/Data/demographic-data/"
chad_demo <- readxl::read_excel(paste(dir_demo, "chad_demographics.xlsx", sep = ""))
mandoul_demo <- readxl::read_excel(paste(dir_demo, "mandoul_demographics.xlsx", sep = ""))
chad_dir <- "C:/Users/putnni/switchdrive/Chad/Data/"
PNLP_chad <- readRDS(paste(chad_dir, "data_PNLP_2014_2022.rds", sep = ""))
plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"

# function to calculate it for other districts
pred_district_pop <- function(growth_rates, initial_pop, n_years){
  pred_pop <- c(initial_pop, rep(NA, n_years))
  for(i in 1:n_years){
    pred_pop[i+1] <- (1 + growth_rates[i]) * initial_pop
    initial_pop <- pred_pop[i]
  }
  return(pred_pop)
}

# Formatting data in long format to match with larger PNLP dataset
format_pop_data <- function(pop_growth, subpop){
  pop_growth <- setNames(data.frame(pop_growth), c("MOISSALA", "GOUNDI", "KOUMRA", "BEDJONDO"))
  pop_growth <- melt(pop_growth) %>% setNames(c("district", "value"))
  
  pop_growth$year <- rep(2023:2014, length(unique(pop_growth$district)))
  pop_growth$variable <- "population"
  pop_growth$subpop <- subpop
  pop_growth$source <- NA
  
  pop_growth$date_ym <- paste(pop_growth$year, "01", sep = "-") 
  pop_growth$date_ymd <- as.Date(paste(pop_growth$year, "01", "01", sep = "-"))
  pop_growth$month <- "JANVIER"
  
  pop_growth <- pop_growth[colnames(PNLP_chad)] # reordering columns
  return(pop_growth)
}


# fit linear model and assess fit
# y1 <- chad_demo$population_total
# y2 <- chad_demo$population_5_and_under
# year <- chad_demo$Year
# year <- 1:11
# 
# lm1 <- lm(y1 ~ year)  
# 
# plot(year, y1, cex = 2, col = "blue", xlab = "Year", ylab = "Population", main = "Linear Fit to Chad Population")
# lines(year, predict(lm1, as.data.frame(x))) # seems to be very linear in this window

# calculate growth rate for Chad - years going backwords (so negative growth rate)
growth_rate_CHAD_total <- diff(chad_demo$population_total) / chad_demo$population_total[1:10]
growth_rate_CHAD_total <- c(growth_rate_CHAD_total, growth_rate_CHAD_total[length(growth_rate_CHAD_total)]) # repeating 2022 value for 2023 because no data available

growth_rate_CHAD_u5 <- diff(chad_demo$population_5_and_under) / chad_demo$population_5_and_under[1:10]
growth_rate_CHAD_u5 <- c(growth_rate_CHAD_u5, growth_rate_CHAD_u5[length(growth_rate_CHAD_u5)]) # repeating 2022 value for 2023 because no data available

growth_rate_CHAD_o5 <- diff(chad_demo$population_over_5) / chad_demo$population_over_5[1:10]
growth_rate_CHAD_o5 <- c(growth_rate_CHAD_o5, growth_rate_CHAD_o5[length(growth_rate_CHAD_o5)]) # repeating 2022 value for 2023 because no data available

# estimating population for each of four districts of interest in Mandoul region
init_MOISSALA_total <- subset(mandoul_demo, (Year == 2023 & District == "Moissala"))$population_total
init_MOISSALA_u5 <- subset(mandoul_demo, (Year == 2023 & District == "Moissala"))$population_5_and_under
init_MOISSALA_o5 <- subset(mandoul_demo, (Year == 2023 & District == "Moissala"))$population_over_5

init_GOUNDI_total <- subset(mandoul_demo, (Year == 2023 & District == "Goundi"))$population_total
init_GOUNDI_u5 <- subset(mandoul_demo, (Year == 2023 & District == "Goundi"))$population_5_and_under
init_GOUNDI_o5 <- subset(mandoul_demo, (Year == 2023 & District == "Goundi"))$population_over_5

init_KOUMRA_total <- subset(mandoul_demo, (Year == 2023 & District == "Koumra"))$population_total
init_KOUMRA_u5 <- subset(mandoul_demo, (Year == 2023 & District == "Koumra"))$population_5_and_under
init_KOUMRA_o5 <- subset(mandoul_demo, (Year == 2023 & District == "Koumra"))$population_over_5

init_BEDJONDO_total <- subset(mandoul_demo, (Year == 2023 & District == "Bedjondo"))$population_total
init_BEDJONDO_u5 <- subset(mandoul_demo, (Year == 2023 & District == "Bedjondo"))$population_5_and_under
init_BEDJONDO_o5 <- subset(mandoul_demo, (Year == 2023 & District == "Bedjondo"))$population_over_5

# Applying function to all age categories

pop_growth_total <- sapply(c(init_MOISSALA_total, init_GOUNDI_total, init_KOUMRA_total, init_BEDJONDO_total),
                           pred_district_pop, growth_rates = growth_rate_CHAD_total, n_years = 9)

pop_growth_u5 <- sapply(c(init_MOISSALA_u5, init_GOUNDI_u5, init_KOUMRA_u5, init_BEDJONDO_u5),
                        pred_district_pop, growth_rates = growth_rate_CHAD_u5, n_years = 9)

pop_growth_o5 <- sapply(c(init_MOISSALA_o5, init_GOUNDI_o5, init_KOUMRA_o5, init_BEDJONDO_o5),
                        pred_district_pop, growth_rates = growth_rate_CHAD_o5, n_years = 9)

# Formatting data in long format to match with larger PNLP dataset
pop_growth_total_2 <- format_pop_data(pop_growth_total, subpop = "all_ages")
pop_growth_u5_2 <- format_pop_data(pop_growth_u5, subpop = "0-4 years")
pop_growth_o5_2 <- format_pop_data(pop_growth_o5, subpop = ">4 years")

# combining age groups and saving data
pop_pred <- rbind(pop_growth_total_2, pop_growth_u5_2, pop_growth_o5_2)
saveRDS(pop_pred, paste(dir_demo, "mandoul_estimated_population", sep = ""))

# plotting population growth over time
pop_plot <- pop_pred %>% ggplot(aes(x = date_ymd, y = value, color = subpop)) +
  geom_point(size = 3) + geom_line() +
  facet_wrap(~district, ncol = 2, scales = "free") + xlab("") + ylab("population") +
  theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12))

ggsave(pop_plot, filename = paste(plot_loc, "pop_estimates.pdf", sep = ""), device = "pdf",
       heigh = 4, width = 7, units = "in")
