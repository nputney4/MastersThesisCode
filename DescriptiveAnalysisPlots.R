################################################################################
### ------------- DESCRIPTIVE ANALYSIS OF NMCP DATA ------------------------ ###
################################################################################
# Loading packages
pacman::p_load(
  rio,          # File import
  here,         # File locator
  skimr,        # get overview of data
  tidyverse,    # data management + ggplot2 graphics 
  gtsummary,    # summary statistics and tests
  rstatix,      # summary statistics and statistical tests
  janitor,      # adding totals and percents to tables
  scales,       # easily convert proportions to percents  
  flextable,     # converting tables to pretty images
  gridExtra
)

################################################################################
### ----------------------- READING IN DATA  ------------------------------- ###
################################################################################
plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"
chad_dir <- "C:/Users/putnni/switchdrive/Chad/Data/"
# Reading in data
PNLP_chad <- as.data.frame(readRDS(paste(chad_dir, "data_PNLP_added_2018.rds", sep = "")))
PNLP_chad <- PNLP_chad %>% filter(variable != "population") # removing population variable

# climate data
dir_clim <- "C:/Users/putnni/switchdrive/Chad/Data/climate-data/"
# reading in climate data
bedjo_met <-  as.data.frame(readRDS(paste(dir_clim, "bedjo_met_2014_2023.rds", sep = "")))
goundi_met <-  as.data.frame(readRDS(paste(dir_clim, "goundi_met_2014_2023.rds", sep = "")))
koumra_met <-  as.data.frame(readRDS(paste(dir_clim, "koumra_met_2014_2023.rds", sep = "")))
moiss_met <-  as.data.frame(readRDS(paste(dir_clim, "moiss_met_2014_2023.rds", sep = "")))

# demographic data
dir_demo <- "C:/Users/putnni/switchdrive/Chad/Data/demographic-data/"
pop_estimates <- readRDS(paste(dir_demo, "mandoul_estimated_population", sep = ""))

################################################################################
### ------------------- SEPARATE CASES BY AGE GROUP ------------------------ ###
################################################################################
# selecting only districts with complete data (examined elsewhere)
districts_to_keep <- c("BEDJONDO", "GOUNDI", "KOUMRA", "MOISSALA")
PNLP_chad <- PNLP_chad[PNLP_chad$district %in% districts_to_keep,]

# grouping cases for those 5 and under
u5_cases_by_district <- PNLP_chad %>% 
  filter(variable == "consult.malaria.suspected" & (subpop == "1-4 years" | subpop == "0-11 months")) %>%
  group_by(date_ymd, district) %>% summarize(cases_u5 = sum(value)) %>% 
  mutate(year = as.numeric(format(date_ymd, "%Y")))


# grouping cases for those over 5
o5_cases_by_district <- PNLP_chad %>% 
  filter(variable == "consult.malaria.suspected" & (subpop == "5-14 years" | subpop == ">15 years")) %>%
  group_by(date_ymd, district) %>% summarize(cases_o5 = sum(value)) %>% 
  mutate(year = as.numeric(format(date_ymd, "%Y")))

# grouping cases across all age groups
total_cases_by_district <- PNLP_chad %>% 
  filter(variable == "consult.malaria.suspected" & (subpop %in% c("0-11 months","1-4 years", "5-14 years", ">15 years"))) %>%
  group_by(date_ymd, district) %>% summarize(cases_total = sum(value)) %>% 
  mutate(year = as.numeric(format(date_ymd, "%Y")))

################################################################################
### ------------------- ADDING INCIDENCE COLUMN ---------------------------- ###
################################################################################
u5_cases_by_district <- u5_cases_by_district %>% 
  left_join(pop_estimates %>% filter(subpop == "0-4 years") %>% 
              dplyr::select(-date_ymd), by = c("district", "year"))  %>% 
  mutate(inc_u5 = cases_u5 / value) %>% dplyr::select(c(date_ymd, district, year, cases_u5, inc_u5))

o5_cases_by_district <- o5_cases_by_district %>% 
  left_join(pop_estimates %>% filter(subpop == ">4 years") %>% 
              dplyr::select(-date_ymd), by = c("district", "year"))  %>% 
  mutate(inc_o5 = cases_o5 / value) %>% dplyr::select(c(date_ymd, district, year, cases_o5, inc_o5))

total_cases_by_district <- total_cases_by_district %>% 
  left_join(pop_estimates %>% filter(subpop == "all_ages") %>% 
              dplyr::select(-date_ymd), by = c("district", "year"))  %>% 
  mutate(inc_total = cases_total / value) %>% dplyr::select(c(date_ymd, district, year, cases_total, inc_total))

################################################################################
### ------------------- ADDING SMC STATUS TO CASE DATA --------------------- ###
################################################################################
# creating Yes and No of SMC for Moissala based on smc data
SMC_months <- c("2014-07", "2014-08", "2014-09", "2014-10",
                "2015-07", "2015-08", "2015-09", "2015-10",
                "2016-07", "2016-08", "2016-09", "2016-10",
                "2017-07", "2017-08", "2017-09", "2017-10",
                "2018-07", "2018-08", "2018-09", "2018-10",
                "2020-07", "2020-08", "2020-09", "2020-10",
                "2021-07", "2021-08", "2021-09", "2021-10", "2021-11",
                "2022-07", "2022-08", "2022-09", "2022-10", "2022-11")

SMC_months <- as.Date(paste(SMC_months, "01", sep = "-"))

# adding SMC column
add_SMC_column <- function(cases_by_district, SMC_months) {
  cases_by_district$SMC <- rep("No", nrow(cases_by_district))
  for (i in 1:nrow(cases_by_district)) {
    if ((cases_by_district$date_ymd[i] %in% SMC_months) &
        (cases_by_district$district[i] == "MOISSALA")) {
      cases_by_district$SMC[i] <- "Yes"
    }
  }
  return(cases_by_district)
}

u5_cases_by_district <- add_SMC_column(u5_cases_by_district, SMC_months)
o5_cases_by_district$SMC <- "No"
total_cases_by_district <- add_SMC_column(total_cases_by_district, SMC_months)

################################################################################
### ------------------- TOTAL MONTHLY CASES ACROSS DISTRICTS --------------- ###
################################################################################
color <- "#4a695a"

# bar graph of total cases across all districts
month_inc_plot <- total_cases_by_district %>% ggplot(aes(x = as.Date(date_ymd), y = inc_total * 1000, alpha = SMC)) +
  geom_bar(stat = "identity", position = "dodge") + scale_alpha_manual(values = c(0.7, 1)) +
  facet_wrap(~district, ncol = 1) +
  ylab("Monthly Cases per 1000 Population") + xlab("") + theme_bw() + 
  theme(plot.title = element_text(hjust=0.5, size=12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12))

ggsave(month_inc_plot, filename = paste(plot_loc, "month_incidence_by_district.pdf", sep = ""), device = "pdf",
       heigh = 6, width = 14, units = "in")


################################################################################
### ------------------- TOTAL MONTHLY CASES ACROSS DISTRICTS --------------- ###
################################################################################
