################################################################################
### --------- EXTRACTING RELEVANT INFO FROM SMC SYSTEMATIC REVIEW ---------- ###
################################################################################

# paper found at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10793029/

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(gt)
library(webshot2)

# Create the data frame
data <- data.frame(
  Study = c("Ambe et al.", "Dicko et al.", "Konate et al.", "Sesay et al.", "Tagbor et al.", "Ndiaye et al."),
  Country = c("Nigeria", "Mali", "Burkina Faso", "The Gambia", "Ghana", "Senegal"),
  Age_group = c("3 - 59 months", "3 - 59 months", "3 - 59 months", "6 - 59 months", "3 - 59 months", "3 - 119 months"),
  Comparison = c("SP+AQ vs placebo", "SP+AQ vs placebo", "SP+AQ vs placebo", "SP+AQ vs placebo", "SP+AQ vs placebo", "SP+AQ vs nothing"),
  No_of_cycles = c(4, 3, 3, 3, 5, 5),
  No_enrolled = c("204 SMC, 195 control", "1509 SMC, 1509 placebo", "1509 SMC, 1505 placebo", "639 SMC, 638 placebo", "741 SMC, 749 control", "2245 SMC, 2301 control"),
  Coverage = c("73%", "> 95%", "70-80%", "> 95%", "70%", "> 90%"),
  Rainy_Season = c("Jun - Sept.", "July - Oct.", "July - Oct.", "July - Sept.", "May - Nov.", "May - Oct."),
  Months_given = c("Aug., Sept., Oct., Nov.", "Aug., Sept., Oct.", "Aug., Sept., Oct.", "Sept., Oct., Nov.", "July, Aug., Sept., Oct., Nov.", "July, Aug., Sept., Oct., Nov."),
  Risk_Ratio_U5 = c("0.31 [0.16, 0.61]", "0.17 [0.14, 0.20]", "0.29 [0.26, 0.32]", "0.34 [0.06, 1.89]", "0.62 [0.41, 0.93]", "0.18 [0.15, 0.21]"),
  Risk_Ratio_O5 = c(NA, NA, NA, NA, NA, "0.17 [0.15, 0.20]")
)

# Create the table using gt
table_gt <- gt(data) %>%
  tab_header(
    title = ""
  ) %>%
  cols_label(
    Study = "Study",
    Country = "Country",
    Age_group = "Age Group",
    Comparison = "Comparison",
    No_of_cycles = "No. of Cycles",
    No_enrolled = "No. Enrolled",
    Coverage = "Coverage",
    Rainy_Season = "Rainy Season",
    Months_given = "Months Given",
    Risk_Ratio_U5 = "Risk or Rate Ratio (<5 years)",
    Risk_Ratio_O5 = "Risk or Rate Ratio (>5 years)"
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  fmt_missing(
    columns = everything(),
    missing_text = "-"
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "lightgray"),
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "black",
      weight = px(1)
    ),
    locations = cells_body(columns = everything())
  )

# Save the table as HTML
gtsave(table_gt, paste(plot_loc, "smc_review_table.html", sep = ""))
