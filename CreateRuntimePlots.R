plot_loc <- "C:/Users/putnni/Documents/MastersThesisTPH/MastersThesisGit/plots/"
# Create the data frame
data <- data.frame(
  n_particles = c(1, 5, 10, 50),
  Time = c(37.85855, 69.65125, 110.06985, 433.83124)
)

theme_plots <- theme(plot.title = element_text(hjust=0.5, size=20),
                     axis.text.x = element_text(size = 14),
                     axis.title.x = element_text(size = 16),
                     axis.title.y = element_text(size = 16),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size = 14)) 


# Create the line graph
p <- ggplot(data, aes(x = n_particles, y = Time)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  labs(title = "Execution Time vs Number of Particles",
       x = "Number of Particles",
       y = "Time (secs)") +
  theme_minimal() + theme_plots

# Save the plot
ggsave(paste(plot_loc, "execution_time_plot.png", sep = ""), plot = p, width = 8, height = 6)