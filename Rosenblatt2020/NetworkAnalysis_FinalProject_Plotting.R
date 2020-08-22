#### Presettings ####

# Load packages

library(ggplot2)
library(tidyverse)
library(dplyr)


# Import dataframes

data_between <- read.csv("data_between.csv", sep = ";")
data_degree <- read.csv("data_degree.csv", sep = ";")
data_eigen <- read.csv("data_eigen.csv", sep = ";")
data_page <- read.csv("data_page.csv", sep = ";")


#### Transformation of data #### 

# Define function to transform data

transform_data <- function(data) {
  data_trans <- data %>%
    group_by(p) %>% 
    summarize(mean_Nie = mean(robustness_Nie), 
              mean_Ros = mean(robustness_Ros),
              std_Nie = sd(robustness_Nie),
              std_Ros = sd(robustness_Ros))
  return(data_trans)
}

# Generate transformed dataframes for each centrality measure using function

data_between_trans <- transform_data(data_between)
data_degree_trans <- transform_data(data_degree)
data_eigen_trans <- transform_data(data_eigen)
data_page_trans <- transform_data(data_page)


#### Robustness_Nie + Robustness_Ros (Same scale) ####

# Define function to plot both robustness measures on the same scale

plot_same_scale <- function(data_trans, title) {
  # Convert dataframe in a format that is suitable for plotting with ggplot
  data_trans <- data_trans %>%
    mutate(Type = "NULL")
  data_plot <- rbind(
  data_trans %>%
    select(p, mean = mean_Nie, std = std_Nie, Type) %>%
    mutate(Type = "Niemeyer"), 
  data_trans %>%
    select(p, mean= mean_Ros, std = std_Ros, Type) %>%
    mutate(Type = "Rosenblatt"))
  # Generate plot 
  ggplot(data_plot, aes(x = p, y = mean, color = Type)) + 
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin = mean-std, ymax = mean+std), 
                  width = 0.02) + 
    scale_y_continuous("Mean robustness") +
    scale_x_continuous("Error level p") +
    scale_color_manual(values = c("#C3A951", "#5683B0")) +
    labs(title = title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}


plot_same_scale <- function(data_trans, title) {
  # Convert dataframe in a format that is suitable for plotting with ggplot
  data_trans <- data_trans %>%
    mutate(Type = "NULL")
  data_plot <- rbind(
    data_trans %>%
      select(p, mean = mean_Nie, std = std_Nie, Type) %>%
      mutate(Type = "Niemeyer"), 
    data_trans %>%
      select(p, mean= mean_Ros, std = std_Ros, Type) %>%
      mutate(Type = "Rosenblatt"))
  # Generate plot 
  ggplot(data_plot, aes(x = p, y = mean, color = Type)) + 
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin = mean-std, ymax = mean+std), 
                  width = 0.02) + 
    scale_y_continuous("Robustness of Centrality Measure", 
                       sec.axis = sec_axis(~.*1, name = "Difference in Outbreak Size vs. Random Immunization")) +
    scale_x_continuous("Proportion of missing nodes p") +
    scale_color_manual(values = c("#C3A951", "#5683B0")) +
    labs(title = title) +
    theme_bw() +
    theme(axis.title.y.left = element_text(colour = "#C3A951"),
          axis.title.y.right = element_text(colour = "#5683B0"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}


# Generate plots for each centrality measure using function


plot_same_scale(data_between_trans, "Betweenness Centrality")
plot_same_scale(data_degree_trans, "Degree Centrality")
plot_same_scale(data_eigen_trans, "Eigenvalue Centrality")
plot_same_scale(data_page_trans, "Pagerank")


#### Robustness_Nie + Robustness_Ros (Different scale) ####

# Define function to plot both robustness measures on different scales

plot_different_scale <- function(data_trans, title) {
  # Define normalizing constant a for robustness measure by Rosenblatt
  a <- max(data_trans$mean_Ros)
  # Convert dataframe in a format that is suitable for plotting with ggplot
  data_trans <- data_trans %>%
    mutate(Type = "NULL")
  data_plot <- rbind(
    data_trans %>%
      select(p, mean_norm = mean_Nie, std_norm = std_Nie, Type) %>%
      mutate(Type = "Niemeyer"), 
    data_trans %>%
      mutate(mean_Ros = mean_Ros/a, std_Ros = std_Ros/a) %>%
      select(p, mean_norm = mean_Ros, std_norm = std_Ros, Type) %>%
      mutate(Type = "Rosenblatt")) 
  
  # Generate plot 
  ggplot(data_plot, aes(x = p, y = mean_norm, color = Type)) + 
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin = mean_norm-std_norm, ymax = mean_norm+std_norm), 
                  width = 0.02) + 
    scale_y_continuous("Robustness of Centrality Measure", 
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       sec.axis = sec_axis(~.*a, name = "Difference in Outbreak Size vs. Random Immunization",
                                           breaks = c(0.05, 0.10, 0.15, 0.20))) +
    scale_x_continuous("Proportion of missing nodes p") +
    scale_color_manual(values = c("#C3A951", "#5683B0")) +
    labs(title = title) +
    theme_bw() +
    theme(axis.title.y.left = element_text(colour = "#C3A951"),
          axis.title.y.right = element_text(colour = "#5683B0"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

# Generate plots for each centrality measure using function

plot_different_scale(data_between_trans, "Betweenness Centrality (different scale)")
plot_different_scale(data_degree_trans, "Degree Centrality (different scale)")
plot_different_scale(data_eigen_trans, "Eigenvalue Centrality (different scale)")
plot_different_scale(data_page_trans, "Pagerank (different scale)")


#### Robustness_Nie (all centrality measures) ####

# Generate transformed dataframe

transform_Nie <- function(data_trans, centrality) {
  data_trans <- data_trans %>%
    mutate(Centrality = centrality) %>%
  select(p, mean = mean_Nie, std = std_Nie, Centrality)
}

data_Nie <- rbind(transform_Nie(data_between_trans, "Betweenness"),
                  transform_Nie(data_degree_trans, "Degree"),
                  transform_Nie(data_eigen_trans, "Eigenvalue"),
                  transform_Nie(data_page_trans, "Pagerank"))

# Generate plot

ggplot(data_Nie, aes(x = p, y = mean, color = Centrality)) + 
  geom_point() + geom_line() +
  scale_y_continuous("Robustness of Centrality Measure") +
  scale_x_continuous("Proportion of missing nodes p") +
  scale_color_manual(values = c("#C3A951", "#5683B0", "#8C6B88", "#66ACA7")) +
  labs(title = "Robustness Niemeyer (All Centralities)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


#### Robustness_Ros (all centrality measures) ####

# Generate transformed dataframe

transform_Ros <- function(data_trans, centrality) {
  data_trans <- data_trans %>%
    mutate(Centrality = centrality) %>%
    select(p, mean = mean_Ros, std = std_Ros, Centrality)
}

data_Ros <- rbind(transform_Ros(data_between_trans, "Betweenness"),
                  transform_Ros(data_degree_trans, "Degree"),
                  transform_Ros(data_eigen_trans, "Eigenvalue"),
                  transform_Ros(data_page_trans, "Pagerank"))

# Generate plot

ggplot(data_Ros, aes(x = p, y = mean, color = Centrality)) + 
  geom_point() + geom_line() +
  scale_y_continuous("Difference in Outbreak Size vs. Random Immunization") +
  scale_x_continuous("Proportion of missing nodes p") +
  scale_color_manual(values = c("#C3A951", "#5683B0", "#8C6B88", "#66ACA7")) +
  labs(title = "Robustness Rosenblatt (All Centralities)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))




