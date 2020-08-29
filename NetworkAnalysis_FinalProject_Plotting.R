#### Presettings ####

# Load required packages

library(ggplot2)
library(tidyverse)
library(dplyr)


# Import dataframes containing results from outbreak simulations

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
# transform_data

data_between_trans <- transform_data(data_between)
data_degree_trans <- transform_data(data_degree)
data_eigen_trans <- transform_data(data_eigen)
data_page_trans <- transform_data(data_page)


#### Comparison of the original measures ####

# Define function to plot the measures proposed by Martin and Niemeyer (2019)
# and Rosenblatt et al. (2020)

plot <- function(data_trans, title) {
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
    scale_y_continuous("Difference in Outbreak Size vs. Random Immunization", 
                       sec.axis = sec_axis(~.*1, name = "Robustness of Centrality Measure")) +
    scale_x_continuous("Proportion of missing nodes p") +
    scale_color_manual(values = c("red", "blue")) +
    labs(title = title) +
    theme_bw() +
    theme(axis.title.y.left = element_text(colour = "blue", size = 16),
          axis.title.y.right = element_text(colour = "red", size = 16),
          axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
}

# Generate plots for each centrality measure using function plot

plot(data_between_trans, "Betweenness Centrality")
plot(data_degree_trans, "Degree Centrality")
plot(data_eigen_trans, "Eigenvector Centrality")
plot(data_page_trans, "Pagerank")


#### Traditional Robustness ####

# Define function to transform data 

transform_Nie <- function(data_trans, centrality) {
  data_trans <- data_trans %>%
    mutate(Centrality = centrality) %>%
    select(p, mean = mean_Nie, std = std_Nie, Centrality)
}

# Generate transformed dataframe using function transform_Nie

data_Nie <- rbind(transform_Nie(data_between_trans, "Betweenness"),
                  transform_Nie(data_degree_trans, "Degree"),
                  transform_Nie(data_eigen_trans, "Eigenvector"),
                  transform_Nie(data_page_trans, "Pagerank"))

# Generate plot with traditional robustness for all centrality measures

ggplot(data_Nie, aes(x = p, y = mean, color = Centrality)) + 
  geom_point() + geom_line() +
  scale_y_continuous("Robustness of Centrality Measure") +
  scale_x_continuous("Proportion of missing nodes p") +
  scale_color_manual(values = c("red", "blue", "green", "yellow3")) +
  labs(title = "Traditional robustness") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = rel(0)),
        legend.text = element_text(size = 14))


#### Immunization Effectiveness & Immunization Robustness ####

# Define function to transform data 

transform_Ros <- function(data_trans, centrality) {
  init <- data_trans$mean_Ros[1]
  data_trans <- data_trans %>%
    mutate(Centrality = centrality,
           diff_abs = mean_Ros - init) %>%
    select(p, mean = mean_Ros, std = std_Ros, Centrality, diff_abs)
}

# Generate transformed dataframe using function transform_Ros

data_Ros <- rbind(transform_Ros(data_between_trans, "Betweenness"),
                  transform_Ros(data_degree_trans, "Degree"),
                  transform_Ros(data_eigen_trans, "Eigenvector"),
                  transform_Ros(data_page_trans, "Pagerank"))

# Generate plot with immunization effectiveness for all centrality measures

ggplot(data_Ros, aes(x = p, y = mean, color = Centrality)) + 
  geom_point() + geom_line() +
  scale_y_continuous("Difference in Outbreak Size vs. Random Immunization") +
  scale_x_continuous("Proportion of missing nodes p") +
  scale_color_manual(values = c("red", "blue", "green", "yellow3")) +
  labs(title = "Immunization Effectiveness") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = rel(0)),
        legend.text = element_text(size = 14))

# Generate plot with immunization robustness for all centrality measures

ggplot(data_Ros, aes(x = p, y = diff_abs, color = Centrality)) + 
  geom_point() + geom_line() +
  scale_y_continuous("Difference in Effectiveness to error-free network") +
  scale_x_continuous("Proportion of missing nodes p") +
  scale_color_manual(values = c("red", "blue", "green", "yellow3")) +
  labs(title = "Immunization Robustness") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = rel(0)),
        legend.text = element_text(size = 14))
