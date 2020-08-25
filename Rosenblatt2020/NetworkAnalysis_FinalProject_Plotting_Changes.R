data <- data_between_trans

init <- data$mean_Ros[1]
prev <- data$mean_Ros[-9]

# Show difference to (true) network (p = 0)
data_init <- data %>%
  mutate(diff_abs = mean_Ros - init,
         diff_pct = (mean_Ros - init) / init) 

ggplot(data_init, aes(x = p, y = diff_abs)) + 
  geom_point() + 
  geom_line() +
  scale_y_continuous("Abs. Difference to network with error-level p = 0") +
  scale_x_continuous("Error level p") 

ggplot(data_init, aes(x = p, y = diff_pct)) + 
  geom_point() + 
  geom_line() +
  scale_y_continuous("Percentage-Difference to network with error-level p = 0") +
  scale_x_continuous("Error level p") 

# Show difference to previous
data_prev <- data[-1,] %>%
  mutate(diff_abs = mean_Ros - prev,
         diff_pct = (mean_Ros - prev) / prev) 

ggplot(data_prev, aes(x = p, y = diff_abs)) + 
  geom_point() + 
  geom_line() +
  scale_y_continuous("Abs. Difference to previous error-level network") +
  scale_x_continuous("Error level p") 

ggplot(data_prev, aes(x = p, y = diff_pct)) + 
  geom_point() + 
  geom_line() +
  scale_y_continuous("Percentage-Difference to previous error-level networ") +
  scale_x_continuous("Error level p")
