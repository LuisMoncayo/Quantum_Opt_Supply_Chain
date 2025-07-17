library(readxl)
library(ggplot2)
library(tidyr)

#/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_T5_QAOAO/Summary
# folder_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_T5_QAOAO/Summary"
# files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
# 
# data_list <- lapply(files, read_excel)
# 
# 
# # Assume data_list contains multiple data frames
# df_combined <- data.frame(matrix(ncol = length(data_list), nrow = nrow(data_list[[1]])))
# 
# # Extract a specific column (e.g., "Value") from each data frame
# for (i in seq_along(data_list)) {
#   df_combined[[i]] <- data_list[[i]]$Nu_Opt_Sol  # Change "Value" to the desired column name
# }
# 
# # Rename columns based on file names (optional)
# colnames(df_combined) <- basename(files)  # files = list of file paths
# 
# # View the final data frame
# #head(df_combined)
# 
# df_long <- pivot_longer(df_combined, cols = everything(), names_to = "Experiment", values_to = "Nu_Opt_Sol")
# 
# 
# ggplot(df_long, aes(x = Experiment, y = Nu_Opt_Sol, fill = Experiment)) +
#   geom_boxplot(width = 0.2, alpha = 0.7) +
#   geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
#   theme_minimal() +
#   labs(title = "Box Plot with Jittered Points")
# 
# 
#     ################################################################################
# # Extract a specific column (e.g., "Value") from each data frame
# for (i in seq_along(data_list)) {
#   df_combined[[i]] <- data_list[[i]]$CPU_Time  # Change "Value" to the desired column name
# }
# 
# # Rename columns based on file names (optional)
# colnames(df_combined) <- basename(files)  # files = list of file paths
# 
# # View the final data frame
# #head(df_combined)
# 
# df_long <- pivot_longer(df_combined, cols = everything(), names_to = "Experiment", values_to = "CPU_Time")
# 
# 
# ggplot(df_long, aes(x = Experiment, y = CPU_Time, fill = Experiment)) +
#   geom_boxplot(width = 0.2, alpha = 0.7) +
#   geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
#   theme_minimal() +
#   labs(title = "Box Plot with Jittered Points")
################################################################################
################################################################################
# Single file
setwd("/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_T5_QAOAO/Solutions/p_[20, 20, 5000, 0.5, 0.025]/")
single <- 'qaoa_p_6_e_1_bal.xlsx'
df <- read_excel(single, col_names = TRUE)
df$Obj_Fun <- rowSums(df[, c("y_1", "y_2", "y_3")], na.rm = TRUE)
df$Shot <- seq_len(nrow(df))
 
ggplot(df, aes(x = Shot, y = Obj_Fun, color = Type)) +
  geom_point(size = 0.5) +  # Increase point size
  labs(title = "Scatter Plot of ID vs Value", x = "ID", y = "Value") +
  theme_minimal()
