library(readxl)
library(ggplot2)
library(tidyr)


folder_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_T5_QAOAO/Summary"  # Replace with the actual folder path
files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)

data_list <- lapply(files, read_excel)
# # Empty list to store results
# stats_list <- list()
# 
# # Loop through each Excel file and compute statistics
# for (file in files) {
#   df <- read_excel(file)  # Read the Excel file
#   data_values <- df[[2]]  #2 = CPU time; 3 = Opt Sol
#   
#   # Compute statistics
#   stats_list[[basename(file)]] <- c(
#     mean = mean(data_values, na.rm = TRUE),
#     min = min(data_values, na.rm = TRUE),
#     max = max(data_values, na.rm = TRUE),
#     lq = quantile(data_values, 0.25, na.rm = TRUE),  
#     uq = quantile(data_values, 0.75, na.rm = TRUE),
#     median = median(data_values, na.rm = TRUE)
#   )
# }
# 
# # Convert results to a data frame
# stats_df <- as.data.frame(do.call(rbind, stats_list))
# rownames(stats_df) <- seq_len(nrow(stats_df))
# stats_df <- data.frame(x = seq_len(nrow(stats_df)), stats_df, row.names = NULL)
# 
# # Print the summary table
# print(stats_df)

# Read the Excel file
shots <- data.frame(read_excel('/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_T5_QAOAO/Solutions/p_[20, 20, 5000, 0.5, 0.025]/qaoa_p_6_e_14_bal.xlsx'))
shots$Obj_fuc <- shots$y_1+shots$y_2+shots$y_3
shots$nu_shot <- 1:nrow(shots)


# Create the plot
ggplot(shots, aes(x = nu_shot, y = Obj_fuc)) +
  geom_line(aes(group = 1), color = "black", linewidth = 0.1) +  # Continuous line for all points
  geom_point(data = subset(shots, Type == "TRUE"), aes(color = Type, shape = Type), size = 2) +  # Only "TRUE" points
  scale_color_manual(values = c("TRUE" = "red")) +  # Assign color only to "TRUE"
  scale_shape_manual(values = c("TRUE" = 16)) +  # Assign shape only to "TRUE"
  labs(title = "X-Y Plot with Only 'TRUE' Points Highlighted", x = "X values", y = "Y values") +
  theme_minimal()

# Select only two columns (e.g., ID and Value)
all_shots <- shots[, c("nu_shot", "Obj_fuc")]
write.table(all_shots, "/Users/luismoncayo/Dropbox/Writing/journal_papers/QUBO_SALB/Plots_Tikz/output.csv", sep = ",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8", quote = FALSE)

# Select rows where Type is "True", keeping only ID and Value columns
selected_df <- subset(shots, Type == "TRUE", select = c(nu_shot, Obj_fuc))
write.table(selected_df, "/Users/luismoncayo/Dropbox/Writing/journal_papers/QUBO_SALB/Plots_Tikz/trues.csv", sep = ",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8", quote = FALSE)



