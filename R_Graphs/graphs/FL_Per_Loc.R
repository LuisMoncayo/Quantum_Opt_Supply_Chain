library(readxl)
library(ggplot2)
library(tidyr)
# 
folder_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L5_Slack/Summary"
#folder_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L5_Unbalanced/Summary"
files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)

data_list <- lapply(files, read_excel)

# Assume data_list contains multiple data frames
df_combined <- data.frame(matrix(ncol = length(data_list), nrow = nrow(data_list[[1]])))
# 
# # Extract a specific column (e.g., "Value") from each data frame
for (i in seq_along(data_list)) {
  df_combined[[i]] <- data_list[[i]]$Nu_Opt_Sol  # Change "Value" to the desired column name
}

# Rename columns based on file names (optional)
colnames(df_combined) <- basename(files)  # files = list of file paths

# View the final data frame
#head(df_combined)

df_long <- pivot_longer(df_combined, cols = everything(), names_to = "Experiment", values_to = "Nu_Opt_Sol")


ggplot(df_long, aes(x = Experiment, y = Nu_Opt_Sol, fill = Experiment)) +
  geom_boxplot(width = 0.2, alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  theme_minimal() +
  labs(title = "Box Plot with Jittered Points")
# 
# 
# ################################################################################
# 
#Extract a specific column (e.g., "Value") from each data frame
for (i in seq_along(data_list)) {
  df_combined[[i]] <- data_list[[i]]$Nu_Opt_Sol  # Change "Value" to the desired column name
  avg <- mean(data_list[[i]]$Nu_Opt_Sol)   # average
  std_dev <- sd(data_list[[i]]$Nu_Opt_Sol)     # standard deviation
  min_val <- min(data_list[[i]]$Nu_Opt_Sol)
  max_val <- max(data_list[[i]]$Nu_Opt_Sol)
  cat("Mean:", round(avg,2), "SD:", round(std_dev,2), "Min:", round(min_val,2), "Max:", round(max_val,2), "\n")
}

# Rename columns based on file names (optional)
colnames(df_combined) <- basename(files)  # files = list of file paths

# View the final data frame
#head(df_combined)

df_long <- pivot_longer(df_combined, cols = everything(), names_to = "Experiment", values_to = "CPU_Time")

ggplot(df_long, aes(x = Experiment, y = CPU_Time, fill = Experiment)) +
  geom_boxplot(width = 0.2, alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  theme_minimal() +
  labs(title = "Box Plot with Jittered Points")

##------------------------------------------------------------------------------
cols <- lapply(data_list, function(x) x$Nu_Opt_Sol)
#cols <- lapply(data_list, function(x) x$CPU_Time)
result_matrix <- do.call(cbind, cols)
colnames(result_matrix) <- paste0("s", 1:length(data_list))

summary_table <- apply(result_matrix, 2, function(x) {
  c(
    mean = mean(x),
    min = min(x),
    max = max(x),
    lq = quantile(x, 0.25),
    uq = quantile(x, 0.75),
    median = median(x)
  )
})

result <- data.frame(
  s = 1:length(data_list),
  round(t(summary_table), 4)
)
print(result)

# Save the summary table to a text file
write.table(result, file = "/Users/luismoncayo/Downloads/summary_output_1.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
print(result)
# 
# ##------------------------------------------------------------------------------


# ## The box plot



# ## The box plot per open facilities -----------

folder_path <- '/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L5_Slack/Solutions/p_[20, 20, 5000, 0.5]'
#folder_path <- '/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L5_Unbalanced/Solutions/p_[30, 30, 8000, 0.5, 0.01]'

files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
data_list <- lapply(files, read_excel)

experiment_1 <- data.frame(data_list[[1]])
names(experiment_1)[2] <- "Values_Obj_Func"
experiment_TRUE_1 <- experiment_1[experiment_1$Type_Sol == "TRUE",]
d_1 <- as.data.frame(table(experiment_TRUE_1$Values_Obj_Func))

# 
# merged_df <- d_1
# for (i in 2:length(data_list)) {
#   experiment <- data.frame(data_list[[i]])
#   names(experiment)[2] <- "Values_Obj_Func"
#   experiment_TRUE <- experiment[experiment$Type_Sol == "TRUE",]
#   d_i <- as.data.frame(table(experiment_TRUE$Values_Obj_Func))
#   merged_df <- merge(merged_df, d_i, by = "Var1", all = TRUE)
# }
# merged_df[is.na(merged_df)] <- 0
# 
# # Get the integer values from the first column as characters
# col_names <- as.character(merged_df[[1]])
# 
# # Transpose the rest of the data
# df_t <- t(merged_df[ , -1])
# 
# # Convert to data frame and set column names
# df_t <- as.data.frame(df_t)
# colnames(df_t) <- col_names
# rownames(df_t) <- NULL
# 
# # Convert wide to long format
# df_long <- pivot_longer(df_t, cols = everything(), names_to = "Variable", values_to = "Value")
# 
# # Plot with ggplot
# ggplot(df_long, aes(x = Variable, y = Value)) +
#   geom_boxplot(fill = "skyblue") +
#   labs(title = basename(folder_path), x = "Open Facilities", y = "Frequency") +
#   theme_minimal()
# 
# copied_df <- df_long
# copied_df[[1]] <- as.integer(copied_df[[1]])
# 
# unique_values <- unique(copied_df$Variable )
# print(unique_values)
# 
# empty_df <- data.frame(matrix(nrow = 0, ncol = 0))
# for (num in unique_values) {
#   subset_df <- copied_df[copied_df$Variable == num, ]
#   x <- subset_df$Value
#   summary_c <- c(
#     mean = mean(x),
#     min = min(x),
#     max = max(x),
#     lq = quantile(x, 0.25),
#     uq = quantile(x, 0.75),
#     median = median(x)
#   )
#   empty_df <- rbind(empty_df,summary_c)
#   #print(summary_c)
# }
# colnames(empty_df) <- c("mean","min","max","lq","uq","median")
# 
# print(round(empty_df,0))










