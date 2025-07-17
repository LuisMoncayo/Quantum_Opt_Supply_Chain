library(readxl)
library(ggplot2)
library(tidyr)
# /Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_LB_T5_Unbalanced/Summary
folder_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L5_Unbalanced/Summary"
files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)

data_list <- lapply(files, read_excel)

df_combined <- data.frame(matrix(ncol = length(data_list), nrow = nrow(data_list[[1]])))

for (i in seq_along(data_list)) {
  df_combined[[i]] <- data_list[[i]]$Nu_Opt_Sol
}

colnames(df_combined) <- basename(files)
df_long <- pivot_longer(df_combined, cols = everything(), names_to = "Experiment", values_to = "Nu_Opt_Sol")

ggplot(df_long, aes(x = Experiment, y = Nu_Opt_Sol, fill = Experiment)) +
  geom_boxplot(width = 0.2, alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  theme_minimal() +
  labs(title = "Box Plot with Jittered Points")

#Extract a specific column (e.g., "Value") from each data frame
for (i in seq_along(data_list)) {
  df_combined[[i]] <- data_list[[i]]$Nu_Opt_Sol  # Change "Value" to the desired column name Nu_Opt_Sol, CPU_Time
  avg <- mean(data_list[[i]]$CPU_Time)   # average
  std_dev <- sd(data_list[[i]]$CPU_Time)     # standard deviation
  min_val <- min(data_list[[i]]$CPU_Time)
  max_val <- max(data_list[[i]]$CPU_Time)
  cat("Mean:", round(avg,2), "SD:", round(std_dev,2), "Min:", round(min_val,2), "Max:", round(max_val,2), "\n")
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
write.table(result, file = "/Users/luismoncayo/Downloads/summary_cpu.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
print(result)