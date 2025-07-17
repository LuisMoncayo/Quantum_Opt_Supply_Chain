library(readxl)
library(ggplot2)
library(tidyr)
# 
#folder_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L5_Slack/Summary"
folder_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L6_Unbalanced/Summary"
files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)

names_all_exper <- print(files)

data_list <- lapply(files, read_excel)

# Assume data_list contains multiple data frames
df_combined <- data.frame(matrix(ncol = length(data_list), nrow = nrow(data_list[[1]])))
# 
# Extract a specific column (e.g., "Value") from each data frame
for (i in seq_along(data_list)) {
  df_combined[[i]] <- data_list[[i]]$Nu_Opt_Sol  # Change "Value" to the desired column name
}

# Rename columns based on file names (optional)
colnames(df_combined) <- basename(files)  # files = list of file paths
percentages_all <- df_combined

percentages_all[,1] <- (percentages_all[,1])/8000
percentages_all[,2] <- (percentages_all[,2])/8000
percentages_all[,3] <- (percentages_all[,3])/8000
percentages_all[,4] <- (percentages_all[,4])/8000
percentages_all[,5] <- (percentages_all[,5])/8000

pct_per_expe <- data.frame(
  Mean = round(sapply(percentages_all, mean) * 100, digits = 2),
  SD = round(sapply(percentages_all, sd) * 100, digits = 2)
)

#### second part #####
# Extract the list inside the square brackets
param_lists <- regmatches(names_all_exper, regexpr("\\[.*?\\]", names_all_exper))
cat(paste(param_lists, collapse = "\n"))

#base_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L5_Slack/Solutions/p_"
base_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L6_Unbalanced/Solutions/p_"
for (param in param_lists) {
  folder_path <- paste0(base_path, param)#full_path
  
  files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
  data_list <- lapply(files, read_excel)
  
  experiment_1 <- data.frame(data_list[[1]])
  names(experiment_1)[2] <- "Values_Obj_Func"
  experiment_TRUE_1 <- experiment_1[experiment_1$Type_Sol == "TRUE",]
  d_1 <- as.data.frame(table(experiment_TRUE_1$Values_Obj_Func))
  
  merged_df <- d_1
  for (i in 2:length(data_list)) {
    experiment <- data.frame(data_list[[i]])
    names(experiment)[2] <- "Values_Obj_Func"
    experiment_TRUE <- experiment[experiment$Type_Sol == "TRUE",]
    d_i <- as.data.frame(table(experiment_TRUE$Values_Obj_Func))
    if (nrow(d_i) == 0){
      next
    }
    merged_df <- merge(merged_df, d_i, by = "Var1", all = TRUE)
  }
  
  merged_df[is.na(merged_df)] <- 0

  # Step 1: Remove the first column
  df_no_var1 <- merged_df[, -1]
  row_means <- rowMeans(df_no_var1, na.rm = TRUE)

  # Step 2: Compute column sums
  col_sums <- colSums(df_no_var1, na.rm = TRUE)

  # Step 3: Get the first row
  first_row <- df_no_var1[1, ]

  # Step 4: Compute percentages
  first_row_pct <- (first_row / col_sums) * 100

  # View the result
  #print(round(first_row_pct, 2))

  # Get the first row as a numeric vector
  values <- as.numeric(first_row_pct[1, ])

  # Compute mean and standard deviation
  mean_val <- mean(values, na.rm = TRUE)
  sd_val <- sd(values, na.rm = TRUE)

  # Print result
  #cat(sprintf("Mean\tSD\n%.2f\t%.2f\n", mean_val, sd_val))

  ### printings
  print(param)
  cat(sprintf("Mean\tSD\n%.4f\t%.4f\n", mean_val, sd_val))
  print(d_1)
  #cat(full_path, "\n")
  print(row_means)
  print("++++++***********+++++++++++")
}
print(pct_per_expe)
write.table(pct_per_expe, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





