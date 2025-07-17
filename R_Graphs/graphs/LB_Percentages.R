library(readxl)
library(ggplot2)
library(tidyr)
#

folder_path <- "/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_LB_T9_Unbalanced/Summary"
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

param_lists <- regmatches(names_all_exper, regexpr("\\[.*?\\]", names_all_exper))
cat(paste(param_lists, collapse = "\n"))

for ( i in seq_along(param_lists)){
  parts <- gsub("\\[|\\]", "", param_lists[i]) |> strsplit(", ") |> unlist()
  n_shots <- as.integer(parts[3])
  
  per_para <- (df_combined[[i]]/n_shots)*100
  
  mean_val_all <- mean(per_para)
  sd_val_all <- sd(per_para)
  
  print(param_lists[i])
  print(per_para)
  print(n_shots)
  cat(sprintf("Mean_all\tSD_all\n%.2f\t%.2f\n", mean_val_all, sd_val_all))
  print("---------------------------------------")
}

cat(paste(param_lists, collapse = "\n"))










