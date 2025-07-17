library(dplyr)
library(readxl)


path_to_file <- '/Users/luismoncayo/Dropbox/Python/Balancing_QUBO_U/R_Graphs/Outputs_Balancing_FL_L5_Slack/Solutions/p_[10, 10, 3000, 1]/p_[10, 10, 3000, 1]_e_1_FL.xlsx'
df <- data.frame(read_excel(path_to_file))

# Count how many times each row appears
counts <- as.data.frame(table(df))
# Filter to keep only rows that actually appear
counts <- counts[counts$Freq > 0, ]
counts_sorted <- counts[order(counts$Freq, decreasing = TRUE), ]
counts_sorted <- counts_sorted %>% mutate(Sol = row_number())

ggplot(counts_sorted, aes(x = Sol, y = Freq)) +
  geom_point(color = "lightgray", size = 1) +  # All points in light gray
  geom_point(data = subset(counts_sorted, Type_Sol == TRUE), aes(x = Sol, y = Freq), color = "red", size = 1.5) +  # Only TRUE points
  labs(x = "Solution", y = "Counts") +
  theme_minimal()

df_TRUES <- counts[counts$Type_Sol == "TRUE", ]
df_TRUES <- df_TRUES %>% mutate(Sol = row_number())

ggplot(df_TRUES, aes(x = Sol, y = Freq)) +
  #geom_line(color = "blue") +
  geom_point(color = "red", size = 1) +
  labs(title = "X-Y Line and Points", x = "Energy", y = "Frequency") +
  theme_minimal()








