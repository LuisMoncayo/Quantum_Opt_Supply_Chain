# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
import pandas as pd
from glob import glob
import re
import os
import matplotlib.pyplot as plt
import numpy as np

### Files Name in Folder 
# folder = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Results/LB_Unba_Noisy_IBM/LB_5/"
# files = os.listdir(folder)
# print("All files:")
# for f in files:
#     print(f)
##################################################


##################
# -> FL Slack IBM
##################
# folder = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Results/FL_Slack_Noisy_IBM/FL_7"
# files = sorted(glob(f"{folder}/facility_location_results_noisy_run*.csv"))

# rows = []
# for f in files:
#     df = pd.read_csv(f)
#     run_id = int(f.split("run")[-1].split(".csv")[0])
#     df_feas = df[df["feasible"].astype(str).str.lower() == "feasible"]
#     sum_feas = df_feas["counts"].sum()
#     sum_opt_feas = df_feas[df_feas["energy_M"] == 2]["counts"].sum()
#     rows.append({"run": run_id, "total_feasible": sum_feas, "optimal_feasible": sum_opt_feas})

# summary = pd.DataFrame(rows).sort_values("run")
# #summary.to_csv(f"{folder}/summary_feasible_counts_weighted.csv", index=False)

# # Replace with your full log file path
# path = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Results/FL_Slack_Noisy_IBM/FL_7/summary.txt"

# with open(path) as f:
#     log = f.read()

# elapsed = [float(x) for x in re.findall(r"Elapsed time: ([\d\.]+) seconds", log)]

# elapsed_df = pd.DataFrame({
#     "run": range(1, len(elapsed) + 1),
#     "elapsed_time_sec": elapsed
# })

# all_data = pd.merge(summary, elapsed_df, on="run", how="inner")
# print(all_data)

# # Box Plot
# plt.figure(figsize=(6, 4))
# plt.boxplot(all_data["total_feasible"], patch_artist=True,
#             boxprops=dict(facecolor="lightblue", color="black"),
#             medianprops=dict(color="red", linewidth=2))
# plt.ylabel("Optimal Feasible")
# plt.title("Distribution of Optimal Feasible Solutions per Run")
# plt.grid(True, axis="y", linestyle="--", alpha=0.6)
# plt.tight_layout()
# plt.show()

# # Example data column
# # all_data = pd.DataFrame({"optimal_feasible": [...]})

# data = all_data["total_feasible"].dropna()

# # Compute statistics
# average = np.mean(data)
# lower_whisker = np.min(data)
# upper_whisker = np.max(data)
# lower_quartile = np.percentile(data, 25)
# upper_quartile = np.percentile(data, 75)
# median = np.percentile(data, 50)

# # Choose plotting position and label
# position = 16
# label = "P4C1"

# # Build LaTeX line in your preferred order
# latex_line = (
#     f"\\addplot[boxplot style,boxplot prepared colored="
#     f"{{{average:.2f}/{lower_whisker:.2f}/{upper_whisker:.2f}/"
#     f"{lower_quartile:.2f}/{upper_quartile:.2f}/{median:.2f} "
#     f"at {position} with {label}}}] coordinates {{}};"
# )

# print(latex_line)

# avg = all_data["elapsed_time_sec"].mean()
# std = all_data["elapsed_time_sec"].std()

# print(f"Average elapsed time: {avg:.3f} sec")
# print(f"Standard deviation:  {std:.3f} sec")
##############################################################################
# -> FL Unbalanced IBM
##################
# folder = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Results/FL_Unba_Noisy_IBM/FL_8"
# files = sorted(glob(f"{folder}/facility_location_results_noisy_run*.csv"))

# rows = []
# for f in files:
#     df = pd.read_csv(f)
#     run_id = int(f.split("run")[-1].split(".csv")[0])
#     df_feas = df[df["feasible"].astype(str).str.lower() == "feasible"]
#     sum_feas = df_feas["counts"].sum()
#     sum_opt_feas = df_feas[df_feas["energy_M"] == 2]["counts"].sum()
#     rows.append({"run": run_id, "total_feasible": sum_feas, "optimal_feasible": sum_opt_feas})

# summary = pd.DataFrame(rows).sort_values("run")
# #summary.to_csv(f"{folder}/summary_feasible_counts_weighted.csv", index=False)

# # Replace with your full log file path
# path = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Results/FL_Unba_Noisy_IBM/FL_8/summary.txt"

# with open(path) as f:
#     log = f.read()

# elapsed = [float(x) for x in re.findall(r"Elapsed time: ([\d\.]+) seconds", log)]

# elapsed_df = pd.DataFrame({
#     "run": range(1, len(elapsed) + 1),
#     "elapsed_time_sec": elapsed
# })

# all_data = pd.merge(summary, elapsed_df, on="run", how="inner")
# print(all_data)

# plt.figure(figsize=(6, 4))
# plt.boxplot(all_data["total_feasible"], patch_artist=True,
#             boxprops=dict(facecolor="lightblue", color="black"),
#             medianprops=dict(color="red", linewidth=2))
# plt.ylabel("Optimal Feasible")
# plt.title("Distribution of Optimal Feasible Solutions per Run")
# plt.grid(True, axis="y", linestyle="--", alpha=0.6)
# plt.tight_layout()
# plt.show()

# # Example data column
# # all_data = pd.DataFrame({"optimal_feasible": [...]})

# data = all_data["total_feasible"].dropna()

# # Compute statistics
# average = np.mean(data)
# lower_whisker = np.min(data)
# upper_whisker = np.max(data)
# lower_quartile = np.percentile(data, 25)
# upper_quartile = np.percentile(data, 75)
# median = np.percentile(data, 50)

# # Choose plotting position and label
# position = 16
# label = "P4C1"

# # Build LaTeX line in your preferred order
# latex_line = (
#     f"\\addplot[boxplot style,boxplot prepared colored="
#     f"{{{average:.2f}/{lower_whisker:.2f}/{upper_whisker:.2f}/"
#     f"{lower_quartile:.2f}/{upper_quartile:.2f}/{median:.2f} "
#     f"at {position} with {label}}}] coordinates {{}};"
# )

# print(latex_line)

# avg = all_data["elapsed_time_sec"].mean()
# std = all_data["elapsed_time_sec"].std()

# print(f"Average elapsed time: {avg:.3f} sec")
# print(f"Standard deviation:  {std:.3f} sec")

###############################################################################

###############################################################################
# LB Unbalanced
###############################################################################

#Folder path
folder = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Results/LB_Unba_Noisy_IBM/LB_8"

# Collect all files that match pattern, skip .DS_Store and summary.txt
files = sorted([
    f for f in glob(f"{folder}/LB5_QAOA_Run*_Shots_*.csv")
    if os.path.isfile(f) and not f.endswith((".DS_Store", "summary.txt"))
])

rows = []
for f in files:
    df = pd.read_csv(f)
    run_id = int(f.split("Run")[1].split("_")[0])
    feas = df[df["feasible"] == True]
    sol_opt = feas["shot_id"].sum()
    rows.append({"run": run_id, "sol_opt": sol_opt})

summary = pd.DataFrame(rows).sort_values("run")

# out_path = f"{folder}/LB5_QAOA_summary_sol_opt.csv"
# summary.to_csv(out_path, index=False)
print(summary)


# Replace with your actual log file path
path = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Results/LB_Unba_Noisy_IBM/LB_8/summary.txt"

with open(path) as f:
    text = f.read()

# Extract run and elapsed time
runs = [int(x) for x in re.findall(r"RUN (\d+)", text)]
elapsed = [float(x) for x in re.findall(r"Elapsed: ([\d.]+) sec", text)]

# Build dataframe
df = pd.DataFrame({"run": runs, "elapsed_time_sec": elapsed}).sort_values("run")

# Merge both dataframes on 'run'
merged = pd.merge(summary, df, on="run", how="inner").sort_values("run")

print(merged)


plt.figure(figsize=(6, 4))
plt.boxplot(merged["sol_opt"], patch_artist=True,
            boxprops=dict(facecolor="lightblue", color="black"),
            medianprops=dict(color="red", linewidth=2))
plt.ylabel("Optimal Feasible")
plt.title("Distribution of Optimal Feasible Solutions per Run")
plt.grid(True, axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()

# Example data column
# all_data = pd.DataFrame({"optimal_feasible": [...]})

data = merged["sol_opt"].dropna()

# Compute statistics
average = np.mean(data)
lower_whisker = np.min(data)
upper_whisker = np.max(data)
lower_quartile = np.percentile(data, 25)
upper_quartile = np.percentile(data, 75)
median = np.percentile(data, 50)

# Choose plotting position and label
position = 16
label = "P4C1"

# Build LaTeX line in your preferred order
latex_line = (
    f"\\addplot[boxplot style,boxplot prepared colored="
    f"{{{average:.2f}/{lower_whisker:.2f}/{upper_whisker:.2f}/"
    f"{lower_quartile:.2f}/{upper_quartile:.2f}/{median:.2f} "
    f"at {position} with {label}}}] coordinates {{}};"
)

print(latex_line)

avg = merged["elapsed_time_sec"].mean()
std = merged["elapsed_time_sec"].std()

print(f"Average elapsed time: {avg:.3f} sec")
print(f"Standard deviation:  {std:.3f} sec")









