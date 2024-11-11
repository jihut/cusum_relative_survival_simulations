rm(list = ls())
library(reticulate)

cusum_signal_df_eta_normal <- readRDS("r_example_scripts/simulations/prop_piecewise_nr3/cusum_signal_df_eta_normal.RDS")
cusum_signal_df_eta_instant <- readRDS("r_example_scripts/simulations/prop_piecewise_nr3/cusum_signal_df_eta_instant.RDS")
cusum_signal_df_eta_star <- readRDS("r_example_scripts/simulations/prop_piecewise_nr3/cusum_signal_df_eta_star.RDS")

final_cusum_signal_df <- rbind(cusum_signal_df_eta_instant, cusum_signal_df_eta_normal, cusum_signal_df_eta_star)
final_cusum_signal_df$signal_ratio <- paste(final_cusum_signal_df$signal_ratio, "%", sep = "")

use_virtualenv("/Users/jimmytran/Documents/General Python Workspace/my_env")
final_table <- r_to_py(final_cusum_signal_df)

py$final_table <- final_table
py_run_string("
import pandas as pd; 

pivot_table = final_table.astype({'eta': pd.CategoricalDtype(['eta = 0', 'eta = 2.5', 'eta star = 2.5'], ordered = True)}).pivot(index = ['rho', 'false_prob'], columns = 'eta', values = ['signal_ratio'])

latex = pivot_table.to_latex(
        index=True,
        escape=False,
        sparsify=True,
        multirow=True,
        multicolumn=True,
        multicolumn_format='c'
    )

print(latex)
")
