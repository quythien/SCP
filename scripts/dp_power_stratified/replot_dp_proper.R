#' Replot DP Power - Exact Multi-Threshold Analysis from Raw p-values
#' Uses consolidated plotting functions from code/plot_dr_power.R

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/plot_dr_power.R")

# Generate DP power plot (convenience wrapper)
replotDRPower(test_name = "DP",
              results_dir = "output/dp_power_stratified",
              output_file = "output/figures/dp_power_stratified.pdf")
