#' Replot DR Power - Exact Multi-Threshold Analysis from Raw p-values
#' Uses consolidated plotting functions from code/plot_dr_power.R

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/plot_dr_power.R")

# Generate DR power plot (convenience wrapper)
replotDRPower(test_name = "DR",
              results_dir = "output/dr_power_stratified",
              output_file = "output/figures/dr_power_stratified.pdf")
