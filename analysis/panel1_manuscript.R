#!/usr/bin/env Rscript
######################################
## Manuscript Figure 2: Performance of ML models across input gene lists
## Regenerates the normMCC radar charts (test and cross-validation) and the
## barplot of input-list sizes from the saved source data in test/manuscript_data/.
## Same plotting logic as sections 3 and 4 of analysis/03_mp_performance.R
## (11 input lists, pname = "performance").
## Iria Pose
## 2026
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
devtools::load_all(".")   # o library(EBEx)
library(dplyr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
pname <- "performance"
# Define minmax and normMCC for radar charts
minmax <- c(0.5,0.9,0.1)
metric <- "normMCC"
# Relative paths
data_dir <- "test/manuscript_data"
output_dir <- file.path("test", "panels_manuscript", "panel1")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. LOAD DATA
# Processed metrics (output of section 2 of 03_mp_performance.R)
print_message("Loading source data...")
all_configs_test_df <- readRDS(file.path(data_dir, "fig2_all_configs_test_df.rds"))
all_configs_cv_df <- readRDS(file.path(data_dir, "fig2_all_configs_cv_df.rds"))
# Number of genes of each input list present in the expression matrix
fsm_genes <- readRDS(file.path(data_dir, "fig2_input_lists_gene_counts.rds"))
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. Plot classification performance
# Plot Test
print_message("Plotting classification performance radar charts...")
data_for_plot <- prepare_data_for_radarchart(all_configs_test_df, metric, "estimate", minmax[1], minmax[2])
plot_radarchart(data_for_plot, metric, plot_name = paste0("radarchart_test_", pname), output_dir, minmax[1], minmax[2], minmax[3])
# Cross validation
data_for_plot <- prepare_data_for_radarchart(all_configs_cv_df, metric, "mean", minmax[1], minmax[2])
plot_radarchart(data_for_plot, metric, plot_name = paste0("radarchart_cross_validation_",pname), output_dir, minmax[1], minmax[2], minmax[3])
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. Barplot of input lists
print_message("Plotting barplot of input lists...")
barplot_feature_selection(fsm_genes, paste0(output_dir, "/input_lists_barplot.pdf"))
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. FINISH
print_message("All plots saved in:", output_dir)
