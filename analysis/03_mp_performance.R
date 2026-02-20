#!/usr/bin/env Rscript
######################################
## Manuscript plots: Performance of ML models across configurations
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(EBEx)
library(purrr)
library(dplyr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
# Parameters
input_lists <- c("data_driven", "dea", "mrmr", "disease_related",
                 "disease_related_entire_list", "omnipath_data_driven",
                 "omnipath_disease_related", "omnipath_intersection",
                 "guildify_disease_related", "guildify_data_driven",
                 "omnipath_union")
input_lists <- c("candidate")
pname <- "candidate"
# Define minmax and normMCC for radar charts
minmax <- c(0.5,0.9,0.1)
metric <- "normMCC"
# input_lists <- c("mrmr_30", "mrmr_76",
#                  "disease_related", "data_driven")
# pname <- "mrmr_comparison"
ml_models_to_run_vector <- c("rf", "svm_r", "svm_p", "glm", "knn", "xgb")
# Relative paths
models_results_dir <- "/home/iria/bsc008817/COPD/COPD/COPD"
models_results_dir <- "test"
expression_dir <- "test/analysis/preprocessing"
features_lists_dir <- "test/feature_selection"
output_dir <- file.path("test", "plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. LOAD DATA
print_message("Loading results and SHAP files...")
results_list <- purrr::map(input_lists, function(name) {
  readRDS(file.path(models_results_dir, paste0("results_", name), "models_results.rds"))
}) %>% set_names(input_lists)
# Load expression data
print_message("Loading expression training data...")
load(file.path(expression_dir, "expression_train.Rda"))
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. Process data for plotting
# 2.1 Process Test metrics for ALL configurations
print_message("Processing all Cross-Validation results...")
all_configs_cv_df <- purrr::imap_dfr(results_list, function(models_results, input_name) {

  df <- process_cross_validation_metrics(
    results_models = models_results,
    classifiers = ml_models_to_run_vector
  )
  # Add the column identifying which gene list this is
  df$input_list <- input_name

  return(df)
})
# 2.2 Process Test metrics for ALL configurations
print_message("Processing all Test results...")
all_configs_test_df <- purrr::imap_dfr(results_list, function(models_results, input_name) {

  df <- process_test_metrics(
    results_models = models_results,
    classifiers = ml_models_to_run_vector
  )
  # Add the column identifying which gene list this is
  df$input_list <- input_name

  return(df)
})

# 2.3 Clean names for plotting (using your package functions)
print_message("Cleaning names for plotting...")
all_configs_cv_df <- all_configs_cv_df %>%
  mutate(
    classifier = rename_classifier(classifier),
    input_list = sapply(input_list, rename_input_list) # Optional: if you want pretty list names
  )
all_configs_test_df <- all_configs_test_df %>%
  mutate(
    classifier = rename_classifier(classifier),
    input_list = sapply(input_list, rename_input_list)
  )
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
print_message("Preparing data for barplot of input lists...")
fsm_genes <- data.frame(
  input_list = character(),
  intersection_count = numeric(),
  stringsAsFactors = FALSE
)
# Loop through input lists and calculate intersections
for (genes_list in input_lists) {
  print_message(genes_list)
  # Load the gene list for the current input
  current_genes <- load_gene_list(genes_list, features_lists_dir)
  # Calculate intersection
  intersection <- length(intersect(colnames(expression_train), current_genes))
  print_message("Intersection:", intersection)
  # Append results to fsm_genes
  fsm_genes <- rbind(fsm_genes, data.frame(input_list = genes_list, intersection_count = intersection))
  print_message("---")
}
print_message("Cleaning names for plotting...")
# Clean names for plotting
fsm_genes <- fsm_genes %>%
  mutate(
    input_list = sapply(input_list, rename_input_list)
  )
print_message("Plotting barplot of input lists...")
# Plot
barplot_feature_selection(fsm_genes, paste0(output_dir, "/input_lists_barplot.pdf"))
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. Plot summary ML statistics for all configurations
# 4.1 Cross-validation
# Prepare data for plotting
print_message("Preparing data for summary barplot of classification performance results...")
summary_stats <- prepare_summary_barplot_data(data = all_configs_cv_df, target_metric = "normMCC", value_col = "mean")
# Plot
p_summary <- plot_summary_barplot(summary_data = summary_stats, metric_name = "normMCC", colors = ml_models_colors, base_size = 20)
# Save the plot
ggplot2::ggsave(filename = paste0(output_dir, "/summary_performance_barplot_cv.pdf"), plot = p_summary,width = 9, height = 11)
# 4.2. Test
# Prepare data for plotting
summary_stats <- prepare_summary_barplot_data(data = all_configs_test_df, target_metric = "normMCC", value_col = "estimate")
# Plot
p_summary <- plot_summary_barplot(summary_data = summary_stats, metric_name = "normMCC", colors = ml_models_colors, base_size = 20)
# Save the plot
ggplot2::ggsave(filename = paste0(output_dir, "/summary_performance_barplot_test.pdf"), plot = p_summary,width = 9, height = 11)
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. FINISH
print_message("All plots saved in:", output_dir)
