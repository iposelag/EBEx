#!/usr/bin/env Rscript
######################################
## Manuscript plots: Performance of ML models across configurations
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(devtools)
devtools::load_all()
devtools::document()
library(MyPipelinePkg)
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
ml_models_to_run_vector <- c("rf", "svm_r", "svm_p", "glm", "knn", "xgb")
# Relative paths
models_results_dir <- "/home/iria/bsc008817/COPD/COPD/COPD"
expression_dir <- "test/analysis/preprocessing"
features_lists_dir <- "test/feature_selection"
output_dir <- file.path("test", "plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. LOAD DATA
cat("Loading results and SHAP files...\n")
results_list <- purrr::map(input_lists, function(name) {
  readRDS(file.path(models_results_dir, paste0("results_", name), "results_models.rds"))
}) %>% set_names(input_lists)
# Load expression data
cat("Loading expression training data...\n")
load(file.path(expression_dir, "expression_train.Rda"))
cat("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. Process data for plotting
# 2.1 Process Test metrics for ALL configurations
cat("Processing all Cross-Validation results...\n")
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
cat("Processing all Test results...\n")
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
cat("Cleaning names for plotting...\n")
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
cat("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. Plot classification performance
# Define minmax and normMCC for radar charts
minmax <- c(0.5,0.9,0.1)
metric <- "normMCC"
# Plot Test
cat("Plotting classification performance radar charts...\n")
data_for_plot <- prepare_data_for_radarchart(all_configs_test_df, metric, "estimate", minmax[1], minmax[2])
plot_radarchart(data_for_plot, metric, plot_name = "radarchart_test_candidate", output_dir, minmax[1], minmax[2], minmax[3])
# Cross validation
data_for_plot <- prepare_data_for_radarchart(all_configs_cv_df, metric, "mean", minmax[1], minmax[2])
plot_radarchart(data_for_plot, metric, plot_name = "radarchart_cross_validation_candidate", output_dir, minmax[1], minmax[2], minmax[3])
cat("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. Barplot of input lists
cat("Preparing data for barplot of input lists...\n")
fsm_genes <- data.frame(
  input_list = character(),
  intersection_count = numeric(),
  stringsAsFactors = FALSE
)
# Loop through input lists and calculate intersections
for (genes_list in input_lists) {
  cat(genes_list, "\n")
  # Load the gene list for the current input
  current_genes <- load_gene_list(genes_list, features_lists_dir)
  # Calculate intersection
  intersection <- length(intersect(colnames(expression_train), current_genes))
  cat("Intersection:", intersection, "\n")
  # Append results to fsm_genes
  fsm_genes <- rbind(fsm_genes, data.frame(input_list = genes_list, intersection_count = intersection))
  cat("---\n")
}
cat("Cleaning names for plotting...\n")
# Clean names for plotting
fsm_genes <- fsm_genes %>%
  mutate(
    input_list = sapply(input_list, rename_input_list)
  )
cat("Plotting barplot of input lists...\n")
# Plot
barplot_feature_selection(fsm_genes, paste0(output_dir, "/input_lists_barplot.pdf"))
cat("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. Plot summary ML statistics for all configurations
# 4.1 Cross-validation
# Prepare data for plotting
cat("Preparing data for summary barplot of classification performance results...\n")
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
cat("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. FINISH
cat("All plots saved in:", output_dir, "\n")
