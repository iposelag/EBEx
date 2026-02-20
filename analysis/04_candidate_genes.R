#!/usr/bin/env Rscript
######################################
## EBEx Pipeline: Ranking of genes based on aggregated SHAP values to obtain candidate genes
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(EBEx)
library(purrr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
input_lists <- c("data_driven", "dea", "mrmr", "disease_related",
                 "disease_related_entire_list", "omnipath_data_driven",
                 "omnipath_disease_related", "omnipath_intersection",
                 "guildify_disease_related", "guildify_data_driven",
                 "omnipath_union")
# Relative paths
directory_to_load <- "/home/iria/bsc008817/COPD/COPD/COPD"
output_dir <- file.path("test", "candidate_genes")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Assign decided thresholds after previous analysis
density_threshold <- 0.010938
mclust_groups <- c(8, 9)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. LOAD DATA
print_message("Loading results and SHAP files...")
results_list <- purrr::map(input_lists, function(name) {
  readRDS(file.path(directory_to_load, paste0("results_", name), "results_models.rds"))
}) %>% purrr::set_names(input_lists)
results_shap <- purrr::map(input_lists, function(name) {
  readRDS(file.path(directory_to_load, paste0("results_", name), "results_shap.rds"))
}) %>% purrr::set_names(input_lists)
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. AGGREGATION
print_message("Running aggregation pipeline...")
aggregated_data <- run_aggregation_pipeline(results_list, results_shap, threshold = 0.75, candidate_results_dir = output_dir)
# Save intermediate aggregation summary
write.csv(aggregated_data$summary,
          file.path(output_dir, "aggregated_shap_summary.csv"), row.names = FALSE)
write.csv(aggregated_data$matrix,
          file.path(output_dir, "aggregated_shap_matrix.csv"), row.names = FALSE)
print_message("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. DENSITY RANKING
print_message("Performing Density analysis...")
density_analysis <- compute_density_analysis(aggregated_data$summary, column = "aggregated_max")
# Save Density Plot
pdf(file.path(output_dir, "ranking_density_diagnostic.pdf"), width = 10, height = 8)
plot_density_diagnostic(density_analysis, density_threshold)
dev.off()
final_genes_density <- get_candidate_genes(density_analysis, threshold = density_threshold)
write.csv(final_genes_density,
          file.path(output_dir, "candidate_genes_density.csv"), row.names = FALSE)
print_message("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. MCLUST RANKING
print_message("Performing Mclust analysis...")
mclust_analysis <- run_mclust_analysis(aggregated_data$summary, column = "aggregated_max", G = 1:9)
# Save Mclust Plot
pdf(file.path(output_dir, "ranking_mclust_diagnostic.pdf"), width = 10, height = 12)
plot_mclust_diagnostic(mclust_analysis, sel_groups = mclust_groups)
dev.off()
# Save Mclust summary
final_genes_mclust <- select_genes_mclust(mclust_analysis, sel_groups = mclust_groups)
write.csv(final_genes_mclust,
          file.path(output_dir, "candidate_genes_mclust.csv"), row.names = FALSE)
print_message("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 6. FINISH
print_message("Pipeline finished. Results saved in:", output_dir)
