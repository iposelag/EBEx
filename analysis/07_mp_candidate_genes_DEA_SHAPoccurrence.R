#!/usr/bin/env Rscript
######################################
## Analysis of candidate genes in DEA and SHAP vs Occurrence plot
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
devtools::load_all()
devtools::document()
library(MyPipelinePkg)
library(dplyr); library(ggplot2); library(ggrepel)

## --------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
p_val_thresh <- 0.01
lfc_thresh <- log2(1.5)
dea_results_dir <- "test/analysis/DEA"
candidate_results_dir <- "test/candidate_genes"
output_dir <- "test/plots"

## --------------------------------------------------------------------------------------------------------------------
# 2. LOAD DATA
cat("Loading data...\n")
deg_results <- read.csv(file.path(dea_results_dir, "deg_results.csv"))
aggregated_shap_data <- read.csv(file.path(candidate_results_dir, "aggregated_shap_summary.csv"))
canonical_genes <- scan(file.path(candidate_results_dir, "canonical_genes.txt"), what = "", sep = ",")
extended_genes <- scan(file.path(candidate_results_dir, "extended_genes.txt"), what = "", sep = ",")
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 3. INTEGRATE DEA RESULTS WITH CANDIDATE GENES
cat("Integrating DEA results with candidate genes...\n")
# 3.1. Sort DEA resuls by adj.P.Val and generate a rank
deg_all_ordered <- deg_results %>%
  arrange(adj.P.Val) %>%
  mutate(rank_DEA = row_number())
# 3.2. Join DEA results and aggregated results
aggregated_dea_ranked <- aggregated_shap_data %>%
  left_join(deg_all_ordered, by = c("variable_name" = "gene_symbol")) %>%
    mutate(
        passed_DEA = .data$adj.P.Val < p_val_thresh & abs(.data$logFC) > lfc_thresh,
        gene_status = dplyr::case_when(
            .data$variable_name %in% canonical_genes ~ "Canonical",
            .data$variable_name %in% extended_genes ~ "Extended",
            TRUE ~ "Other"
        )
    )

# Save data
cat("Saving integrated data...\n")
save_helper_csv(aggregated_dea_ranked, "candidate_genes_ranked_in_dea.csv", "test", "candidate_genes")
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 4. PREPARE DATA FOR PLOTTING
cat("Preparing data for plotting...\n")
plot_data <- aggregated_dea_ranked %>%
  mutate(
    is_DEG = adj.P.Val < p_val_thresh & abs(logFC) > lfc_thresh,
    status = case_when(
      variable_name %in% canonical_genes ~ "Canonical",
      variable_name %in% extended_genes ~ "Extended",
      is_DEG ~ "DEG",
      TRUE ~ "Other"
    )
  )
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 5. PLOT MA PLOT AND VOLCANO PLOT
cat("Plotting MA and Volcano plots...\n")
p_ma <- plot_dea_candidates(plot_data, type = "MA")
p_volcano <- plot_dea_candidates(plot_data, type = "Volcano")
# Save the plots
ggplot2::ggsave(file.path(output_dir, "DEA_MA_candidates.pdf"), p_ma, width = 12, height = 8)
ggplot2::ggsave(file.path(output_dir, "DEA_Volcano_candidates.pdf"), p_volcano, width = 12, height = 8)
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 6. PLOT SHAP VS OCCURRENCE FOR CANDIDATE GENES
cat("Plotting SHAP vs Occurrence for candidate genes...\n")
plot_data <- plot_data %>%
  filter(variable_name != "FANCC")
p_scatter <- plot_shap_occurrence(data = plot_data, y_threshold = 0.0107, title = "Candidate Genes: Importance vs Robustness")
# Save the plot
ggplot2::ggsave(file.path(output_dir, "maxSHAP_vs_Occurrence_candidates.pdf"), p_scatter, width = 15, height = 12)
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 7. FINISH
cat("All plots saved in:", output_dir, "\n")