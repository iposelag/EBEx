#!/usr/bin/env Rscript
######################################
## Enrichment of candidate genes in the phenotypic variables
## Iria Pose
## 19 nov 2025
######################################

## --------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(devtools)
devtools::load_all()
devtools::document()
library(MyPipelinePkg)
library(dplyr)
library(tidyplots)

## --------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
enrichment_results_dir <- "test/analysis/enrichment_candidate_genes"
output_dir <- "test/plots"
raw_data_dir <- "test/raw_data/"
candidate_results_dir <- "test/candidate_genes"

## --------------------------------------------------------------------------------------------------------------------
# 1. LOAD DATA
print_message("Loading data...")
load(file.path(raw_data_dir, "phenotypic.Rda"))
candidate_genes <- read.csv(file.path(candidate_results_dir, "candidate_genes_density.csv"))$variable_name
canonical_genes <- load_gene_list("canonical_genes", file.path(candidate_results_dir))
extended_genes <- load_gene_list("extended_genes", file.path(candidate_results_dir))
shap_matrix <- read.csv(file.path(candidate_results_dir, "aggregated_shap_matrix.csv"), row.names = 1)
print_message("________________________________")

## --------------------------------------------------------------------------------------------------------------------
# 2. PREPARE DATA
print_message("Preparing data for enrichment analysis...")
# Adjust SHAP matrix sample names
samples_shap_matrix <- colnames(shap_matrix)
samples_shap_matrix <- gsub("^(COPD_|CTRL_)", "", samples_shap_matrix)
colnames(shap_matrix) <- samples_shap_matrix
shap_matrix_df <- as.data.frame(t(shap_matrix))
# Filter and clean phenotypes
pheno_clean <- phenotypic_ctics %>%
  dplyr::select(dis_condition, sex, smoker, age, emphysema, pred_dlco, GOLD_stage, pneumocystis_colonization, pred_fev_postbd, pred_fev_prebd, pred_fvc_postbd, pred_fvc_prebd) %>%
  dplyr::filter(dis_condition == "COPD") %>%
  dplyr::mutate(
    across(c(sex, smoker, GOLD_stage), as.factor),
    across(c(age, emphysema, pred_dlco, pneumocystis_colonization, pred_fev_postbd, pred_fev_prebd, pred_fvc_postbd, pred_fvc_prebd), as.numeric),
    smoker_grouped = case_when(
      smoker == "1-Current" ~ "Current",
      TRUE ~ "Other"
    )
  )
# Synchronize samples between SHAP and Phenotype
common_samples <- intersect(rownames(shap_matrix_df), rownames(pheno_clean))
shap_sub <- shap_matrix_df[common_samples, candidate_genes]
pheno_sub <- pheno_clean[common_samples, ]
print_message("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 3. RUN ENRICHMENT
print_message("Running enrichment analysis...")
vars_interest <- c("sex", "smoker", "smoker_grouped","age", 
                    "emphysema", "pred_dlco", "GOLD_stage", "pneumocystis_colonization",
                    "pred_fev_postbd", "pred_fev_prebd", "pred_fvc_postbd", "pred_fvc_prebd")
enrich_results <- run_phenotype_enrichment(
  shap_matrix = shap_sub,
  phenotype_data = pheno_sub,
  genes = candidate_genes,
  variables = vars_interest
)
# Save results
save_helper_csv(enrich_results, "enrichment_full_results.csv", enrichment_results_dir)
# Create enrichment matrix pvalues
enriched_matrix <- get_enrichment_matrix(enrich_results)
# Save enrichment matrix
save_helper_csv(enriched_matrix, "enrichment_matrix_padj.csv", enrichment_results_dir)
# Create enrichment logical matrix
enriched_logical_matrix <- get_enrichment_logical_matrix(enrich_results)
# List of genes enriched per variable
enriched_per_variable <- enriched_logical_matrix %>%
  pivot_longer(-gene, names_to = "variable", values_to = "enriched") %>%
  filter(enriched == TRUE) %>%
  group_by(variable) %>%
  summarise(enriched_genes = list(unique(gene))) %>%
  ungroup()
print_message("Enriched genes per variable:\n")
print_message(enriched_per_variable)
enriched_logical_matrix <- enriched_logical_matrix %>%
  dplyr::mutate(
    canonical_gene = ifelse(gene %in% canonical_genes, TRUE, FALSE),
    extended_gene = ifelse(gene %in% extended_genes, TRUE, FALSE)
  )
# Save enrichment matrix
save_helper_csv(enriched_logical_matrix, "enrichment_logical_matrix.csv", enrichment_results_dir)
print_message("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 4. Visualization (Tidyplots)
print_message("Creating enrichment plots...")
# Create long version for plotting
plot_data <- shap_sub %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  tidyr::pivot_longer(-sample_id, names_to = "Gene", values_to = "mean_shap") %>%
  left_join(tibble::rownames_to_column(pheno_sub, "sample_id"), by = "sample_id")
# Emphysema
enriched_emphysema <- enrich_results %>% 
  filter(variable == "emphysema", p_adj < 0.05) %>% 
  pull(gene)
plot_data %>%
  filter(Gene %in% enriched_emphysema) %>%
  tidyplot(x = emphysema, y = mean_shap) %>%
  add_data_points() %>%
  split_plot(by = Gene) %>%
  save_plot(file.path(output_dir, "enrichment_emphysema.pdf"))
# Predicted DLCO
enriched_dlco <- enrich_results %>% 
  filter(variable == "pred_dlco", p_adj < 0.05) %>% 
  pull(gene)
plot_data %>%
  filter(Gene %in% enriched_dlco) %>%
  tidyplot(x = pred_dlco, y = mean_shap) %>%
  add_data_points() %>%
  split_plot(by = Gene) %>%
  save_plot(file.path(output_dir, "enrichment_pred_dlco.pdf"))
# Predicted FEV1 post bronchodilator
enriched_fev1_postbd <- enrich_results %>% 
  filter(variable == "pred_fev_postbd", p_adj < 0.05) %>% 
  pull(gene)
plot_data %>%
  filter(Gene %in% enriched_fev1_postbd) %>%
  tidyplot(x = pred_fev_postbd, y = mean_shap) %>%
  add_data_points() %>%
  split_plot(by = Gene) %>%
  save_plot(file.path(output_dir, "enrichment_pred_fev1_postbd.pdf"))
# Predicted FEV1 pre bronchodilator
enriched_fev1_prebd <- enrich_results %>% 
  filter(variable == "pred_fev_prebd", p_adj < 0.05) %>% 
  pull(gene)
plot_data %>%
  filter(Gene %in% enriched_fev1_prebd) %>%
  tidyplot(x = pred_fev_prebd, y = mean_shap) %>%
  add_data_points() %>%
  split_plot(by = Gene) %>%
  save_plot(file.path(output_dir, "enrichment_pred_fev1_prebd.pdf"))

## --------------------------------------------------------------------------------------------------------------------
# 5. FINISH
print_message("Enrichment analysis completed. Results saved in:", output_dir)
print_message("Plots saved in:", output_dir)
print_message("________________________________\n")
