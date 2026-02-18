#!/usr/bin/env Rscript
######################################
## Subgroup analysis based on SHAP values and phenotypic variables
## Iria Pose
## 19 nov 2025
######################################

## --------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(devtools)
devtools::load_all()
devtools::document()
library(MyPipelinePkg)
library(tibble)
library(ComplexHeatmap)

## --------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
heatmap_data <- "expression"
raw_data_dir <- "test/raw_data/"
preprocessing_dir <- "test/analysis/preprocessing"
candidate_results_dir <- "test/candidate_genes"
output_dir <- "test/plots"
subgroups_dir <- "test/analysis/subgroups"
set.seed(1234)

## --------------------------------------------------------------------------------------------------------------------
# 1. LOAD DATA
print_message("Loading data...")
load(file.path(raw_data_dir, "phenotypic.Rda"))
load(file.path(preprocessing_dir, "expression_train.Rda"))
candidate_genes <- read.csv(file.path(candidate_results_dir, "candidate_genes_density.csv"))$variable_name
normalized_shap_list <- readRDS(file.path(candidate_results_dir, "normalized_shap_list.Rds"))
print_message("________________________________")

## --------------------------------------------------------------------------------------------------------------------
# 2. Prepare data
print_message("Preparing data for clustering and enrichment analysis...")
# 2.1 Generate SHAP Gene x Sample matrix
print_message(paste0("Generating ", heatmap_data, " matrix..."))
if(heatmap_data == "shap") {
    aggregated_shap_matrix <- aggregate_gene_sample_matrix(normalized_shap_list, absolute = FALSE)
    # Filter SHAP matrix to keep only candidate genes
    heatmap_matrix <- aggregated_shap_matrix %>%
    filter(variable_name %in% candidate_genes) %>%
    column_to_rownames(var = "variable_name") %>%
    as.matrix()
    # Adjust SHAP matrix sample names
    samples_heatmap_matrix <- colnames(heatmap_matrix)
    samples_heatmap_matrix <- gsub("^(COPD_|CTRL_)", "", samples_heatmap_matrix)
    colnames(heatmap_matrix) <- samples_heatmap_matrix
    #
}else if(heatmap_data == "expression") {
    # Filter expression matrix to keep only candidate genes
    genes_to_keep <- intersect(colnames(expression_train), candidate_genes)
    heatmap_matrix <- expression_train[, genes_to_keep]
    heatmap_matrix <- t(heatmap_matrix)
    samples_heatmap_matrix <- colnames(heatmap_matrix)
}
# Scale heatmap matrix by rows (genes)
heatmap_matrix <- t(scale(t(heatmap_matrix))) # scale() scale by columns, so we transpose before and after
# Cap the values to a range of -3 to 3 for better visualization
heatmap_matrix[heatmap_matrix > 3] <- 3
heatmap_matrix[heatmap_matrix < -3] <- -3
default_colors <- colorRamp2(seq(min(heatmap_matrix), max(heatmap_matrix), length.out = 3), c("blue", "white", "red"))

# 2.2 Prepare phenotypic data for annotation
print_message("Preparing phenotypic data for annotation...")
# Filter and clean phenotypes
annotation_data <- phenotypic_ctics %>%
  dplyr::select(dis_condition, sex, smoker, age, emphysema, pred_dlco, GOLD_stage, pred_fev_postbd, pred_fvc_postbd) %>%
#   dplyr::filter(dis_condition == "COPD") %>%
  dplyr::mutate(
    across(c(sex, smoker, GOLD_stage), as.factor),
    across(c(age, emphysema, pred_dlco, pred_fev_postbd, pred_fvc_postbd), as.numeric),
    smoker_grouped = case_when(
      smoker == "1-Current" ~ "Current",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::mutate(sample_id = rownames(.)) %>%
  subset(sample_id %in% samples_heatmap_matrix)
# Reorder annotation_data to match heatmap_matrix columns
annotation_data <- annotation_data[match(colnames(heatmap_matrix), annotation_data$sample_id), ]
print_message("________________________________")

## --------------------------------------------------------------------------------------------------------------------
# 3. CLUSTERING
print_message("Performing clustering analysis...")
cl_res <- calculate_sample_clustering(heatmap_matrix, method = "euclidean", k = 5)
clusters <- split(names(cl_res$clusters), cl_res$clusters); print(clusters)
print_message("________________________________")

## ---------------------------------------------------------------------------------------------------
# 4. Run enrichment analysis
print_message("Running enrichment analysis for clusters vs phenotypes...")
enrich_cl <- run_cluster_enrichment(cl_res$clusters, annotation_data)
filtered_results <- enrich_cl[enrich_cl$p_value < 0.05, ]
save_helper_csv(filtered_results, paste0("enrichment_results_clusters_vs_phenotypes_", heatmap_data, ".csv"), subgroups_dir)
print_message("________________________________")

## --------------------------------------------------------------------------------------------------------------------
# 5. Plot heatmap with clusters and annotations
print_message("Plotting heatmap with clusters and annotations...")
# 5.1 Define annotation colors
# Categorical variables
annotation_colors <- list(GOLD_stage = GOLD_stage, dis_condition = dis_condition, sex = sex,
                          smoker = smoker)
# Continuous variables
col_age <- colorRamp2(c(min(annotation_data$age, na.rm=TRUE), max(annotation_data$age, na.rm=TRUE)), c("white", "black"))
col_emphysema <- colorRamp2(c(min(annotation_data$emphysema, na.rm=TRUE), max(annotation_data$emphysema, na.rm=TRUE)), c("white", "red"))
col_pred_dlco <- colorRamp2(c(min(annotation_data$pred_dlco, na.rm=TRUE), max(annotation_data$pred_dlco, na.rm=TRUE)), c("white", "blue"))
col_pred_fev_postbd <- colorRamp2(c(min(annotation_data$pred_fev_postbd, na.rm=TRUE), max(annotation_data$pred_fev_postbd, na.rm=TRUE)), c("white", "green"))
col_pred_fvc_postbd <- colorRamp2(c(min(annotation_data$pred_fvc_postbd, na.rm=TRUE), max(annotation_data$pred_fvc_postbd, na.rm=TRUE)), c("white", "purple"))
# 5.2 Create Annotation Heatmap
ha <- HeatmapAnnotation(
  'Disease condition' = annotation_data$dis_condition,
  # 'GOLD stage' = annotation_data$GOLD_stage,
  'Sex' = annotation_data$sex,
  'Smoker' = annotation_data$smoker,
  'Age'  = annotation_data$age,
  'Emphysema' = annotation_data$emphysema,
  'predicted DLCO' = annotation_data$pred_dlco,
  'predicted FEV1' = annotation_data$pred_fev_postbd,
  'predicted FVC' = annotation_data$pred_fvc_postbd,
  col = list(
    # 'GOLD stage' = GOLD_stage,
    'Disease condition' = dis_condition,
    'Sex' = sex,
    'Smoker' = smoker,
    'Age' = col_age,
    'Emphysema' = col_emphysema,
    'predicted DLCO' = col_pred_dlco,
    'predicted FEV1' = col_pred_fev_postbd,
    'predicted FVC' = col_pred_fvc_postbd
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(2, "mm"),
  annotation_name_gp= gpar(fontsize = 8),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 8),
                                 title_gp = gpar(fontsize = 9))
)
# 5.3 Heatmap
p1 <- Heatmap(
  name = paste0("Scaled ",heatmap_data, " data"),
  heatmap_matrix,
  na_col = "white",
  col = default_colors,
  cluster_columns = cl_res$hclust,  # Pass the precomputed clustering
  column_split = 5,  # Optional: split columns into 5 clusters
  row_split = 3,
  row_names_gp = gpar(fontsize = 3.8),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 5), # Ajusta el tamaño aquí (por ejemplo, 5)
  top_annotation = ha,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 9),
    title_gp = gpar(fontsize = 9)
  )
)
# Open a PDF device
pdf(file.path(output_dir, paste0("subgroups_heatmap_euclidean_", heatmap_data, ".pdf")), width = 14, height = 14) # Adjust width and height as needed
# Draw the heatmap
draw(p1)
# Close the device
dev.off()

## --------------------------------------------------------------------------------------------------------------------
# 6. FINISH
print_message("Clustering and enrichment analysis completed. Results saved in:", subgroups_dir)
print_message("Plots saved in:", output_dir)
print_message("________________________________")
