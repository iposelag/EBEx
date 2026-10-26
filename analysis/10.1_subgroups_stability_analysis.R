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
library(EBEx)
library(tibble)
library(ComplexHeatmap)
library(mclust)       # adjustedRandIndex
library(ggplot2)
library(ggalluvial)
library(fpc)

## --------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
heatmap_data <- "expression"
raw_data_dir <- "test/raw_data/"
preprocessing_dir <- "test/analysis/preprocessing"
candidate_results_dir <- "test/candidate_genes"
output_dir <- "test/plots"
subgroups_dir <- "test/analysis/subgroups"
set.seed(1234)


stability_dir <- file.path(subgroups_dir, "stability", heatmap_data)
dir.create(stability_dir, recursive = TRUE, showWarnings = FALSE)
k_ref <- 5
k_range <- 2:10
min_cluster_size <- 1       # same threshold as the clinical enrichment tests
alluvial_k <- c(2, 4, 7)    # k-means solutions shown in the alluvial plot
# Code -> manuscript cluster labels (expression clustering). NULL keeps the code labels
paper_labels <- if(heatmap_data == "expression") c("1" = "5", "2" = "3", "3" = "1", "4" = "2", "5" = "4") else NULL
# Cluster colours as in Fig. 3e, named by manuscript label, e.g. c("1" = "#...", "2" = "#..."). NULL = ggplot defaults
cluster_colors <- c("1" = "#ff595e", "2" = "#ffca3a", "3" = "#8ac926", "4" = "#1982c4", "5" = "#6a4c93")

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

## --------------------------------------------------------------------------------------------------------------------
# 4. CLUSTERING STABILITY ACROSS METHODS
print_message("Assessing clustering stability across methods...")
# 4.1 Reference clustering (samples x genes, same scaled & capped matrix as section 3)
print_message("Preparing reference clustering...")
sample_matrix <- t(heatmap_matrix)
cl_ref <- cl_res$clusters[rownames(sample_matrix)]
dist_euclidean <- dist(sample_matrix, method = "euclidean")
hc_ref <- hclust(dist_euclidean, method = "complete")
# Check that the reference reproduces calculate_sample_clustering()
if(adjustedRandIndex(cutree(hc_ref, k = k_ref), cl_ref) < 1) {
  print_message("WARNING: hclust(complete) does not reproduce cl_res$clusters, check calculate_sample_clustering()")
}
# Relabel clusters with the manuscript numbering
if(!is.null(paper_labels)) {
  cl_ref <- setNames(paper_labels[as.character(cl_ref)], names(cl_ref))
}
# Clusters large enough to be tested clinically (singletons are excluded from the plots)
clusters_tested <- names(which(table(cl_ref) >= min_cluster_size))
print(table(cl_ref))

# 4.2 Alternative partitions for each k
print_message("Computing alternative partitions...")
alternative_trees <- list(
  hc_ward_euclidean = hclust(dist_euclidean, method = "ward.D2"),
  hc_complete_manhattan = hclust(dist(sample_matrix, method = "manhattan"), method = "complete"),
  hc_complete_correlation = hclust(as.dist(1 - cor(t(sample_matrix))), method = "complete")
)
alternative_partitions <- list()
set.seed(1234)
alternative_partitions$kmeans_euclidean <- lapply(k_range, function(k) {
  kmeans(sample_matrix, centers = k, nstart = 100, iter.max = 100)$cluster
})
for(method_name in names(alternative_trees)) {
  alternative_partitions[[method_name]] <- lapply(k_range, function(k) cutree(alternative_trees[[method_name]], k = k))
}
alternative_partitions <- lapply(alternative_partitions, setNames, paste0("k", k_range))
method_labels <- c(
  kmeans_euclidean = "k-means\n(Euclidean)",
  hc_ward_euclidean = "Ward\n(Euclidean)",
  hc_complete_manhattan = "Complete\n(Manhattan)",
  hc_complete_correlation = "Complete\n(1 - Pearson)"
)

# 4.3 Top-level split of the dendrogram (k = 2) vs k-means (k = 2), and its association with disease condition
print_message("Comparing top-level split with k-means...")
branch <- cutree(hc_ref, k = 2)
kmeans_k2 <- alternative_partitions$kmeans_euclidean[["k2"]]
branch_table <- table(HC_k2 = branch, KMEANS_k2 = kmeans_k2)
branch_concordance <- max(sum(diag(branch_table)), sum(branch_table[cbind(1:2, 2:1)])) / sum(branch_table)
print(branch_table)
print_message(paste0("Top-level split: ARI = ", round(adjustedRandIndex(branch, kmeans_k2), 2),
                     "; concordant samples = ", round(100 * branch_concordance, 1), "%"))
condition_table <- table(branch = branch, condition = annotation_data$dis_condition)
print(condition_table)
print(chisq.test(condition_table))

# 4.4 Agreement of the full partition (reported in the text only)
print_message("Calculating ARI of the full partition...")
ari_results <- do.call(rbind, lapply(names(alternative_partitions), function(method_name) {
  data.frame(
    method = method_name,
    k = k_range,
    ARI = sapply(alternative_partitions[[method_name]], adjustedRandIndex, y = cl_ref)
  )
}))
ari_summary <- ari_results %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(ARI_min = round(min(ARI), 2), ARI_max = round(max(ARI), 2))
print(ari_summary)
write.csv(ari_results, file.path(stability_dir, "ari_vs_reference.csv"), row.names = FALSE)

# 4.5 Cluster-level reproducibility: best Jaccard of each reference cluster across k
# (Hennig 2007: > 0.75 stable, 0.6-0.75 pattern, < 0.5 dissolved)
print_message("Calculating cluster-level Jaccard...")
jaccard_index <- function(a, b) length(intersect(a, b)) / length(union(a, b))
best_jaccard_match <- function(reference, partition) {
  reference_clusters <- split(names(reference), reference)
  partition_clusters <- split(names(partition), partition)
  do.call(rbind, lapply(names(reference_clusters), function(cluster_id) {
    jaccard_values <- sapply(partition_clusters, jaccard_index, a = reference_clusters[[cluster_id]])
    data.frame(
      ref_cluster = cluster_id,
      n = length(reference_clusters[[cluster_id]]),
      best_jaccard = max(jaccard_values),
      matched_cluster = names(jaccard_values)[which.max(jaccard_values)]
    )
  }))
}
jaccard_results <- do.call(rbind, lapply(names(alternative_partitions), function(method_name) {
  do.call(rbind, lapply(seq_along(k_range), function(i) {
    partition <- setNames(alternative_partitions[[method_name]][[i]], rownames(sample_matrix))
    data.frame(method = method_name, k = k_range[i], best_jaccard_match(cl_ref, partition))
  }))
}))
write.csv(jaccard_results, file.path(stability_dir, "cluster_jaccard_by_k.csv"), row.names = FALSE)
# Best match across all k for each reference cluster and method
jaccard_summary <- jaccard_results %>%
  dplyr::filter(ref_cluster %in% clusters_tested) %>%
  dplyr::group_by(method, ref_cluster, n) %>%
  dplyr::slice_max(best_jaccard, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(ref_cluster, method)
print(jaccard_summary)
write.csv(jaccard_summary, file.path(stability_dir, "cluster_best_jaccard_across_k.csv"), row.names = FALSE)

# Contingency tables vs k-means (internal check, not reported)
for(k in alluvial_k) {
  print(table(HC = cl_ref, KMEANS = alternative_partitions$kmeans_euclidean[[paste0("k", k)]]))
}

# 4.6 Panel A: alluvial plot, reference clusters -> k-means at selected k
print_message("Plotting alluvial diagram...")
alluvial_samples <- names(cl_ref)[cl_ref %in% clusters_tested]
alluvial_data <- data.frame(sample = alluvial_samples, hc_cluster = cl_ref[alluvial_samples], HC = cl_ref[alluvial_samples])
for(k in alluvial_k) {
  alluvial_data[[paste0("KM_k", k)]] <- as.character(alternative_partitions$kmeans_euclidean[[paste0("k", k)]][alluvial_samples])
}
partition_levels <- c("HC", paste0("KM_k", alluvial_k))
partition_names <- c(paste0("Hierarchical\n(k = ", k_ref, ")"), paste0("k-means\n(k = ", alluvial_k, ")"))
alluvial_long <- alluvial_data %>%
  tidyr::pivot_longer(-c(sample, hc_cluster), names_to = "partition", values_to = "cluster") %>%
  dplyr::mutate(
    partition = factor(partition, levels = partition_levels, labels = partition_names),
    hc_cluster = factor(hc_cluster, levels = sort(clusters_tested))
  )
p_alluvial <- ggplot(alluvial_long, aes(x = partition, stratum = cluster, alluvium = sample, fill = hc_cluster)) +
  geom_flow(alpha = 0.6) +
  geom_stratum(fill = "grey95", colour = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  labs(x = NULL, y = "Samples", fill = "Cluster") +
  theme_bw() +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_text(size = 16),   # "Samples"
    axis.text    = element_text(size = 14),   # nombres de las particiones y números
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )
if(!is.null(cluster_colors)) {
  p_alluvial <- p_alluvial + scale_fill_manual(values = cluster_colors, na.translate =FALSE)
}
ggsave(file.path(stability_dir, "panelA_alluvial_hc_vs_kmeans.pdf"), p_alluvial, width = 8, height = 5)
ggsave(file.path(stability_dir, "panelA_alluvial_hc_vs_kmeans.png"), p_alluvial, width = 8, height = 5, dpi = 500)

# 4.7 Panel B: best Jaccard per cluster and method (k at which it is reached)
print_message("Plotting cluster-level Jaccard...")
jaccard_plot_data <- jaccard_summary %>%
  dplyr::mutate(
    method = factor(method_labels[method], levels = method_labels),
    cluster_label = paste0("Cluster ", ref_cluster, " (n = ", n, ")"),
    cell_label = sprintf("%.2f\n(k = %d)", best_jaccard, k)
  )
cluster_order <- unique(jaccard_plot_data$cluster_label[order(as.numeric(jaccard_plot_data$ref_cluster))])
jaccard_plot_data$cluster_label <- factor(jaccard_plot_data$cluster_label, levels = rev(cluster_order))
p_jaccard <- ggplot(jaccard_plot_data, aes(x = method, y = cluster_label, fill = best_jaccard)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = cell_label), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1), name = "Best\nJaccard") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text       = element_text(size = 10),   # métodos y clusters
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  )
ggsave(file.path(stability_dir, "panelB_jaccard_by_method.pdf"), p_jaccard, width = 6, height = 4)
ggsave(file.path(stability_dir, "panelB_jaccard_by_method.png"), p_jaccard, width = 6, height = 4, dpi = 500)
print_message("________________________________")



################3 NEW heatmap (better visualization)
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
  simple_anno_size = unit(4, "mm"),
  gap = unit(0.5, "mm"),                     # pequeña separación entre pistas
  na_col = "grey80",                         # explícito, para que coincida con la leyenda
  annotation_name_gp= gpar(fontsize = 10),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 9),
                                 title_gp = gpar(fontsize = 10, fontface = "bold"))
)
# 5.3 Heatmap
p1 <- Heatmap(
  name = paste0("Scaled ",heatmap_data, " data"),
  heatmap_matrix,
  na_col = "white",
  col = default_colors,
  cluster_columns = cl_res$hclust,  # Pass the precomputed clustering
  column_split = 5,  # Optional: split columns into 5 clusters
  gap = unit(0.5, "mm"),                     # pequeña separación entre pistas
  row_split = 3,
  row_names_gp = gpar(fontsize = 3.8),
  column_dend_height = unit(10, "mm"),   # dendrograma más bajo
  show_column_names = FALSE,
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
