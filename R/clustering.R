#' @importFrom stats dist hclust cutree
#' @importFrom proxy dist
#' @importFrom factoextra fviz_eig
NULL

#' Perform Sample Clustering
#'
#' @description
#' This function computes a distance matrix based on the provided method
#' (either "euclidean" or "cosine") and performs hierarchical clustering
#' on the samples. It then cuts the resulting dendrogram into a specified
#' number of clusters (k) and returns both the hierarchical clustering object
#' and the cluster assignments for each sample.
#'
#' @param mat Matrix. Expression or SHAP matrix (genes as rows, samples as columns).
#' @param method Character. "euclidean" or "cosine".
#' @param k Integer. Number of clusters to cut.
#'
#' @return A list containing the hierarchical clustering object and the cluster assignments.
#'
#' @export
calculate_sample_clustering <- function(mat, method = "cosine", k = 5) {

  set.seed(1234)
  print_message("Computing", method, "distance...")
  if (method == "cosine") {
    d <- proxy::dist(t(mat), method = "cosine")
  } else {
    d <- stats::dist(t(mat), method = "euclidean")
  }
  hc <- stats::hclust(d)
  clusters <- stats::cutree(hc, k = k)

  return(list(hclust = hc, clusters = clusters))
}

#' Enrichment Analysis of Clusters vs Phenotypes
#'
#' @description
#' Performs "Cluster vs Rest" statistical tests for each phenotypic variable.
#' For numeric variables, it tests for normality and applies either a t-test or Wilcoxon test.
#' For categorical variables, it applies a chi-squared test. The function returns
#' a dataframe summarizing the results, including adjusted p-values for multiple testing.
#' Clusters with fewer than 3 samples are automatically skipped to ensure statistical validity.
#' The results include the type of test used, p-values, test statistics, and any enriched classes for categorical variables.
#'
#' @param cluster_vec Named integer vector. Output from `calculate_sample_clustering$clusters`.
#' @param phenotype_data Dataframe. Phenotypic variables (samples as rows).
#'
#' @return A dataframe summarizing the enrichment results for each cluster and phenotypic variable.
#'
#' @export
run_cluster_enrichment <- function(cluster_vec, phenotype_data) {
  
  # 1. Extract unique clusters and prepare results dataframe
  unique_clusters <- sort(unique(cluster_vec))
  results_all <- data.frame()
  # 2. Subset phenotype data to samples in cluster_vec
  common_samples <- intersect(names(cluster_vec), rownames(phenotype_data))
  pdata_sub <- phenotype_data[common_samples, ]
  # 3. Samples in each cluster
  cluster_sub <- cluster_vec[common_samples]
  cluster_sizes <- table(cluster_sub)
  # 4. Iterate over clusters and perform tests
  for (cl_id in unique_clusters) {
    # 4.1 Skip clusters with less than 3 samples
    n_samples <- cluster_sizes[as.character(cl_id)]
    if (n_samples < 3) {
      print_message("Skipping Cluster:", cl_id, "(only", n_samples, "samples)")
      next
    }
    # 4.2. Processing Cluster
    print_message("Processing Cluster:", cl_id)
    # Create binary variable "This Cluster vs The Rest"
    current_pheno <- pdata_sub %>%
      dplyr::mutate(target_cluster = ifelse(cluster_sub == cl_id, 1, 0))
    # 4.2.1 Iterate over each phenotypic variable (excluding target_cluster)
    vars <- setdiff(colnames(pdata_sub), "sample_id")
    for (v in vars) {
      if (is.numeric(current_pheno[[v]])) {
        # Test for Numeric Variables
        normality <- tryCatch(stats::shapiro.test(current_pheno[[v]])$p.value > 0.05, error = function(e) FALSE)
        if (!normality) {
          test <- stats::wilcox.test(current_pheno[[v]] ~ current_pheno$target_cluster)
          method <- "Wilcoxon"
        } else {
          test <- stats::t.test(current_pheno[[v]] ~ current_pheno$target_cluster)
          method <- "T-test"
        }
        res <- data.frame(Variable = v, Cluster = cl_id, Test = method, 
                          p_value = test$p.value, Statistic = test$statistic, Enriched_Class = NA)
      } else {
        # Test for Categorical Variables
        tab <- table(current_pheno$target_cluster, current_pheno[[v]])
        if (nrow(tab) < 2 || ncol(tab) < 2) next
        test <- stats::chisq.test(tab, simulate.p.value = (min(tab) < 5))
        # Identify enriched class (where observed > expected in cluster 1)
        direction <- test$observed[2, ] > test$expected[2, ]
        enriched <- paste(colnames(tab)[direction], collapse = ", ")  
        res <- data.frame(Variable = v, Cluster = cl_id, Test = "Chi-squared", 
                          p_value = test$p.value, Statistic = test$statistic, Enriched_Class = enriched)
      }
      # 4.3. Append results
      results_all <- rbind(results_all, res)
    }
  }
  # 5. Adjust p-values for multiple testing
  results_all$p_adj <- stats::p.adjust(results_all$p_value, method = "BH")

  return(results_all)
}

#' Perform PCA
#'
#' @description
#' Computes a Principal Component Analysis (PCA) on the provided data matrix
#' (samples as rows, features as columns), scaling and centering it beforehand.
#' Useful for inspecting global structure and batch effects. Optionally saves
#' a scree plot (eigenvalues/variance explained per component) to disk.
#'
#' @param data Matrix or dataframe. Samples as rows, features (e.g. genes) as columns.
#' @param num_components Integer. Number of principal components to keep.
#' @param save_eig_plot Logical. If `TRUE`, saves a scree plot of the first
#'   `eig_ncp` components via `factoextra::fviz_eig()`. Default `FALSE`.
#' @param eig_ncp Integer. Number of components to show in the scree plot
#'   (passed to `factoextra::fviz_eig(ncp = eig_ncp)`). Default 15.
#' @param output_dir Character. Directory where the scree plot is saved when
#'   `save_eig_plot = TRUE`. Created if it doesn't exist. Required in that case.
#' @param plot_filename Character. File name for the saved scree plot.
#'   Default `"pca_eigenvalues.png"`.
#'
#' @return A list containing `pca_data` (a dataframe with the retained principal
#'   components, samples as row names) and `pca_summary` (the `summary.prcomp`
#'   object, useful for inspecting variance explained per component).
#'
#' @export
perform_pca <- function(data, num_components = 30,
                         save_eig_plot = FALSE, eig_ncp = 15,
                         output_dir = NULL, plot_filename = "pca_eigenvalues.png") {
 
  print_message("Computing PCA...")
  pca <- stats::prcomp(data, scale. = TRUE, center = TRUE)
  max_comp <- min(num_components, ncol(pca$x))
 
  if (save_eig_plot) {
    if (is.null(output_dir)) stop("output_dir must be provided when save_eig_plot = TRUE")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  print_message("Saving PCA eigenvalue plot to", file.path(output_dir, plot_filename))

  eig_plot <- factoextra::fviz_eig(pca, ncp = eig_ncp, addlabels = TRUE) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 18),
      axis.title   = ggplot2::element_text(size = 16),   # nombres de los ejes
      axis.text    = ggplot2::element_text(size = 14)    # números/etiquetas de los ejes
    )
  # Tamaño de los % sobre las barras
  eig_plot$layers[[length(eig_plot$layers)]]$aes_params$size <- 4.5

  base_name <- file.path(output_dir, sub("\\.[^.]+$", "", plot_filename))  # nombre sin extensión

  ggplot2::ggsave(paste0(base_name, ".svg"), plot = eig_plot,
                  width = 6, height = 5, device = "pdf")

  ggplot2::ggsave(paste0(base_name, ".jpeg"), plot = eig_plot,
                  width = 6, height = 5, dpi = 600, quality = 100, bg = "white")
  }
 
  return(list(
    pca_data = data.frame(pca$x[, seq_len(max_comp)]),
    pca_summary = summary(pca)
  ))
}
 

#' Join PCA Coordinates with Sample Metadata
#'
#' @description
#' Joins the PCA coordinates (as returned in `perform_pca()$pca_data`) with
#' sample-level phenotypic/metadata, matching by `sample_id`.
#'
#' @param pca_data Dataframe. PCA coordinates with samples as row names, as
#'   returned in `perform_pca()$pca_data`.
#' @param metadata Dataframe. Sample metadata, must contain a `sample_id` column.
#'
#' @return A dataframe with the PCA coordinates joined to the metadata.
#'
#' @export
join_pca_with_metadata <- function(pca_data, metadata) {

  pca_df <- pca_data %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::inner_join(metadata, by = "sample_id")

  return(pca_df)
}

## ----------------------------------------------------------------------------------------------------------------------------------------

#' Perform UMAP
#'
#' @description
#' Computes a UMAP embedding on the provided data matrix (samples as rows,
#' features as columns), scaling it beforehand.
#'
#' @param data Matrix or dataframe. Samples as rows, features (e.g. genes) as columns.
#'
#' @return A `umap` object, as returned by `umap::umap()`.
#'
#' @export
perform_umap <- function(data) {

  print_message("Computing UMAP...")
  umap_data <- data %>%
    base::scale() %>%
    umap::umap()

  return(umap_data)
}

#' Join UMAP Coordinates with Sample Metadata
#'
#' @description
#' Extracts the UMAP layout (UMAP1/UMAP2) from a `umap` object (as returned by
#' `perform_umap()`) and joins it with sample-level phenotypic/metadata,
#' matching by `sample_id`.
#'
#' @param umap_data List. A `umap` object, as returned by `perform_umap()`.
#' @param metadata Dataframe. Sample metadata, must contain a `sample_id` column.
#'
#' @return A dataframe with the UMAP1/UMAP2 coordinates joined to the metadata.
#'
#' @export
join_umap_with_metadata <- function(umap_data, metadata) {

  umap_df <- umap_data$layout %>%
    as.data.frame()
  colnames(umap_df) <- c("UMAP1", "UMAP2")

  umap_df <- umap_df %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::inner_join(metadata, by = "sample_id")

  return(umap_df)
}