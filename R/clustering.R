#' @importFrom stats dist hclust cutree
#' @importFrom proxy dist
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