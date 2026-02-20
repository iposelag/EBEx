#' @importFrom mclust mclustBIC Mclust mclustICL
#' @importFrom graphics plot hist abline par
#' @importFrom grDevices pdf dev.off
NULL

#' Run Mclust Analysis for Gene Ranking
#'
#' @description
#' This function fits Gaussian Mixture Models to the aggregated SHAP values to identify 
#' underlying groups of genes.
#'
#' @param data Dataframe. The gene summary output from aggregation.
#' @param column Character. The column name to analyze (default "aggregated_max").
#' @param G Numeric vector. The number of clusters to test (default 1:9).
#'
#' @return A list containing the Mclust model, BIC, ICL, used data, and column name.
#'
#' @export
run_mclust_analysis <- function(data, column = "aggregated_max", G = 1:9) {
  
  # 1. Extract the column of interest
  X <- data[[column]]
  if (is.null(X)) stop("Column not found in data.")
  # 2. Fit Mclust model and compute BIC and ICL
  print_message("Computing BIC and fitting Mclust model...")
  BIC <- mclust::mclustBIC(X, G = G)
  mod1 <- mclust::Mclust(X, x = BIC, G = G)
  # 3. Compute ICL for the same model
  print_message("Computing ICL...")
  ICL <- mclust::mclustICL(X, G = G)
  # 4. Summary output for the console
  print(summary(mod1))
  print_message("Cluster sizes:")
  print(table(mod1$classification))

  return(list(
    model = mod1,
    BIC = BIC,
    ICL = ICL,
    data = data,
    column_used = column
  ))
}

#' Plot Mclust Diagnostics
#'
#' @description
#' Visualizes BIC, ICL and the resulting clusters relative to the SHAP scores.
#'
#' @param mclust_list List. Output from `run_mclust_analysis`.
#' @param sel_groups Numeric vector. Optional: clusters to highlight with a threshold line.
#' @param output_file Character. Optional: path to save as PDF.
#'
#' @return NULL. Generates plots for diagnostics.
#'
#' @export
plot_mclust_diagnostic <- function(mclust_list, sel_groups = NULL, output_file = NULL) {
  
  # Extract data and model
  X <- mclust_list$data[[mclust_list$column_used]]
  mod1 <- mclust_list$model
  # Save original graphical parameters
  # This ensures we don't "break" the user's display
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  # Set up PDF output if specified
  if (!is.null(output_file)) {
    grDevices::pdf(output_file, width = 12, height = 10)
  }
  # Plot
  # Layout: BIC, ICL, and Classification
  graphics::par(mfrow = c(3, 1))
  # 1. Plot BIC
  graphics::plot(mclust_list$BIC, main = "Mclust BIC")
  # 2. Plot ICL
  graphics::plot(mclust_list$ICL, main = "Mclust ICL")
  # 3. Plot Classification vs Scores
  graphics::plot(mod1, what = "classification", main = "Cluster Assignment")
  # Optional: Histogram with threshold line
  if (!is.null(sel_groups)) {
    graphics::par(mfrow = c(1, 1))
    threshold <- min(X[mod1$classification %in% sel_groups])
    graphics::hist(X, breaks = 100, main = "Threshold based on selected clusters")
    graphics::abline(v = threshold, col = "red", lwd = 2, lty = 2)
  }
  # Close PDF device if it was opened
  if (!is.null(output_file)) {
    grDevices::dev.off()
    print_message("Mclust diagnostic saved to:", output_fil)
  }

  return(invisible(NULL))
}

#' Select Candidate Genes based on Mclust Clusters
#'
#' @description
#' This function identifies candidate genes based on the clusters defined by the Mclust model.
#' It determines a threshold based on the minimum score of the selected clusters and filters the genes
#' accordingly. The resulting candidate genes are sorted by their SHAP scores and returned as a dataframe.
#' The function also prints the defined threshold and the number of selected genes to the console.
#'
#' @param mclust_list List. Output from `run_mclust_analysis`.
#' @param sel_groups Numeric vector. The cluster IDs chosen as "high SHAP".
#'
#' @return A sorted dataframe of selected genes.
#'
#' @export
select_genes_mclust <- function(mclust_list, sel_groups) {
  
  # Extract data and model
  data <- mclust_list$data
  X <- data[[mclust_list$column_used]]
  mod1 <- mclust_list$model
  # Find the threshold: the minimum score among the selected clusters
  threshold <- min(X[mod1$classification %in% sel_groups])
  # Filter genes
  candidates <- data[X >= threshold, ]
  candidates <- candidates[order(-candidates[[mclust_list$column_used]]), ]
  print_message(sprintf("Threshold defined at: %.6f\n", threshold))
  print_message(sprintf("Selected %d genes from clusters: %s\n", 
              nrow(candidates), paste(sel_groups, collapse = ", ")))
  
  return(candidates)
}