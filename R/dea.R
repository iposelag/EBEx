#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable decideTests vennDiagram plotMD
#' @importFrom stats as.formula model.matrix
#' @importFrom grDevices png dev.off
#' @importFrom graphics hist
NULL

#' Run Differential Expression Analysis using Limma
#'
#' @description
#' This function performs differential expression analysis using the Limma package. It constructs a design matrix
#' based on the specified grouping variable, fits a linear model to the expression data, and computes contrasts
#' to identify differentially expressed genes. The results are saved as both RDS and CSV files in the specified output directory.
#' The function returns the fitted model object containing the DEA results.
#'
#' @param expr_matrix Matrix or dataframe. Expression data (genes as rows, samples as columns).
#' @param metadata Dataframe. Sample metadata.
#' @param gene_names Dataframe or vector. Gene symbols or annotations.
#' @param output_dir Character. Directory to save results.
#' @param group_col Character. Name of the column in metadata for grouping.
#' @param contrast Character. Contrast strings (e.g., "COPD - CTRL").
#' @param design_formula Formula. Optional custom design formula.
#' @param colnames_design_formula Character vector. Optional names for design matrix columns.
#'
#' @return MArrayLM object. The fitted model containing DEA results.
#'
#' @export
run_dea_limma <- function(expr_matrix, metadata, gene_names, output_dir, group_col, contrast, 
                          design_formula = NULL, colnames_design_formula = NULL) {
  
  stopifnot(group_col %in% colnames(metadata))
  expr_matrix <- as.matrix(expr_matrix)
  # 1. Build Design Matrix
  if (is.null(design_formula)) {
    formula_obj <- stats::as.formula(paste0("~ 0 + ", group_col))
    mod <- stats::model.matrix(formula_obj, data = metadata)
    colnames(mod) <- levels(factor(metadata[[group_col]]))
  } else {
    mod <- stats::model.matrix(design_formula, data = metadata)
    if (!is.null(colnames_design_formula)) colnames(mod) <- colnames_design_formula
  }
  # 2. Fit Linear Model
  cat("Fitting linear model...\n")
  fit0 <- limma::lmFit(expr_matrix, mod)
  # 3. Contrast Matrix
  contrast.matrix <- limma::makeContrasts(contrasts = contrast, levels = mod)
  # 4. Compute DEGs
  fit <- limma::contrasts.fit(fit0, contrast.matrix)
  fit <- limma::eBayes(fit)
  # Add gene annotations
  fit$genes <- gene_names
  # 5. Save results using package helpers
  if (!is.null(output_dir)) {
    # Save RDS fit
    save_helper_rds(fit, "fit.rds", output_dir, "")
    # Save TopTable CSV
    deg_results <- limma::topTable(fit, coef = 1, number = Inf, adjust.method = "fdr", sort.by = "p")
    save_helper_csv(deg_results, "deg_results.csv", output_dir, "")
  }

  return(fit)
}

#' Post-DEA Analysis and Diagnostics
#'
#' @description
#' This function takes the fitted model object from the DEA analysis and performs post-analysis diagnostics.
#' It generates several plots to visualize the results, including a Venn diagram of DEGs, a histogram of p-values,
#' and a mean-difference (MD) plot. The function also returns a matrix indicating which genes are differentially
#' expressed based on the specified p-value and log-fold change thresholds. The generated plots are saved in a
#' "DEA_plots" subdirectory within the specified output directory.
#'
#' @param fit_object MArrayLM object. Output from `run_dea_limma`.
#' @param pv Numeric. P-value threshold.
#' @param lfc Numeric. Log-fold change threshold.
#' @param coef Integer. Coefficient to analyze.
#' @param output_dir Character. Directory to save plots.
#'
#' @return The results Matrix with DEGs.
#'
#' @export
dea_analysis <- function(fit_object, pv = 0.05, lfc = 0, coef = 1, output_dir) {
  
  # 1. Decide DEGs based on thresholds
  results_decide <- limma::decideTests(fit_object, adjust.method = "fdr", p.value = pv, lfc = lfc)
  # 2. Plots
  # Venn Diagram
  grDevices::png(file.path(output_dir, "dea_venn.png"), width = 800, height = 800)
  limma::vennDiagram(results_decide, include = c("up", "down"), main = "DEGs Venn Diagram")
  grDevices::dev.off()
  # Histogram
  grDevices::png(file.path(output_dir, "dea_pvalue_hist.png"), width = 800, height = 600)
  graphics::hist(fit_object$p.value[, coef], main = "Distribution of p-values", 
                 xlab = "P-value", col = "gray", breaks = 50)
  grDevices::dev.off()
  # MD Plot
  grDevices::png(file.path(output_dir, "dea_MDplot.png"), width = 800, height = 600)
  limma::plotMD(fit_object, coef = coef, status = results_decide[, coef], 
                main = "Mean-Difference Plot", hl.col = c("red", "blue"))
  grDevices::dev.off()
  
  return(results_decide)
}