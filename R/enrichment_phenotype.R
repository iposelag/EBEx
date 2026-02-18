#' @importFrom stats shapiro.test cor.test t.test wilcox.test aov kruskal.test p.adjust
#' @importFrom broom tidy
#' @importFrom dplyr mutate filter group_by summarize select n ungroup left_join rename across
#' @importFrom tidyr unnest pivot_wider pivot_longer
#' @importFrom purrr map
#' @importFrom rlang .data !! sym
NULL

#' Perform Statistical Test for Gene-Phenotype Association
#'
#' @description
#' This function takes a dataframe containing mean SHAP values for a gene and a phenotypic variable,
#' checks for normality, and automatically selects the appropriate statistical test (correlation for numeric
#' variables, t-test or ANOVA for normally distributed groups, and non-parametric tests otherwise).
#'
#' @param data Dataframe containing 'mean_shap' and the target variable.
#' @param variable Character. Name of the phenotypic variable column.
#'
#' @return A tidy dataframe with test results, including p-values, test type, and normality information.
#'
#' @export
perform_pheno_test <- function(data, variable) {

  # 1. Clean data: remove rows with NAs in the variable or in SHAP
  clean_data <- data[!is.na(data[["mean_shap"]]) & !is.na(data[[variable]]), ]
  n_obs <- nrow(clean_data)
  # 2. Safety check: Need at least 3 observations for shapiro and correlation
  if (n_obs < 3) {
    return(data.frame(
      p.value = NA, 
      variable = variable, 
      test_type = "Insufficient data", 
      normality_p_value = NA
    ))
  }
  # 3. Check normality of SHAP values
  shap_normality <- stats::shapiro.test(clean_data[["mean_shap"]])
  is_normal <- shap_normality$p.value >= 0.05
  
  # 4. Determine variable type and run test
  if (is.numeric(clean_data[[variable]])) {
    # 4.1. Correlation
    method <- if (is_normal) "pearson" else "spearman"
    test_result <- stats::cor.test(clean_data[["mean_shap"]], clean_data[[variable]], 
                                   method = method, use = "complete.obs")
    test_type <- paste0(method, " (Normal: ", is_normal, ")")
  } else {
    # 4.2.Comparison of groups
    n_groups <- length(unique(stats::na.omit(clean_data[[variable]])))
    if (n_groups == 2) {
      test_result <- if (is_normal) stats::t.test(mean_shap ~ get(variable), data = clean_data)
                     else stats::wilcox.test(mean_shap ~ get(variable), data = clean_data)
      test_type <- if (is_normal) "t-test" else "Mann-Whitney U"
    } else {
      test_result <- if (is_normal) stats::aov(mean_shap ~ get(variable), data = clean_data)
                     else stats::kruskal.test(mean_shap ~ get(variable), data = clean_data)
      test_type <- if (is_normal) "ANOVA" else "Kruskal-Wallis"
    }
  }
  
  # 5. Clean results
  res <- broom::tidy(test_result) %>%
    dplyr::mutate(
      variable = variable, 
      test_type = test_type,
      normality_p_value = shap_normality$p.value
    )

  return(res)
}

#' Run Full Phenotype Enrichment Analysis
#'
#' @description
#' This function takes a SHAP matrix (samples x genes) and a phenotype dataframe,
#' performs the `perform_pheno_test` for each gene and specified phenotypic variable, 
#' and compiles the results into a tidy dataframe. It also adjusts p-values for multiple
#' testing within each variable and filters out any residual terms or NAs before 
#' returning the final results. The output includes the gene, variable, test type, 
#' p-values, and adjusted p-values for each gene-phenotype association tested.
#'
#' @param shap_matrix Matrix or Dataframe. Samples as rows, genes as columns.
#' @param phenotype_data Dataframe. Phenotypic variables for the same samples.
#' @param genes Character vector. List of genes to test.
#' @param variables Character vector. List of phenotypic variables to test.
#'
#' @return A tidy dataframe summarizing the results of the phenotype enrichment analysis, including p-values and adjusted p-values for each gene-phenotype association.
#'
#' @export
run_phenotype_enrichment <- function(shap_matrix, phenotype_data, genes, variables) {
  
  # 0. Datafreme for results
  all_results <- list()
  # 1. Run tests for each gene and variable
  for (gene in genes) {
    # 1.1. Merge gene SHAP with phenotypes
    merged_data <- phenotype_data %>%
      dplyr::mutate(sample_id = rownames(phenotype_data)) %>%
      dplyr::left_join(
        dplyr::tibble(
          sample_id = rownames(shap_matrix),
          mean_shap = shap_matrix[[gene]]
        ),
        by = "sample_id"
      )
    # 1.2. Run tests for each variable
    gene_results <- purrr::map_dfr(variables, ~ perform_pheno_test(merged_data, .x))
    all_results[[gene]] <- gene_results %>% dplyr::mutate(gene = gene)
  }
  final_df <- dplyr::bind_rows(all_results) %>%
    dplyr::filter(!grepl("Residual", .data$term, ignore.case = TRUE)) %>%
    dplyr::filter(!is.na(.data$p.value))
  # 2. Adjust P-values
  final_df <- final_df %>%
    dplyr::group_by(.data$variable) %>%
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) %>%
    dplyr::ungroup()
  
  return(final_df)
}

#' Create Enrichment Logical Matrix
#'
#' @description
#' This function takes the output from `run_phenotype_enrichment` and creates a logical
#' matrix indicating which genes are significantly associated with each phenotypic variable
#' based on adjusted p-values. The resulting matrix has genes as rows and phenotypic variables
#' as columns, with TRUE indicating significant association and FALSE otherwise.
#'
#' @param enrichment_df Output from `run_phenotype_enrichment`.
#' @param alpha Numeric. Significance threshold.
#'
#' @return A logical matrix indicating significant gene-phenotype associations.
#'
#' @export
get_enrichment_logical_matrix <- function(enrichment_df, alpha = 0.05) {

  matrix <- enrichment_df %>%
    dplyr::select(.data$gene, .data$variable, .data$p_adj) %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$p_adj) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ !is.na(.) & . < alpha))

  return(matrix)
}

#' Create Enrichment Matrix padjusted values
#'
#' @description
#' This function takes the output from `run_phenotype_enrichment` and creates a matrix of 
#' adjusted p-values for each gene-phenotype association. The resulting matrix has genes 
#' as rows and phenotypic variables as columns, with the values representing the adjusted 
#' p-values for each association.
#'
#' @param enrichment_df Output from `run_phenotype_enrichment`.
#' @param alpha Numeric. Significance threshold.
#'
#' @return A matrix of adjusted p-values for gene-phenotype associations.
#'
#' @export
get_enrichment_matrix <- function(enrichment_df, alpha = 0.05) {

  matrix <- enrichment_df %>%
    dplyr::select(.data$gene, .data$variable, .data$p_adj) %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$p_adj)

  return(matrix)
}