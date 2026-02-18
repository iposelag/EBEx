#' @importFrom dplyr filter select mutate group_by summarize arrange n desc ungroup bind_rows
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map imap_dfr
#' @importFrom magrittr %>%
#' @importFrom rlang .data !! sym
NULL

#' Filter best scenarios (input list + classifier) based on performance metrics
#'
#' @description
#' This function identifies which combinations of input gene lists and classifiers meet
#' a specified performance threshold based on both cross-validation and test metrics.
#'
#' @param results_models_list A named list of performance objects.
#' @param metric Character. Metric to filter by (e.g., "normMCC").
#' @param threshold Numeric. Minimum value for inclusion.
#'
#' @return A data frame of valid (input_list, classifier) pairs.
#'
#' @export
filter_model_performance <- function(results_models_list, metric = "normMCC", threshold = 0.75) {
  
  valid_pairs <- purrr::imap_dfr(results_models_list, function(res, list_name) {
    # Check Cross-Validation
    cv_metrics <- process_cross_validation_metrics(res, names(res$cross_validation))
    valid_cv <- cv_metrics %>% 
      dplyr::filter(.data$metric == !!metric, .data$mean > !!threshold) %>%
      dplyr::select(.data$classifier)
    # Check Test
    test_metrics <- process_test_metrics(res, names(res$test))
    valid_test <- test_metrics %>% 
      dplyr::filter(.data$metric == !!metric, .data$estimate > !!threshold) %>%
      dplyr::select(.data$classifier)
    # Join valid classifiers from both CV and Test
    valid_clfs <- unique(c(valid_cv$classifier, valid_test$classifier))
    if(length(valid_clfs) > 0) {
      return(data.frame(input_list = list_name, classifier = valid_clfs))
    } else {
      return(NULL)
    }
  })

  return(valid_pairs)
}

#' Normalize SHAP Values
#'
#' @description
#' This function takes a dataframe of SHAP values and normalizes them within each sample.
#' The normalization is done by dividing the mean SHAP value of each variable by the sum 
#' of absolute mean SHAP values for that sample, resulting in a 'norm_shap' value that 
#' represents the relative importance of each gene within that sample.
#'
#' @param shap_df A dataframe with columns variable_name, contribution, and sample.
#'
#' @return A dataframe with an added 'norm_shap' column.
#'
#' @export
normalize_shap_values <- function(shap_df) {

  # Normalize SHAP values within each sample
  res <- shap_df %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::mutate(norm_shap = .data$mean_shap / sum(abs(.data$mean_shap))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(.data$norm_shap)) %>% # Optional: sort by normalized_shap
    dplyr::filter(.data$norm_shap != 0)
  
  return(res)
}

#' Generate Gene-Level Aggregation
#'
#' @description
#' This function takes a list of normalized SHAP dataframes (one per model) and
#' aggregates the importance of each gene across all models. It calculates several
#' metrics for each gene, including the maximum, mean, median, standard deviation of
#' the mean normalized SHAP values, and the count of how many models included that gene.
#' The resulting summary dataframe is sorted by the maximum aggregated importance, allowing
#' for easy identification of top candidate genes.
#'
#' @param normalized_shap_list A named list of normalized SHAP dataframes.
#'
#' @return A summary dataframe with metrics per gene.
#'
#' @export
aggregate_genes <- function(normalized_shap_list) {

  # Aggregate importance across models for each gene
  res <- purrr::imap_dfr(normalized_shap_list, function(df, name) {
    df %>%
      dplyr::group_by(.data$variable_name) %>%
      dplyr::summarize(mean_norm_shap = mean(abs(.data$norm_shap)), .groups = "drop") %>%
      dplyr::mutate(input_classifier = name)
  }) %>%
    dplyr::group_by(.data$variable_name) %>%
    dplyr::summarize(
      aggregated_max = max(.data$mean_norm_shap, na.rm = TRUE),
      aggregated_mean = mean(.data$mean_norm_shap, na.rm = TRUE),
      aggregated_median = stats::median(.data$mean_norm_shap, na.rm = TRUE),
      sd = stats::sd(.data$mean_norm_shap, na.rm = TRUE),
      occurrence_count = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$aggregated_max))

  return(res)
}

#' Generate Gene-Sample Matrix
#'
#' @description
#' This function takes a list of normalized SHAP dataframes (one per model) and
#' creates a wide matrix where rows represent genes and columns represent samples.
#' The values in the matrix are the mean normalized SHAP values for each gene-sample
#' combination, averaged across all models. The 'absolute' parameter allows the user 
#' to choose whether to use the absolute values of the normalized SHAP values or to 
#' keep their original sign, depending on whether they want to focus on the magnitude 
#' of importance or also consider the direction of the effect.
#'
#' @param normalized_shap_list A named list of normalized SHAP dataframes.
#' @param absolute Logical. If TRUE, use absolute values of normalized SHAP; if FALSE, use original values.
#'
#' @return A wide matrix of Gene x Sample.
#'
#' @export
aggregate_gene_sample_matrix <- function(normalized_shap_list, absolute) {

  # Combine all normalized SHAP dataframes and calculate mean normalized SHAP for each gene-sample pair
  res <- dplyr::bind_rows(normalized_shap_list, .id = "model") %>%
    dplyr::group_by(.data$variable_name, .data$sample) %>%
    # dplyr::summarise(mean_norm_shap = mean(abs(.data$norm_shap)), .groups = "drop") %>%
    dplyr::summarise(
      mean_norm_shap = if (absolute) {
        mean(abs(.data$norm_shap), na.rm = TRUE)
      } else {
        mean(.data$norm_shap, na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = .data$sample, values_from = .data$mean_norm_shap)

  return(res)
}

#' Full Aggregation Pipeline (Example)
#'
#' @description
#' This demonstrates how to use the individual steps to go from raw results to
#' a candidate gene list.
#'
#' @param results_models_list A named list of performance objects.
#' @param results_shap_list A named list of SHAP dataframes corresponding to the models.
#' @param threshold Numeric. Performance threshold for filtering models.
#' @param absolute Logical. Whether to use absolute values for the gene-sample matrix.
#'
#' @export
run_aggregation_pipeline <- function(results_models_list, results_shap_list, threshold = 0.75, absolute = TRUE, candidate_results_dir = NULL) {
  
  # 1. Filter valid list-classifier combinations
  valid_pairs <- filter_model_performance(results_models_list, threshold = threshold)
  # 2. Extract and Normalize SHAP for those pairs
  normalized_list <- list()
  for(i in 1:nrow(valid_pairs)) {
    inp <- valid_pairs$input_list[i]
    clf <- valid_pairs$classifier[i]
    raw_shap <- results_shap_list[[inp]][[clf]]$shap_values
    normalized_list[[paste(inp, clf, sep="_")]] <- normalize_shap_values(raw_shap)
  }
  # Save normalized SHAP list for reference
  if(!is.null(candidate_results_dir)) {
    print_message("Saving list of normalized SHAP values for valid models...")
    save_helper_rds(normalized_list, "normalized_shap_list.Rds", candidate_results_dir)
  }
  # 3. Generate the summary results
  gene_summary <- aggregate_genes(normalized_list)
  # 4. Generate the sample matrix
  gene_sample_matrix <- aggregate_gene_sample_matrix(normalized_list, absolute)

  return(list(
    summary = gene_summary, 
    matrix = gene_sample_matrix,
    pairs = valid_pairs
  ))
}