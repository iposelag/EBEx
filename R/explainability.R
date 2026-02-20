#' @importFrom DALEXtra explain_tidymodels
#' @importFrom DALEX predict_parts
#' @importFrom dplyr mutate group_by summarize select all_of
#' @importFrom magrittr %>%
#' @importFrom rlang .data !! sym
NULL

#' Generate SHAP Contributions for a Single Sample
#'
#' @description
#' This function generates SHAP contributions for a single sample using a DALEX
#' explainer object. It constructs a sample label based on the first column of the 
#' data and the row name, computes the SHAP values, and optionally saves the results 
#' as an RDS file in a structured directory.
#'
#' @param explainer DALEX explainer object.
#' @param data_to_shap Dataframe containing the samples.
#' @param target_var Character. Name of the target variable column.
#' @param sample_idx Integer. Index of the sample to explain.
#' @param directory_to_save Character. Optional directory to save the RDS file.
#'
#' @return A list containing the sample label and the SHAP values dataframe.
#'
#' @export
generate_shap_contributions <- function(explainer, data_to_shap, target_var, sample_idx, directory_to_save = NULL) {

  # Remove target variable from the observation
  observation <- data_to_shap[sample_idx, !(names(data_to_shap) %in% target_var)]
  # Construct sample label: First column value + row name
  sample_label <- paste0(data_to_shap[sample_idx, 1], "_", rownames(data_to_shap)[sample_idx])
  # Calculate SHAP
  shap_values <- DALEX::predict_parts(
    explainer = explainer,
    new_observation = observation,
    type = "shap",
    B = 1,
    random_state = 1
  )
  results <- list(sample = sample_label, shap_values = shap_values)
  if(!is.null(directory_to_save)){
    # Save into explainability/ModelName/sample_name.rds
    subfolder <- file.path("explainability", as.character(explainer$label))
    save_helper_rds(results, paste0(sample_label, ".rds"), directory_to_save, subfolder)
  }

  return(results)
}

#' Calculate Variable Importance for a Batch of Samples
#'
#' @description
#' This function calculates SHAP values for a batch of samples using a DALEX explainer.
#' It iterates over the specified number of samples, generates SHAP contributions for each,
#' and aggregates the results into a summary dataframe containing mean SHAP values for each
#' variable and sample. The function also measures the time taken for the entire process.
#'
#' @param explainer DALEX explainer object.
#' @param train_data Dataframe used for training.
#' @param target_var Character. Name of the target variable column.
#' @param n_samples Integer. Number of samples to process.
#' @param directory_to_save Character. Optional directory to save individual results.
#'
#' @return A list containing the summary dataframe of mean SHAP values and the time taken for the process.
#'
#' @export
calculate_variable_importance_batch <- function(explainer, train_data, target_var, n_samples, directory_to_save = NULL) {
  
  ini <- Sys.time()
  all_shaps <- data.frame()
  # Calculate SHAP values for each sample
  for (i in 1:n_samples) {
    results_shap <- generate_shap_contributions(explainer, train_data, target_var, i, directory_to_save)
    p <- results_shap$shap_values %>%
      dplyr::mutate(sample = results_shap$sample)
    all_shaps <- rbind(all_shaps, p)
  }
  # Aggregate mean contributions per variable and sample
  summary_df <- all_shaps %>% 
    dplyr::group_by(.data$variable_name, .data$sample) %>%
    dplyr::summarize(mean_shap = mean(.data$contribution), .groups = "drop")

  return(list(
    shap_values = summary_df,
    time = Sys.time() - ini
  ))
}

#' Calculate SHAP values for all ML models
#'
#' @description
#' This function iterates over a list of machine learning models, calculates SHAP values for each model using the
#' `calculate_variable_importance_batch` function, and compiles the results into a list. It also ensures that
#' the necessary package for explaining tidymodels is available and handles the renaming of model labels for better
#' readability in the output.
#'
#' @param ml_models_to_run Character vector of model names.
#' @param results_models List of fitted models (output from ML_models).
#' @param expression_data Dataframe. The data used for training.
#' @param target_var Character. The target variable column name.
#' @param directory_to_save Character. Optional directory path.
#'
#' @return A list containing SHAP values for each model and the time taken for the entire process.
#'
#' @export
calculate_all_shap <- function(ml_models_to_run, results_models, expression_data, target_var, directory_to_save = NULL) {
  
  # 0. Check that workflows is available for DALEX
  if (!requireNamespace("workflows", quietly = TRUE)) {
    stop("El paquete 'workflows' es necesario para explicar modelos de tidymodels.")
  }
  all_results_shap <- list()
  # 1. Iterate over models and calculate SHAP values
  for(classif in ml_models_to_run) {
    # 1.1 Rename classifiers by using the utility function to get the pretty name
    classif_pretty_name <- rename_classifier(classif)
    print_message("Calculating SHAP for model:", classif_pretty_name)
    # 1.2 Get the fitted model and calculate SHAP values
    fit <- results_models[[classif]]$fit
    target_var_numeric <- as.numeric(expression_data[[target_var]]) - 1
    # 1.3 Create DALEX explainer for tidymodels and calculate SHAP values
    explainer <- DALEXtra::explain_tidymodels(
      fit,
      data = expression_data %>% dplyr::select(-dplyr::all_of(target_var)),
      y = target_var_numeric,
      label = classif_pretty_name # Use the renamed label here
    )
    n_rows <- nrow(expression_data)
    # 1.4 Calculate SHAP values for the batch of samples
    results_shap <- calculate_variable_importance_batch(
      explainer, 
      expression_data, 
      target_var = target_var, 
      n_samples = n_rows, 
      directory_to_save = directory_to_save
    )
    # 1.5 Store results in the list
    all_results_shap[[classif]] <- results_shap
  }

  return(all_results_shap)
}