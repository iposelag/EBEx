#' @importFrom recipes recipe step_BoxCox step_corr step_normalize all_numeric_predictors all_outcomes all_predictors
#' @importFrom themis step_downsample
#' @importFrom parsnip rand_forest svm_rbf svm_poly logistic_reg multinom_reg nearest_neighbor boost_tree set_mode set_engine fit
#' @importFrom workflows workflow add_recipe add_model update_model
#' @importFrom tune tune tune_bayes select_best extract_parameter_set_dials fit_resamples control_bayes finalize_model conf_mat_resampled collect_metrics collect_predictions
#' @importFrom yardstick metric_set mcc bal_accuracy accuracy sens spec precision roc_auc mn_log_loss conf_mat
#' @importFrom tidyr unnest
#' @importFrom dials finalize
#' @importFrom bundle bundle
#' @importFrom dplyr select all_of filter bind_cols bind_rows rename group_by summarise mutate
#' @importFrom stats as.formula predict
#' @importFrom rlang !! sym .data
#' @importFrom magrittr %>%
NULL

## --------------------------------------------------------------------------------------------------------------------
# PREPROCESSING RECIPE

#' Create Preprocessing Recipe
#'
#' @description
#' This function creates a preprocessing recipe for machine learning models using the `recipes` package.
#' It includes steps for optional Box-Cox transformation, correlation filtering to remove highly correlated predictors,
#' downsampling to address class imbalance, and optional normalization of predictors. 
#'
#' @param data Data frame containing the training data.
#' @param target_var Character string of the target variable column name.
#' @param transformation Character string specifying the type of transformation to apply (e.g., "boxcox"). Default is NULL (no transformation). 
#' @param normalize Logical. Whether to apply normalization to the predictors. Default is FALSE.
#'
#' @return A `recipes` recipe object that can be used in a machine learning workflow.
#'
#' @keywords internal
create_recipe <- function(data, target_var, transformation = NULL, normalize = FALSE) {

  # 1. Start with a basic recipe
  recipe_obj <- recipes::recipe(stats::as.formula(paste(target_var, "~ .")), data = data)
  # 2. Add optional Box-Cox transformation
  if (!is.null(transformation) && transformation == "boxcox") {
    recipe_obj <- recipe_obj %>% recipes::step_BoxCox(recipes::all_numeric_predictors())
  }
  # 3. Add correlation filtering and downsampling
  recipe_obj <- recipe_obj %>% 
    recipes::step_corr(recipes::all_numeric_predictors(), threshold = 0.85) %>%
    themis::step_downsample(recipes::all_outcomes())
  # 4. Add optional normalization
  if (normalize) {
    recipe_obj <- recipe_obj %>% recipes::step_normalize(recipes::all_predictors())
  }

  return(recipe_obj)
}

## --------------------------------------------------------------------------------------------------------------------
# MODEL FUNCTIONS

#' Random Forest Model Training
#'
#' @description
#' This function trains a Random Forest model using the `parsnip` and `workflows` packages. 
#' It performs hyperparameter tuning using Bayesian optimization with the `tune` package, 
#' evaluates the model using resampling, and saves the fitted model and resampling results as 
#' bundles in a specified directory.
#'
#' @param data Data frame containing the training data.
#' @param target_var Character string of the target variable column name.
#' @param folds Resampling object (e.g., vfold_cv) for cross-validation.
#' @param predictors Character vector of predictor variable names.
#' @param metrics Metrics object for model evaluation.
#' @param control Control object for tuning and resampling.
#' @param directory_to_save Character string specifying the directory to save the model bundles. Default is NULL (do not save).
#'
#' @return A list containing the time taken for training, resampling results, and the fitted model object.
#'
#' @export
rf_model <- function(data, target_var, folds, predictors, metrics, control, directory_to_save = NULL) {

  set.seed(1234)
  ini_rf <- Sys.time()
  # Specify the model with tunable parameters
  spec <- parsnip::rand_forest(mtry = tune::tune(), trees = 1000, min_n = tune::tune()) %>%
    parsnip::set_mode("classification") %>%
    parsnip::set_engine("ranger", oob.error = TRUE, importance = "impurity", seed = 1234)
  # Create the workflow with the recipe and model specification
  wf <- workflows::workflow() %>%
    workflows::add_recipe(create_recipe(data, target_var)) %>%
    workflows::add_model(spec)
  # Extract and finalize parameters for tuning
  params <- spec %>% tune::extract_parameter_set_dials() %>% dials::finalize(predictors)
  # Perform Bayesian optimization tuning
  res <- tune::tune_bayes(wf, resamples = folds, param_info = params, initial = 7, iter = 25, metrics = metrics, control = control)
  best <- tune::select_best(res, metric = "mcc")
  # Finalize the model with the best parameters and fit it to the entire training data
  final_spec <- tune::finalize_model(spec, best)
  final_fit <- wf %>% workflows::update_model(final_spec) %>% parsnip::fit(data = data)
  # Evaluate the final model using resampling
  rs <- final_fit %>% tune::fit_resamples(resamples = folds, metrics = metrics, control = control)
  # Save the fitted model and resampling results as bundles if a directory is specified
  if(!is.null(directory_to_save)){
    rf_bundle_fit <- bundle::bundle(final_fit)
    rf_bundle_rs <- bundle::bundle(rs)
    # Save BOTH bundles in ONE file: rf_model.Rda
    save_helper_multi_rda(
      obj_list = list(rf_bundle_fit, rf_bundle_rs),
      name_list = list("rf_bundle_fit", "rf_bundle_rs"),
      filename = "rf_model.Rda",
      directory = directory_to_save,
      subfolder = "ML_models"
    )
  }

  return(list(time = Sys.time() - ini_rf, resamples = rs, fit = final_fit))
}

#' SVM Radial Model
#'
#' @description
#' This function trains a Support Vector Machine with Radial Basis Function kernel using the `parsnip` 
#' and `workflows` packages. It performs hyperparameter tuning using Bayesian optimization with the 
#' `tune` package, evaluates the model using resampling, and saves the fitted model and resampling 
#' results as bundles in a specified directory.
#'
#' @param data Data frame containing the training data.
#' @param target_var Character string of the target variable column name.
#' @param folds Resampling object (e.g., vfold_cv) for cross-validation.
#' @param metrics Metrics object for model evaluation.
#' @param control Control object for tuning and resampling.
#' @param directory_to_save Character string specifying the directory to save the model bundles. Default is NULL (do not save).
#'
#' @return A list containing the time taken for training, resampling results, and the fitted model object.
#'
#' @export
svm_r_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {

  set.seed(1234)
  ini <- Sys.time()
  # Specify the SVM model with tunable parameters
  spec <- parsnip::svm_rbf(cost = tune::tune(), rbf_sigma = tune::tune()) %>%
    parsnip::set_engine("kernlab") %>%
    parsnip::set_mode("classification")
  # Create the workflow with the recipe and model specification
  wf <- workflows::workflow() %>%
    workflows::add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
    workflows::add_model(spec)
  # Perform Bayesian optimization tuning
  res <- tune::tune_bayes(wf, resamples = folds, param_info = tune::extract_parameter_set_dials(spec), initial = 7, iter = 25, metrics = metrics, control = control)
  # Select the best parameters and finalize the model
  best <- tune::select_best(res, metric = "mcc")
  final_spec <- tune::finalize_model(spec, best)
  final_fit <- wf %>% workflows::update_model(final_spec) %>% parsnip::fit(data = data)
  rs <- final_fit %>% tune::fit_resamples(resamples = folds, metrics = metrics, control = control)
  # Save the fitted model and resampling results as bundles if a directory is specified
  if(!is.null(directory_to_save)) {
    # Create the two bundles
    svm_r_bundle_fit <- bundle::bundle(final_fit)
    svm_r_bundle_rs <- bundle::bundle(rs)
    # Save both bundles in one file: svm_r_model.Rda
    save_helper_multi_rda(
      obj_list = list(svm_r_bundle_fit, svm_r_bundle_rs),
      name_list = list("svm_r_bundle_fit", "svm_r_bundle_rs"),
      filename = "svm_r_model.Rda",
      directory = directory_to_save,
      subfolder = "ML_models"
    )
  }

  return(list(time = Sys.time() - ini, resamples = rs, fit = final_fit))
}

#' SVM Polynomial Model
#'
#' @description
#' This function trains a Support Vector Machine with Polynomial kernel using the `parsnip` and `workflows` packages.
#' It performs hyperparameter tuning using Bayesian optimization with the `tune` package, evaluates the model using resampling,
#' and saves the fitted model and resampling results as bundles in a specified directory.
#'
#' @param data Data frame containing the training data.
#' @param target_var Character string of the target variable column name.
#' @param folds Resampling object (e.g., vfold_cv) for cross-validation.
#' @param metrics Metrics object for model evaluation.
#' @param control Control object for tuning and resampling.
#' @param directory_to_save Character string specifying the directory to save the model bundles. Default is NULL (do not save).
#'
#' @return A list containing the time taken for training, resampling results, and the fitted model object.
#'
#' @export
svm_p_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {

  set.seed(1234)
  ini <- Sys.time()
  # Specify the SVM model with tunable parameters
  spec <- parsnip::svm_poly(cost = tune::tune(), degree = tune::tune()) %>%
    parsnip::set_engine("kernlab") %>%
    parsnip::set_mode("classification")
  # Create the workflow with the recipe and model specification
  wf <- workflows::workflow() %>%
    workflows::add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
    workflows::add_model(spec)
  # Perform Bayesian optimization tuning
  res <- tune::tune_bayes(wf, resamples = folds, param_info = tune::extract_parameter_set_dials(spec), initial = 7, iter = 25, metrics = metrics, control = control)
  # Select the best parameters and finalize the model
  best <- tune::select_best(res, metric = "mcc")
  final_spec <- tune::finalize_model(spec, best)
  final_fit <- wf %>% workflows::update_model(final_spec) %>% parsnip::fit(data = data)
  rs <- final_fit %>% tune::fit_resamples(resamples = folds, metrics = metrics, control = control)
  # Save the fitted model and resampling results as bundles if a directory is specified
  if(!is.null(directory_to_save)) {
    # Create the two bundles
    svm_p_bundle_fit <- bundle::bundle(final_fit)
    svm_p_bundle_rs <- bundle::bundle(rs)
    # Save both bundles in one file: svm_p_model.Rda
    save_helper_multi_rda(
      obj_list = list(svm_p_bundle_fit, svm_p_bundle_rs),
      name_list = list("svm_p_bundle_fit", "svm_p_bundle_rs"),
      filename = "svm_p_model.Rda",
      directory = directory_to_save,
      subfolder = "ML_models"
    )
  }

  return(list(time = Sys.time() - ini, resamples = rs, fit = final_fit))
}

#' GLM Binary Model
#'
#' @description
#' This function trains a Generalized Linear Model (GLM) for binary classification using the 
#' `parsnip` and `workflows` packages. It performs hyperparameter tuning using Bayesian optimization
#' and evaluates the model using resampling. The fitted model and resampling results can be saved as bundles in a specified directory.
#'
#' @param data Data frame containing the training data.
#' @param target_var Character string of the target variable column name.
#' @param folds Resampling object (e.g., vfold_cv) for cross-validation
#' @param metrics Metrics object for model evaluation.
#' @param control Control object for tuning and resampling.
#' @param directory_to_save Character string specifying the directory to save the model bundles. Default is NULL (do not save).
#'
#' @return A list containing the time taken for training, resampling results, and the fitted model object.
#'
#' @export
glm_binary_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {

  set.seed(1234)
  ini <- Sys.time()
  # Specify the GLM model with tunable parameters
  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = tune::tune()) %>%
    parsnip::set_engine("glmnet") %>%
    parsnip::set_mode("classification")
  # Create the workflow with the recipe and model specification
  wf <- workflows::workflow() %>%
    workflows::add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
    workflows::add_model(spec)
  # Perform Bayesian optimization tuning
  res <- tune::tune_bayes(wf, resamples = folds, param_info = tune::extract_parameter_set_dials(spec), initial = 7, iter = 25, metrics = metrics, control = control)
  # Select the best parameters and finalize the model
  best <- tune::select_best(res, metric = "mcc")
  final_spec <- tune::finalize_model(spec, best)
  final_fit <- wf %>% workflows::update_model(final_spec) %>% parsnip::fit(data = data)
  rs <- final_fit %>% tune::fit_resamples(resamples = folds, metrics = metrics, control = control)
  # Save the fitted model and resampling results as bundles if a directory is specified
  if(!is.null(directory_to_save)) {
    # Create the two bundles
    glm_bundle_fit <- bundle::bundle(final_fit)
    glm_bundle_rs <- bundle::bundle(rs)
    # Save both bundles in one file: glm_model.Rda
    save_helper_multi_rda(
      obj_list = list(glm_bundle_fit, glm_bundle_rs),
      name_list = list("glm_bundle_fit", "glm_bundle_rs"),
      filename = "glm_model.Rda",
      directory = directory_to_save,
      subfolder = "ML_models"
    )
  }

  return(list(time = Sys.time() - ini, resamples = rs, fit = final_fit))
}

#' GLM Multiclass Model
#'
#' @description
#' This function trains a Generalized Linear Model (GLM) for multiclass classification using the 
#' `parsnip` and `workflows` packages. It performs hyperparameter tuning using Bayesian optimization 
#' and evaluates the model using resampling. The fitted model and resampling results can be saved 
#' as bundles in a specified directory.
#' 
#' @param data Data frame containing the training data.
#' @param target_var Character string of the target variable column name.
#' @param folds Resampling object (e.g., vfold_cv) for cross-validation
#' @param metrics Metrics object for model evaluation.
#' @param control Control object for tuning and resampling.
#' @param directory_to_save Character string specifying the directory to save the model bundles. Default is NULL (do not save).
#'
#' @return A list containing the time taken for training, resampling results, and the fitted model object.
#'
#' @export
glm_multiclass_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {

  set.seed(1234)
  ini <- Sys.time()
  # Specify the GLM model with tunable parameters
  spec <- parsnip::multinom_reg(penalty = tune::tune(), mixture = tune::tune()) %>%
    parsnip::set_engine("glmnet") %>%
    parsnip::set_mode("classification")
  # Create the workflow with the recipe and model specification
  wf <- workflows::workflow() %>%
    workflows::add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
    workflows::add_model(spec)
  # Perform Bayesian optimization tuning
  res <- tune::tune_bayes(wf, resamples = folds, param_info = tune::extract_parameter_set_dials(spec), initial = 7, iter = 25, metrics = metrics, control = control)
  # Select the best parameters and finalize the model
  best <- tune::select_best(res, metric = "mcc")
  final_spec <- tune::finalize_model(spec, best)
  final_fit <- wf %>% workflows::update_model(final_spec) %>% parsnip::fit(data = data)
  rs <- final_fit %>% tune::fit_resamples(resamples = folds, metrics = metrics, control = control)
  # Save the fitted model and resampling results as bundles if a directory is specified
  if(!is.null(directory_to_save)) {
    # Create the two bundles
    glm_bundle_fit <- bundle::bundle(final_fit)
    glm_bundle_rs <- bundle::bundle(rs)
    # Save both bundles in one file: glm_model.Rda
    save_helper_multi_rda(
      obj_list = list(glm_bundle_fit, glm_bundle_rs),
      name_list = list("glm_bundle_fit", "glm_bundle_rs"),
      filename = "glm_model.Rda",
      directory = directory_to_save,
      subfolder = "ML_models"
    )
  }

  return(list(time = Sys.time() - ini, resamples = rs, fit = final_fit))
}

#' KNN Model
#' 
#' @description
#' This function trains a K-Nearest Neighbors (KNN) model using the `parsnip` and `workflows` packages. 
#' It performs hyperparameter tuning using Bayesian optimization with the `tune` package, evaluates the 
#' model using resampling, and saves the fitted model and resampling results as bundles in a specified directory.
#' 
#' @param data Data frame containing the training data.
#' @param target_var Character string of the target variable column name.
#' @param folds Resampling object (e.g., vfold_cv) for cross-validation
#' @param metrics Metrics object for model evaluation.
#' @param control Control object for tuning and resampling.
#' @param directory_to_save Character string specifying the directory to save the model bundles. Default is NULL (do not save).
#' 
#' @return A list containing the time taken for training, resampling results, and the fitted model object.
#' 
#' @export
knn_model <- function(data, target_var, folds, metrics, control, directory_to_save = NULL) {

  set.seed(1234)
  ini <- Sys.time()
  # Specify the KNN model with tunable parameters
  spec <- parsnip::nearest_neighbor(neighbors = tune::tune(), dist_power = tune::tune(), weight_func = tune::tune()) %>%
    parsnip::set_engine("kknn") %>%
    parsnip::set_mode("classification")
  # Create the workflow with the recipe and model specification
  wf <- workflows::workflow() %>%
    workflows::add_recipe(create_recipe(data, target_var, transformation = "boxcox", normalize = TRUE)) %>%
    workflows::add_model(spec)
  # Perform Bayesian optimization tuning
  res <- tune::tune_bayes(wf, resamples = folds, param_info = tune::extract_parameter_set_dials(spec), initial = 7, iter = 25, metrics = metrics, control = control)
  # Select the best parameters and finalize the model
  best <- tune::select_best(res, metric = "mcc")
  final_spec <- tune::finalize_model(spec, best)
  final_fit <- wf %>% workflows::update_model(final_spec) %>% parsnip::fit(data = data)
  rs <- final_fit %>% tune::fit_resamples(resamples = folds, metrics = metrics, control = control)
  # Save the fitted model and resampling results as bundles if a directory is specified
  if(!is.null(directory_to_save)) {
    # Create the two bundles
    knn_bundle_fit <- bundle::bundle(final_fit)
    knn_bundle_rs <- bundle::bundle(rs)
    # Save both bundles in one file: knn_model.Rda
    save_helper_multi_rda(
      obj_list = list(knn_bundle_fit, knn_bundle_rs),
      name_list = list("knn_bundle_fit", "knn_bundle_rs"),
      filename = "knn_model.Rda",
      directory = directory_to_save,
      subfolder = "ML_models"
    )
  }

  return(list(time = Sys.time() - ini, resamples = rs, fit = final_fit))
}

#' XGBoost Model
#'
#' @description
#' This function trains an Extreme Gradient Boosting (XGBoost) model using the `parsnip` and `workflows` packages. 
#' It performs hyperparameter tuning using Bayesian optimization with the `tune` package, evaluates the model 
#' using resampling, and saves the fitted model and resampling results as bundles in a specified directory.
#'
#' @param data Data frame containing the training data.
#' @param target_var Character string of the target variable column name.
#' @param folds Resampling object (e.g., vfold_cv) for cross-validation
#' @param predictors Character vector of predictor variable names.
#' @param metrics Metrics object for model evaluation.
#' @param control Control object for tuning and resampling.
#' @param directory_to_save Character string specifying the directory to save the model bundles. Default is NULL (do not save).
#'
#' @return A list containing the time taken for training, resampling results, and the fitted model object.
#'
#' @export
xgb_model <- function(data, target_var, folds, predictors, metrics, control, directory_to_save = NULL) {

  set.seed(1234)
  ini <- Sys.time()
  # Specify the XGBoost model with tunable parameters
  spec <- parsnip::boost_tree(trees = 1000, tree_depth = tune::tune(), min_n = tune::tune(), loss_reduction = tune::tune(), sample_size = tune::tune(), mtry = tune::tune(), learn_rate = tune::tune()) %>%
    parsnip::set_engine("xgboost") %>%
    parsnip::set_mode("classification")
  # Create the workflow with the recipe and model specification
  wf <- workflows::workflow() %>%
    workflows::add_recipe(create_recipe(data, target_var)) %>%
    workflows::add_model(spec)
  # Extract and finalize parameters for tuning
  params <- spec %>% tune::extract_parameter_set_dials() %>% dials::finalize(predictors)
  res <- tune::tune_bayes(wf, resamples = folds, param_info = params, initial = 7, iter = 25, metrics = metrics, control = control)
  # Select the best parameters and finalize the model
  best <- tune::select_best(res, metric = "mcc")
  final_spec <- tune::finalize_model(spec, best)
  final_fit <- wf %>% workflows::update_model(final_spec) %>% parsnip::fit(data = data)
  rs <- final_fit %>% tune::fit_resamples(resamples = folds, metrics = metrics, control = control)
  # Save the fitted model and resampling results as bundles if a directory is specified
  if(!is.null(directory_to_save)) {
    # Create the two bundles
    xgb_bundle_fit <- bundle::bundle(final_fit)
    xgb_bundle_rs <- bundle::bundle(rs)
    # Save both bundles in one file: xgb_model.Rda
    save_helper_multi_rda(
      obj_list = list(xgb_bundle_fit, xgb_bundle_rs),
      name_list = list("xgb_bundle_fit", "xgb_bundle_rs"),
      filename = "xgb_model.Rda",
      directory = directory_to_save,
      subfolder = "ML_models"
    )
  }

  return(list(time = Sys.time() - ini, resamples = rs, fit = final_fit))
}

## --------------------------------------------------------------------------------------------------------------------
# MAIN ML WRAPPER

#' Run ML Pipeline
#'
#' @description
#' This function serves as a wrapper to run multiple machine learning models on a given dataset.
#' It takes the training data, target variable, and a list of models to run, and executes the 
#' corresponding model training functions. The results from each model, including training time, 
#' resampling results, and fitted model objects, are collected in a list and returned. Optionally,
#' the trained models and resampling results can be saved to a specified directory.
#'
#' @param data Data frame containing the training data.
#' @param target_var Name of the target variable.
#' @param models_to_run Vector of model names to run.
#' @param directory_to_save Directory to save the trained models and results (optional).
#' 
#' @return A list containing the results from each model, including training time, resampling results, and fitted model objects.
#' 
#' @export
ML_models <- function(data, target_var, models_to_run = c("rf", "svm_r", "svm_p", "glm", "knn", "xgb"), directory_to_save = NULL, verbose = FALSE) {

  set.seed(1234)
  # 1. Define the metrics to save for model evaluation
  metrics_to_save <- yardstick::metric_set(
    yardstick::mcc, yardstick::bal_accuracy, yardstick::accuracy, 
    yardstick::sens, yardstick::spec, yardstick::precision, 
    yardstick::roc_auc, yardstick::mn_log_loss
  )
  # 2. Prepare the predictors and resampling folds
  expression_predictors <- dplyr::select(data, -dplyr::all_of(target_var))
  expression_folds <- rsample::vfold_cv(data, strata = !!rlang::sym(target_var), repeats = 10)
  # Save the folds object if a directory is specified
  if(!is.null(directory_to_save)){
    save_helper_rda(expression_folds, "expression_folds", "folds_object.Rda", directory_to_save, "ML_models")
  }
  # 3. Set control parameters for Bayesian tuning
  ctrl_bayes <- tune::control_bayes(no_improve = 15, verbose = verbose, save_pred = TRUE, save_workflow = TRUE, seed = 1234)
  results <- list()
  # 4. Run the specified models and collect results
  if ("rf" %in% models_to_run) {
    results$rf <- rf_model(data, target_var, expression_folds, expression_predictors, metrics_to_save, ctrl_bayes, directory_to_save)
    cat(rename_classifier("rf"), "done\n")
  }
  if ("svm_r" %in% models_to_run) {
    results$svm_r <- svm_r_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
    cat(rename_classifier("svm_r"), "done\n")
  }
  if ("svm_p" %in% models_to_run) {
    results$svm_p <- svm_p_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
    cat(rename_classifier("svm_p"), "done\n")
  }
  if ("glm" %in% models_to_run) {
    if(length(levels(data[[target_var]])) == 2) {
      results$glm <- glm_binary_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
    } else {
      results$glm <- glm_multiclass_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
    }
    cat(rename_classifier("glm"), "done\n")
  }
  if ("knn" %in% models_to_run) {
    results$knn <- knn_model(data, target_var, expression_folds, metrics_to_save, ctrl_bayes, directory_to_save)
    cat(rename_classifier("knn"), "done\n")
  }
  if ("xgb" %in% models_to_run) {
    results$xgb <- xgb_model(data, target_var, expression_folds, expression_predictors, metrics_to_save, ctrl_bayes, directory_to_save)
    cat(rename_classifier("xgb"), "done\n")
  }

  return(results)
}

## --------------------------------------------------------------------------------------------------------------------
# POSTPROCESSING FUNCTIONS

#' @keywords internal
add_norm_mcc_metric <- function(metrics_df, val_col) {
  # Si no hay MCC, devolvemos las métricas tal cual
  if (!"mcc" %in% metrics_df$.metric) return(metrics_df)
  
  # Extraemos la fila de MCC
  mcc_row <- metrics_df[metrics_df$.metric == "mcc", ]
  
  # Creamos la fila normMCC
  norm_mcc_row <- mcc_row
  norm_mcc_row$.metric <- "normMCC"
  norm_mcc_row[[val_col]] <- (mcc_row[[val_col]] + 1) / 2
  
  # Si es CV, el error estándar ya no aplica al valor transformado
  if ("std_err" %in% colnames(norm_mcc_row)) {
    norm_mcc_row$std_err <- NA
  }

  return(dplyr::bind_rows(metrics_df, norm_mcc_row))
}

#' @keywords internal
extract_train_test_results <- function(fit, data, target_var) {
  # Generar predicciones
  pred <- fit %>% stats::predict(data) %>% dplyr::bind_cols(data)
  
  # Matriz de confusión
  matrix <- yardstick::conf_mat(pred, truth = !!rlang::sym(target_var), estimate = .data$.pred_class)
  
  # Métricas (summary devuelve columnas .metric y .estimate)
  metrics <- summary(matrix)
  metrics <- add_norm_mcc_metric(metrics, val_col = ".estimate")
  
  return(list(
    model_metrics = metrics,
    confusion_matrix = matrix,
    misclassified_samples = pred %>% 
      dplyr::filter(!!rlang::sym(target_var) != .data$.pred_class) %>% 
      rownames()
  ))
}

#' @keywords internal
extract_cross_validation_results <- function(resamples, data, target_var) {
  # Predicciones de CV
  pred <- resamples %>% tune::collect_predictions()
  
  # Muestras mal clasificadas
  misclassified <- pred %>% 
    dplyr::filter(!!rlang::sym(target_var) != .data$.pred_class) %>%
    dplyr::group_by(.data$.row) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  # Métricas de CV (collect_metrics devuelve columnas .metric y mean)
  metrics <- tune::collect_metrics(resamples)
  metrics <- add_norm_mcc_metric(metrics, val_col = "mean")
  
  matrix <- tune::conf_mat_resampled(resamples)
  
  return(list(
    model_metrics = metrics,
    confusion_matrix = matrix,
    misclassified_samples = misclassified
  ))
}

#' Extract models results (cross-validation, train, test and timing)
#'
#' @description
#' This function extracts the results from fitted machine learning models, including 
#' cross-validation metrics, training and test metrics, and timing information. It 
#' takes a list of fitted models and the corresponding training and test datasets,
#' and returns a structured list containing the extracted results for each model. 
#' The function relies on internal helper functions to extract the relevant metrics 
#' and confusion matrices for both training and test datasets, as well as the cross-validation results.
#' 
#' @param models_to_run Vector of models names.
#' @param results_models List of fitted models.
#' @param expression_train Training data.
#' @param expression_test Test data.
#' @param target_var Name of the target variable.
#' 
#' @return A list containing the extracted results for each model, including cross-validation metrics, training and test metrics, confusion matrices, misclassified samples, and timing information.
#' 
#' @export
extract_models_results <- function(models_to_run, results_models, expression_train, expression_test, target_var) {

  # 1. Initialization
  times <- list()
  train <- list()
  test <- list()
  cross_validation <- list()
  extracted_results <- list()
  # 2. Loop through models
  for (classif_name in models_to_run) {
    # 2.1. Extract timing, resampling, and fit
    times[[classif_name]] <- results_models[[classif_name]]$time
    resample <- results_models[[classif_name]]$resamples
    fit <- results_models[[classif_name]]$fit
    # 2.2. Extract results for cross-validation, training, and test datasets
    cross_validation[[classif_name]] <- extract_cross_validation_results(resample, expression_train, target_var)
    train[[classif_name]] <- extract_train_test_results(fit, expression_train, target_var)
    test[[classif_name]] <- extract_train_test_results(fit, expression_test, target_var)
  }
  # 3. Save in the final list with your original structure
  extracted_results$cross_validation <- cross_validation
  extracted_results$train <- train
  extracted_results$test <- test
  extracted_results$times <- times

  return(extracted_results)
}