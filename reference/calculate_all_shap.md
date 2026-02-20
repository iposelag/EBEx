# Calculate SHAP values for all ML models

This function iterates over a list of machine learning models,
calculates SHAP values for each model using the
`calculate_variable_importance_batch` function, and compiles the results
into a list. It also ensures that the necessary package for explaining
tidymodels is available and handles the renaming of model labels for
better readability in the output.

## Usage

``` r
calculate_all_shap(
  ml_models_to_run,
  results_models,
  expression_data,
  target_var,
  directory_to_save = NULL
)
```

## Arguments

- ml_models_to_run:

  Character vector of model names.

- results_models:

  List of fitted models (output from ML_models).

- expression_data:

  Dataframe. The data used for training.

- target_var:

  Character. The target variable column name.

- directory_to_save:

  Character. Optional directory path.

## Value

A list containing SHAP values for each model and the time taken for the
entire process.
