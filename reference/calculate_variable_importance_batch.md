# Calculate Variable Importance for a Batch of Samples

This function calculates SHAP values for a batch of samples using a
DALEX explainer. It iterates over the specified number of samples,
generates SHAP contributions for each, and aggregates the results into a
summary dataframe containing mean SHAP values for each variable and
sample. The function also measures the time taken for the entire
process.

## Usage

``` r
calculate_variable_importance_batch(
  explainer,
  train_data,
  target_var,
  n_samples,
  directory_to_save = NULL
)
```

## Arguments

- explainer:

  DALEX explainer object.

- train_data:

  Dataframe used for training.

- target_var:

  Character. Name of the target variable column.

- n_samples:

  Integer. Number of samples to process.

- directory_to_save:

  Character. Optional directory to save individual results.

## Value

A list containing the summary dataframe of mean SHAP values and the time
taken for the process.
