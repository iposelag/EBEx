# Extract models results (cross-validation, train, test and timing)

This function extracts the results from fitted machine learning models,
including cross-validation metrics, training and test metrics, and
timing information. It takes a list of fitted models and the
corresponding training and test datasets, and returns a structured list
containing the extracted results for each model. The function relies on
internal helper functions to extract the relevant metrics and confusion
matrices for both training and test datasets, as well as the
cross-validation results.

## Usage

``` r
extract_models_results(
  models_to_run,
  results_models,
  expression_train,
  expression_test,
  target_var
)
```

## Arguments

- models_to_run:

  Vector of models names.

- results_models:

  List of fitted models.

- expression_train:

  Training data.

- expression_test:

  Test data.

- target_var:

  Name of the target variable.

## Value

A list containing the extracted results for each model, including
cross-validation metrics, training and test metrics, confusion matrices,
misclassified samples, and timing information.
