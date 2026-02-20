# GLM Multiclass Model

This function trains a Generalized Linear Model (GLM) for multiclass
classification using the `parsnip` and `workflows` packages. It performs
hyperparameter tuning using Bayesian optimization and evaluates the
model using resampling. The fitted model and resampling results can be
saved as bundles in a specified directory.

## Usage

``` r
glm_multiclass_model(
  data,
  target_var,
  folds,
  metrics,
  control,
  directory_to_save = NULL
)
```

## Arguments

- data:

  Data frame containing the training data.

- target_var:

  Character string of the target variable column name.

- folds:

  Resampling object (e.g., vfold_cv) for cross-validation

- metrics:

  Metrics object for model evaluation.

- control:

  Control object for tuning and resampling.

- directory_to_save:

  Character string specifying the directory to save the model bundles.
  Default is NULL (do not save).

## Value

A list containing the time taken for training, resampling results, and
the fitted model object.
