# Run ML Pipeline

This function serves as a wrapper to run multiple machine learning
models on a given dataset. It takes the training data, target variable,
and a list of models to run, and executes the corresponding model
training functions. The results from each model, including training
time, resampling results, and fitted model objects, are collected in a
list and returned. Optionally, the trained models and resampling results
can be saved to a specified directory.

## Usage

``` r
ML_models(
  data,
  target_var,
  models_to_run = c("rf", "svm_r", "svm_p", "glm", "knn", "xgb"),
  directory_to_save = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame containing the training data.

- target_var:

  Name of the target variable.

- models_to_run:

  Vector of model names to run.

- directory_to_save:

  Directory to save the trained models and results (optional).

## Value

A list containing the results from each model, including training time,
resampling results, and fitted model objects.
