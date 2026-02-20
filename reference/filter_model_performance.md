# Filter best scenarios (input list + classifier) based on performance metrics

This function identifies which combinations of input gene lists and
classifiers meet a specified performance threshold based on both
cross-validation and test metrics.

## Usage

``` r
filter_model_performance(
  results_models_list,
  metric = "normMCC",
  threshold = 0.75
)
```

## Arguments

- results_models_list:

  A named list of performance objects.

- metric:

  Character. Metric to filter by (e.g., "normMCC").

- threshold:

  Numeric. Minimum value for inclusion.

## Value

A data frame of valid (input_list, classifier) pairs.
