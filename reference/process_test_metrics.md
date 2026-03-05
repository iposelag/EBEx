# Process Test Metrics from Results List:

converts the nested list of test metrics into a long-format data frame

## Usage

``` r
process_test_metrics(results_models, classifiers)
```

## Arguments

- results_models:

  The list containing the results of the ML models, including test
  metrics.

- classifiers:

  A character vector of classifier names to process (e.g., c("rf",
  "svm_r", "knn")).

## Value

A data frame with columns: classifier, metric, estimate.
