# Full Aggregation Pipeline (Example)

This demonstrates how to use the individual steps to go from raw results
to a candidate gene list.

## Usage

``` r
run_aggregation_pipeline(
  results_models_list,
  results_shap_list,
  threshold = 0.75,
  absolute = TRUE,
  candidate_results_dir = NULL
)
```

## Arguments

- results_models_list:

  A named list of performance objects.

- results_shap_list:

  A named list of SHAP dataframes corresponding to the models.

- threshold:

  Numeric. Performance threshold for filtering models.

- absolute:

  Logical. Whether to use absolute values for the gene-sample matrix.
