# Run Mclust Analysis for Gene Ranking

This function fits Gaussian Mixture Models to the aggregated SHAP values
to identify underlying groups of genes.

## Usage

``` r
run_mclust_analysis(data, column = "aggregated_max", G = 1:9)
```

## Arguments

- data:

  Dataframe. The gene summary output from aggregation.

- column:

  Character. The column name to analyze (default "aggregated_max").

- G:

  Numeric vector. The number of clusters to test (default 1:9).

## Value

A list containing the Mclust model, BIC, ICL, used data, and column
name.
