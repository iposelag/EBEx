# Compute Density and Derivatives for Ranking

This function calculates the density of aggregated SHAP values and
computes slopes and window-based differences to help identify a cut-off
(Knee/Elbow method).

## Usage

``` r
compute_density_analysis(data, column = "aggregated_max", window_size = 0.001)
```

## Arguments

- data:

  Dataframe. The gene summary output from aggregation (must contain the
  score column).

- column:

  Character. The column name to analyze (default "aggregated_max").

- window_size:

  Numeric. The continuous window size for difference calculation.

## Value

A list containing 'density_df' (x, y, slope), 'window_results', and the
original 'data'.
