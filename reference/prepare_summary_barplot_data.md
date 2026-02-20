# Prepare data for Summary Barplot

Filters the performance data by a specific metric and calculates mean,
SD, and median per classifier.

## Usage

``` r
prepare_summary_barplot_data(data, target_metric, value_col)
```

## Arguments

- data:

  Dataframe. Processed metrics (long format).

- target_metric:

  Character. The metric to summarize (e.g., "normMCC").

- value_col:

  Character. The name of the column containing the values.

## Value

A summarized dataframe ready for plotting.
