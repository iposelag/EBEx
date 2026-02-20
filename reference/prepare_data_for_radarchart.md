# Prepare data for Radar Chart

Filters metrics and constructs the specific matrix format required by
fmsb, identifying the best performing input list for each classifier.

## Usage

``` r
prepare_data_for_radarchart(results, metric, estimate_column, min = 0, max = 1)
```

## Arguments

- results:

  Dataframe containing processed metrics.

- metric:

  Character. The metric to plot (e.g., "normMCC").

- estimate_column:

  Character. The column name containing the values (e.g., "estimate" or
  "mean").

- min:

  Numeric. Minimum value for the radar scale.

- max:

  Numeric. Maximum value for the radar scale.

## Value

A list containing the radarchart dataframe, best selection names, and
best values.
