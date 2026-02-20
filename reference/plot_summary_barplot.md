# Plot Summary Barplot of Classifier Performance

Generates a barplot with error bars (SD), median points, and SD
annotations.

## Usage

``` r
plot_summary_barplot(
  summary_data,
  metric_name,
  colors = ml_models_colors,
  base_size = 14
)
```

## Arguments

- summary_data:

  Dataframe. Output from `prepare_summary_barplot_data`.

- metric_name:

  Character. Metric name for titles and axes.

- colors:

  Named vector. Colors for each classifier.

- base_size:

  Numeric. Base font size (default 14).

## Value

A ggplot object.
