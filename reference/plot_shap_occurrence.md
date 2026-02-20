# Plot SHAP Magnitude vs Occurrence

This function creates a scatter plot of SHAP magnitude (aggregated_max)
against occurrence count for each gene, colored by gene status. An
optional horizontal line can be added to indicate a threshold for SHAP
magnitude. Marginal density plots are included on the axes to show the
distribution of occurrence counts and SHAP magnitudes.

## Usage

``` r
plot_shap_occurrence(
  data,
  y_threshold = NULL,
  title = "SHAP Magnitude vs Occurrence"
)
```

## Arguments

- data:

  Dataframe containing columns 'occurrence_count', 'aggregated_max', and
  'status'.

- y_threshold:

  Numeric. Optional threshold for SHAP magnitude to be indicated with a
  horizontal line.

- title:

  Character. Title for the plot.

## Value

A ggplot object with marginal density plots.
