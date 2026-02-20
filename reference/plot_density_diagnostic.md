# Plot Density Diagnostic for Threshold Selection

This function takes the output from `compute_density_analysis` and
generates two plots: a density curve of the specified column and a line
plot of the slope differences across windows

## Usage

``` r
plot_density_diagnostic(
  analysis_list,
  suggested_threshold = NULL,
  output_file = NULL
)
```

## Arguments

- analysis_list:

  List. The output from `compute_density_analysis`.

- suggested_threshold:

  Numeric. Optional value to draw a vertical line.

- output_file:

  Character. Optional path to save the plot (e.g., "plot.pdf"). If NULL
  (default), the plot is displayed in the active device.

## Value

A grob object (invisibly) containing the arranged plots.
