# Plot Mclust Diagnostics

Visualizes BIC, ICL and the resulting clusters relative to the SHAP
scores.

## Usage

``` r
plot_mclust_diagnostic(mclust_list, sel_groups = NULL, output_file = NULL)
```

## Arguments

- mclust_list:

  List. Output from `run_mclust_analysis`.

- sel_groups:

  Numeric vector. Optional: clusters to highlight with a threshold line.

- output_file:

  Character. Optional: path to save as PDF.

## Value

NULL. Generates plots for diagnostics.
