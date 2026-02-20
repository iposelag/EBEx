# Post-DEA Analysis and Diagnostics

This function takes the fitted model object from the DEA analysis and
performs post-analysis diagnostics. It generates several plots to
visualize the results, including a Venn diagram of DEGs, a histogram of
p-values, and a mean-difference (MD) plot. The function also returns a
matrix indicating which genes are differentially expressed based on the
specified p-value and log-fold change thresholds. The generated plots
are saved in a "DEA_plots" subdirectory within the specified output
directory.

## Usage

``` r
dea_analysis(fit_object, pv = 0.05, lfc = 0, coef = 1, output_dir)
```

## Arguments

- fit_object:

  MArrayLM object. Output from `run_dea_limma`.

- pv:

  Numeric. P-value threshold.

- lfc:

  Numeric. Log-fold change threshold.

- coef:

  Integer. Coefficient to analyze.

- output_dir:

  Character. Directory to save plots.

## Value

The results Matrix with DEGs.
