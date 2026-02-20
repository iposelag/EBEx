# Plot DEA Candidate Genes

Generates either an MA plot or a Volcano plot for DEA candidate genes,
with points colored by gene status and threshold lines for log fold
change and adjusted p-value.

## Usage

``` r
plot_dea_candidates(
  deg_data,
  type = "MA",
  lfc_thresh = log2(1.5),
  p_val_thresh = 0.01
)
```

## Arguments

- deg_data:

  Dataframe containing DEA results with columns 'AveExpr', 'logFC',
  'adj.P.Val', and 'status'.

- type:

  Character. Type of plot to generate: "MA" or "Vol

- lfc_thresh:

  Numeric. Log fold change threshold for significance (default
  log2(1.5)).

- p_val_thresh:

  Numeric. Adjusted p-value threshold for significance (default 0.01).

## Value

A ggplot object representing the MA or Volcano plot.
