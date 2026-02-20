# Plot Functional Enrichment Dotplot

This function creates a dotplot for functional enrichment results,
showing gene ratio, count, and adjusted p-value for the top terms.

## Usage

``` r
plot_enrichment_dotplot(
  enrichment_df,
  title,
  n_input_genes,
  top_n = 15,
  wrap_width = 30,
  color_limits = c(1, 5),
  size_limits = c(1, 30),
  aspect_ratio = 1,
  output_file = NULL
)
```

## Arguments

- enrichment_df:

  Dataframe from one database of Enrichr results.

- title:

  Character. Title for the plot.

- n_input_genes:

  Integer. Total number of genes used in the enrichment call.

- top_n:

  Integer. Number of top terms to show.

- wrap_width:

  Integer. Width for wrapping long term names.

- color_limits:

  Numeric vector of length 2. Limits for the color scale.

- size_limits:

  Numeric vector of length 2. Limits for the size scale.

- aspect_ratio:

  Numeric. Aspect ratio for the plot.

- output_file:

  Character. Optional path to save the plot (e.g.,
  "enrichment_dotplot.pdf"). If NULL (default), the plot is displayed in
  the active device.

## Value

A ggplot object representing the enrichment dotplot.
