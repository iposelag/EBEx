# Plot Computation Time vs Gene List Size

This function creates a plot of computation time against the number of
genes in the input list for each classifier.

## Usage

``` r
plot_time_vs_genes(time_df, gene_list_counts, output_dir, colors = NULL)
```

## Arguments

- time_df:

  Dataframe containing columns 'time', 'classifier', and 'input_list'.

- gene_list_counts:

  Named vector with the number of genes for each input list.

- output_dir:

  Character. Directory to save the plot.

- colors:

  Named vector. Optional colors for each classifier.

## Value

Invisibly returns the plot object.
