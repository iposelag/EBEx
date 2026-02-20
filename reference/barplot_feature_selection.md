# Plot Barplot of input lists

This function generates a barplot showing the number of genes in each
input list used for feature selection. It uses the package-wide color
scheme for consistency.

## Usage

``` r
barplot_feature_selection(feature_selection_data, output_file = NULL)
```

## Arguments

- feature_selection_data:

  Dataframe with columns 'input_list' and 'intersection_count'.

- output_file:

  Character. Path to save the output PDF.

## Value

Invisibly returns the plot object.
