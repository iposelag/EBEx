# Plot GTEx Expression Heatmap

This function generates a ComplexHeatmap based on median TPM values with
GTEx-style scaling. It includes annotations for tissue macro-categories
and uses a color scheme inspired by GTEx visualizations.

## Usage

``` r
plot_gtex_heatmap(
  tpm_matrix,
  tissue_macro_map,
  tissue_fine_map,
  title = "GTEx Expression"
)
```

## Arguments

- tpm_matrix:

  Matrix. Genes as rows, GTEx tissue IDs as columns.

- tissue_macro_map:

  Named vector mapping tissue IDs to macro-categories.

- tissue_fine_map:

  Named vector mapping tissue IDs to pretty names.

- title:

  Character. Title for the heatmap.

- output_path:

  Character. Path to save the PDF.

## Value

A ComplexHeatmap object.
