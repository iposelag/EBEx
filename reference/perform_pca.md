# Perform PCA

Computes a Principal Component Analysis (PCA) on the provided data
matrix (samples as rows, features as columns), scaling and centering it
beforehand. Useful for inspecting global structure and batch effects.
Optionally saves a scree plot (eigenvalues/variance explained per
component) to disk.

## Usage

``` r
perform_pca(
  data,
  num_components = 30,
  save_eig_plot = FALSE,
  eig_ncp = 15,
  output_dir = NULL,
  plot_filename = "pca_eigenvalues.png"
)
```

## Arguments

- data:

  Matrix or dataframe. Samples as rows, features (e.g. genes) as
  columns.

- num_components:

  Integer. Number of principal components to keep.

- save_eig_plot:

  Logical. If `TRUE`, saves a scree plot of the first `eig_ncp`
  components via
  [`factoextra::fviz_eig()`](https://rpkgs.datanovia.com/factoextra/reference/eigenvalue.html).
  Default `FALSE`.

- eig_ncp:

  Integer. Number of components to show in the scree plot (passed to
  `factoextra::fviz_eig(ncp = eig_ncp)`). Default 15.

- output_dir:

  Character. Directory where the scree plot is saved when
  `save_eig_plot = TRUE`. Created if it doesn't exist. Required in that
  case.

- plot_filename:

  Character. File name for the saved scree plot. Default
  `"pca_eigenvalues.png"`.

## Value

A list containing `pca_data` (a dataframe with the retained principal
components, samples as row names) and `pca_summary` (the
`summary.prcomp` object, useful for inspecting variance explained per
component).
