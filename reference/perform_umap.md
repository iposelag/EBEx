# Perform UMAP

Computes a UMAP embedding on the provided data matrix (samples as rows,
features as columns), scaling it beforehand.

## Usage

``` r
perform_umap(data)
```

## Arguments

- data:

  Matrix or dataframe. Samples as rows, features (e.g. genes) as
  columns.

## Value

A `umap` object, as returned by `umap::umap()`.
