# Join PCA Coordinates with Sample Metadata

Joins the PCA coordinates (as returned in `perform_pca()$pca_data`) with
sample-level phenotypic/metadata, matching by `sample_id`.

## Usage

``` r
join_pca_with_metadata(pca_data, metadata)
```

## Arguments

- pca_data:

  Dataframe. PCA coordinates with samples as row names, as returned in
  `perform_pca()$pca_data`.

- metadata:

  Dataframe. Sample metadata, must contain a `sample_id` column.

## Value

A dataframe with the PCA coordinates joined to the metadata.
