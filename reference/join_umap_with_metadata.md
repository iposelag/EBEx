# Join UMAP Coordinates with Sample Metadata

Extracts the UMAP layout (UMAP1/UMAP2) from a `umap` object (as returned
by
[`perform_umap()`](https://iposelag.github.io/EBEx/reference/perform_umap.md))
and joins it with sample-level phenotypic/metadata, matching by
`sample_id`.

## Usage

``` r
join_umap_with_metadata(umap_data, metadata)
```

## Arguments

- umap_data:

  List. A `umap` object, as returned by
  [`perform_umap()`](https://iposelag.github.io/EBEx/reference/perform_umap.md).

- metadata:

  Dataframe. Sample metadata, must contain a `sample_id` column.

## Value

A dataframe with the UMAP1/UMAP2 coordinates joined to the metadata.
