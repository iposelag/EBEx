# Create Enrichment Logical Matrix

This function takes the output from `run_phenotype_enrichment` and
creates a logical matrix indicating which genes are significantly
associated with each phenotypic variable based on adjusted p-values. The
resulting matrix has genes as rows and phenotypic variables as columns,
with TRUE indicating significant association and FALSE otherwise.

## Usage

``` r
get_enrichment_logical_matrix(enrichment_df, alpha = 0.05)
```

## Arguments

- enrichment_df:

  Output from `run_phenotype_enrichment`.

- alpha:

  Numeric. Significance threshold.

## Value

A logical matrix indicating significant gene-phenotype associations.
