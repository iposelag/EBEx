# Create Enrichment Matrix padjusted values

This function takes the output from `run_phenotype_enrichment` and
creates a matrix of adjusted p-values for each gene-phenotype
association. The resulting matrix has genes as rows and phenotypic
variables as columns, with the values representing the adjusted p-values
for each association.

## Usage

``` r
get_enrichment_matrix(enrichment_df, alpha = 0.05)
```

## Arguments

- enrichment_df:

  Output from `run_phenotype_enrichment`.

- alpha:

  Numeric. Significance threshold.

## Value

A matrix of adjusted p-values for gene-phenotype associations.
