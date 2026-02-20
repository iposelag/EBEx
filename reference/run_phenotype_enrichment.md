# Run Full Phenotype Enrichment Analysis

This function takes a SHAP matrix (samples x genes) and a phenotype
dataframe, performs the `perform_pheno_test` for each gene and specified
phenotypic variable, and compiles the results into a tidy dataframe. It
also adjusts p-values for multiple testing within each variable and
filters out any residual terms or NAs before returning the final
results. The output includes the gene, variable, test type, p-values,
and adjusted p-values for each gene-phenotype association tested.

## Usage

``` r
run_phenotype_enrichment(shap_matrix, phenotype_data, genes, variables)
```

## Arguments

- shap_matrix:

  Matrix or Dataframe. Samples as rows, genes as columns.

- phenotype_data:

  Dataframe. Phenotypic variables for the same samples.

- genes:

  Character vector. List of genes to test.

- variables:

  Character vector. List of phenotypic variables to test.

## Value

A tidy dataframe summarizing the results of the phenotype enrichment
analysis, including p-values and adjusted p-values for each
gene-phenotype association.
