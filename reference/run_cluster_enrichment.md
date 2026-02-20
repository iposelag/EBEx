# Enrichment Analysis of Clusters vs Phenotypes

Performs "Cluster vs Rest" statistical tests for each phenotypic
variable. For numeric variables, it tests for normality and applies
either a t-test or Wilcoxon test. For categorical variables, it applies
a chi-squared test. The function returns a dataframe summarizing the
results, including adjusted p-values for multiple testing. Clusters with
fewer than 3 samples are automatically skipped to ensure statistical
validity. The results include the type of test used, p-values, test
statistics, and any enriched classes for categorical variables.

## Usage

``` r
run_cluster_enrichment(cluster_vec, phenotype_data)
```

## Arguments

- cluster_vec:

  Named integer vector. Output from
  `calculate_sample_clustering$clusters`.

- phenotype_data:

  Dataframe. Phenotypic variables (samples as rows).

## Value

A dataframe summarizing the enrichment results for each cluster and
phenotypic variable.
