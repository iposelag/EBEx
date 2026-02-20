# Generate Gene-Level Aggregation

This function takes a list of normalized SHAP dataframes (one per model)
and aggregates the importance of each gene across all models. It
calculates several metrics for each gene, including the maximum, mean,
median, standard deviation of the mean normalized SHAP values, and the
count of how many models included that gene. The resulting summary
dataframe is sorted by the maximum aggregated importance, allowing for
easy identification of top candidate genes.

## Usage

``` r
aggregate_genes(normalized_shap_list)
```

## Arguments

- normalized_shap_list:

  A named list of normalized SHAP dataframes.

## Value

A summary dataframe with metrics per gene.
