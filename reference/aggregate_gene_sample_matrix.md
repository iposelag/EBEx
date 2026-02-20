# Generate Gene-Sample Matrix

This function takes a list of normalized SHAP dataframes (one per model)
and creates a wide matrix where rows represent genes and columns
represent samples. The values in the matrix are the mean normalized SHAP
values for each gene-sample combination, averaged across all models. The
'absolute' parameter allows the user to choose whether to use the
absolute values of the normalized SHAP values or to keep their original
sign, depending on whether they want to focus on the magnitude of
importance or also consider the direction of the effect.

## Usage

``` r
aggregate_gene_sample_matrix(normalized_shap_list, absolute)
```

## Arguments

- normalized_shap_list:

  A named list of normalized SHAP dataframes.

- absolute:

  Logical. If TRUE, use absolute values of normalized SHAP; if FALSE,
  use original values.

## Value

A wide matrix of Gene x Sample.
